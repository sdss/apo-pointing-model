#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-07-31
# @Filename: tcc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2021-09-26
# @Filename: tcc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

from typing import Callable, cast

from clu.command import Command, CommandStatus
from clu.legacy.tron import TronConnection

from apo_pointing_model import log


__all__ = ["TCCHelper"]


class HALHelper:
    """A helper class to control an actor or piece of hardware."""

    name: str | None = None

    def __init__(self, tron: TronConnection):
        self.tron = tron

    async def _send_command(
        self,
        target: str,
        cmd_str: str,
        raise_on_fail: bool = True,
        **kwargs,
    ):
        """Sends a command to a target."""

        if self.tron is None or self.tron.connected() is False:
            raise RuntimeError("Not connected to Tron. Cannot send commands.")

        cmd = await self.tron.send_command(target, cmd_str, **kwargs)

        if raise_on_fail and cmd.status.did_fail:
            if cmd.status == CommandStatus.TIMEDOUT:
                raise RuntimeError(f"Command '{target} {cmd_str}' timed out.")
            else:
                raise RuntimeError(f"Command '{target} {cmd_str}' failed.")

        return cast(Command, cmd)


class TCCHelper(HALHelper):
    """Helper for the TCC."""

    name = "tcc"

    is_slewing: bool = False

    async def axis_init(self) -> bool:
        """Executes TCC axis init or fails."""

        status = await self._send_command(
            "tcc",
            "axis status",
            time_limit=20.0,
            raise_on_fail=False,
        )
        if status.status.did_fail:
            raise RuntimeError("'tcc status' failed. Is the TCC connected?")

        if self.check_stop_in() is True:
            raise RuntimeError(
                "Cannot tcc axis init because of bad axis status: "
                "Check stop buttons on Interlocks panel."
            )

        sem = await self.mcp_semaphore_ok()
        if sem is False:
            raise RuntimeError("Failed getting the semaphore information.")

        if sem == "TCC:0:0" and self.axes_are_clear():
            log.debug(
                "Axes clear and TCC has semaphore. "
                "No axis init needed, so none sent."
            )
            return True

        log.debug("Sending tcc axis init.")
        axis_init_cmd_str = "axis init"
        if self.below_alt_limit():
            log.warning(
                "Altitude below interlock limit! Only initializing altitude "
                "and rotator: cannot move in az."
            )
            axis_init_cmd_str += " rot,alt"
        axis_init_cmd = await self._send_command(
            "tcc",
            axis_init_cmd_str,
            time_limit=20.0,
            raise_on_fail=False,
        )

        if axis_init_cmd.status.did_fail:
            log.error("Cannot slew telescope: failed tcc axis init.")
            raise RuntimeError("Cannot slew telescope: check and clear interlocks?")

        return True

    async def axis_stop(self, command: Command, axis: str = "") -> bool:
        """Issues an axis stop to the TCC."""

        axis_stop_cmd = await self._send_command(
            "tcc",
            f"axis stop {axis}".strip(),
            time_limit=30.0,
            raise_on_fail=False,
        )

        if axis_stop_cmd.status.did_fail:
            raise RuntimeError(
                "Error: failed to cleanly stop telescope via tcc axis stop."
            )

        self.is_slewing = False
        return True

    def get_bad_axis_bits(self, axes=("az", "alt", "rot"), mask=None):
        """Return the bad status bits for the requested axes."""

        if mask is None:
            mask = self.tron.models["tcc"]["axisBadStatusMask"][0]

        if mask is None:
            log.error("Cannot retrieve TCC axisBadStatusMask.")
            return [1, 1, 1]

        return [(self.tron.models["tcc"][f"{axis}Stat"][3] & mask) for axis in axes]

    def check_stop_in(self, axes=("az", "alt", "rot")) -> bool:
        """Returns `True` if any stop bit is set in the ``<axis>Stat`` TCC keywords.

        The [az,alt,rot]Stat[3] bits show the exact status:
        http://www.apo.nmsu.edu/Telescopes/HardwareControllers/AxisControllers.html#25mStatusBits

        """
        try:
            # 0x2000 is "stop button in"
            return any(self.get_bad_axis_bits(axes=axes, mask=0x2000))
        except TypeError:
            # Some axisStat is unknown (and thus None)
            return False

    async def mcp_semaphore_ok(self):
        """Returns the semaphore if the semaphore is owned by the TCC or nobody."""

        mcp_model = self.tron.models["mcp"]

        sem = mcp_model["semaphoreOwner"]
        if sem is None:
            sem_show_cmd = await self._send_command(
                "mcp",
                "sem.show",
                time_limit=20.0,
                raise_on_fail=False,
            )

            if sem_show_cmd.status.did_fail:
                raise RuntimeError("Cannot get mcp semaphore. Is the MCP alive?")

            sem = mcp_model["semaphoreOwner"]

        if (
            (sem[0] != "TCC:0:0")
            and (sem[0] != "")
            and (sem[0] != "None")
            and (sem[0] is not None)
        ):
            raise RuntimeError(
                f"Cannot axis init: Semaphore is owned by  {sem[0]}. "
                "If you are the owner (e.g., via MCP Menu), release it and try again. "
                "If you are not the owner, confirm that you can steal "
                "it from them, then issue: mcp sem.steal"
            )

        return sem

    def axes_are_clear(self, axes=("az", "alt", "rot")) -> bool:
        """Checks that no bits are set in any axis status field."""

        axes = [ax for ax in axes if ax != "alt" or self.below_alt_limit() is False]
        try:
            tcc_model = self.tron.models["tcc"]
            return all((tcc_model[f"{axis}Stat"][3] == 0) for axis in axes)
        except TypeError:
            # Some axisStat is unknown (and thus None)
            return False

    def below_alt_limit(self) -> bool:
        """Check if we are below the alt=18 limit that prevents init/motion in az."""

        limit = 18
        return self.tron.models["tcc"]["axePos"][1] < limit

    def check_axes_status(self, status: str) -> bool:
        """Returns `True` if all the axes are at ``status``."""

        axes_status = self.tron.models["tcc"]["AxisCmdState"].value
        return all([axis.lower() == status.lower() for axis in axes_status])

    async def do_slew(
        self,
        coords: dict[str, float] | None = None,
        track_command: str | None = None,
        keep_offsets: bool = True,
        offset: bool = False,
        rotwrap: str | None = None,
        callback: Callable | None = None,
    ) -> bool:
        """Correctly handle a slew command, given what parse_args had received.

        Parameters
        ----------
        coords
            The coordinates where to slew. It must be a dictionary with keys
            ``ra, dec, rot`` or ``alt, az, rot``.
        track_command
            A raw TCC ``track`` command (without the ``tcc`` target) to send.
            In this case other arguments like ``keep_offsets`` or ``rotwrap``
            are ignored.
        keep_offsets
            Whether to keep the existing offsets.
        rotwrap
            The type of ``/RotWrap`` to use.
        offset
            If defined, the coordinates will be treated as an offset and not
            absolute positions.

        """

        tcc_model = self.tron.models["tcc"]

        # Check axes.
        result = self.axes_are_clear()
        if not result:
            raise RuntimeError("Some axes are not clear. Cannot continue.")

        # NOTE: TBD: We should limit which offsets are kept.
        keep_args = "/keep=(arc,gcorr,calib,bore)" if keep_offsets else ""
        rotwrap = f"/rotwrap={rotwrap}" if rotwrap else ""

        slew_cmd = None
        if track_command:
            slew_cmd = self._send_command(
                "tcc",
                track_command,
                time_limit=300,
                raise_on_fail=False,
                callback=callback,
            )
        else:
            if coords is None:
                raise RuntimeError("coords needed if track_command is None.")

            if not offset:
                if "ra" in coords and "dec" in coords and "rot" in coords:
                    ra = coords["ra"]
                    dec = coords["dec"]
                    rot = coords["rot"]

                    log.info(
                        "Slewing to (ra, dec, rot) == "
                        f"({ra:.4f}, {dec:.4f}, {rot:g})"
                    )
                    if keep_args:
                        log.warning("keeping all offsets")

                    slew_cmd = self._send_command(
                        "tcc",
                        f"track {ra}, {dec} icrs /rottype=object/rotang={rot:g} "
                        f"{rotwrap} {keep_args}",
                        time_limit=300,
                        raise_on_fail=False,
                        callback=callback,
                    )

                elif "az" in coords and "alt" in coords and "rot" in coords:
                    alt = coords["alt"]
                    az = coords["az"]
                    rot = coords["rot"]

                    log.info(
                        "Slewing to (az, alt, rot) == "
                        f"({az:.4f}, {alt:.4f}, {rot:.4f})"
                    )

                    slew_cmd = self._send_command(
                        "tcc",
                        f"track {az:f}, {alt:f} mount /rottype=mount "
                        f"/rotangle={rot:f} {rotwrap} {keep_args}",
                        time_limit=300,
                        raise_on_fail=False,
                        callback=callback,
                    )

                else:
                    raise RuntimeError("Not enough coordinates information provided.")

            else:
                if "alt" not in coords or "az" not in coords:
                    raise RuntimeError("No alt/az offsets provided.")

                # In arcsec
                alt = coords["alt"] or 0.0
                az = coords["az"] or 0.0
                rot = coords["rot"] or 0.0

                log.info(f"Offseting alt={alt:.3f}, az={az:.3f}")

                slew_cmd = self._send_command(
                    "tcc",
                    f"offset guide {az/3600.:g},{alt/3600.:g},{rot/3600.:g} /computed",
                    time_limit=300,
                    raise_on_fail=False,
                    callback=callback,
                )

        # "tcc track" in the new TCC is only Done successfully when all requested
        # axes are in the "tracking" state. All other conditions mean the command
        # failed, and the appropriate axisCmdState and axisErrCode will be set.
        # However, if an axis becomes bad during the slew, the TCC will try to
        # finish it anyway, so we need to explicitly check for bad bits.

        try:
            self.is_slewing = True
            slew_cmd = await slew_cmd
        finally:
            self.is_slewing = False

        if slew_cmd.status.did_fail:
            str_axis_state = ",".join(tcc_model["axisCmdState"].value)
            str_axis_code = ",".join(tcc_model["axisErrCode"].value)
            log.error(
                f"tcc track command failed with axis states: {str_axis_state} "
                f"and error codes: {str_axis_code}"
            )
            raise RuntimeError("Failed to complete slew: see TCC messages for details.")

        if self.axes_are_clear() is False:
            axis_bits = self.get_bad_axis_bits()
            log.error(
                "TCC slew command ended with some bad bits set: "
                "0x{:x},0x{:x},0x{:x}".format(*axis_bits)
            )
            raise RuntimeError("Failed to complete slew: see TCC messages for details.")

        return True
