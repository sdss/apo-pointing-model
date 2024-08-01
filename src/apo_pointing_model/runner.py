#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-07-31
# @Filename: runner.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import asyncio
import pathlib

import polars
from clu.legacy import TronConnection
from clu.legacy.types.parser import Reply
from pydantic import BaseModel, Field

from sdsstools import get_sjd

from apo_pointing_model import log
from apo_pointing_model.sample import get_random_sample, to_icrs
from apo_pointing_model.tcc import TCCHelper


SCHEMA: dict[str, polars.DataTypeClass] = {
    "alt": polars.Float64,
    "az": polars.Float64,
    "rot": polars.Float64,
    "ra": polars.Float64,
    "dec": polars.Float64,
    "done": polars.Boolean,
    "failed": polars.Boolean,
    "mjd": polars.Int32,
    "gimg": polars.Int32,
    "n_cameras_solved": polars.Int32,
    "ra_bore": polars.Float64,
    "dec_bore": polars.Float64,
    "offset_ra": polars.Float64,
    "offset_dec": polars.Float64,
    "offset_rot": polars.Float64,
    "tai_ref": polars.Float64,
    "ptcorr_azcorr": polars.Float64,
    "ptcorr_altcorr": polars.Float64,
    "ptcorr_xpos": polars.Float64,
    "ptcorr_ypos": polars.Float64,
    "ptdata_azphys": polars.Float64,
    "ptdata_altphys": polars.Float64,
    "ptdata_azmount": polars.Float64,
    "ptdata_altmount": polars.Float64,
    "ptdata_rotphys": polars.Float64,
}


class PointingData(BaseModel):
    """Pointing model data."""

    alt: float
    az: float
    rot: float
    ra: float | None = None
    dec: float | None = None
    done: bool = False
    failed: bool = False
    mjd: float = Field(default=get_sjd("APO"))
    gimg: int | None = None
    n_cameras_solved: int = 0
    tai_ref: float | None = None
    ra_bore: float | None = None
    dec_bore: float | None = None
    offset_ra: float | None = None
    offset_dec: float | None = None
    offset_rot: float | None = None
    ptcorr_azcorr: float | None = None
    ptcorr_altcorr: float | None = None
    ptcorr_xpos: float | None = None
    ptcorr_ypos: float | None = None
    ptdata_azphys: float | None = None
    ptdata_altphys: float | None = None
    ptdata_azmount: float | None = None
    ptdata_altmount: float | None = None
    ptdata_rotphys: float | None = None


async def get_pointing_data(
    npoints: int | None,
    output_file: str | pathlib.Path,
    overwrite: bool = False,
    write_csv: bool = True,
    alt_range: tuple[float, float] = (28, 85),
    az_range: tuple[float, float] = (0, 359),
    write_log: bool = True,
):
    """Collects pointing model data."""

    output_file = pathlib.Path(output_file)
    if write_log:
        log.start_file_logger(str(output_file.with_suffix(".log")), rotating=False)

    if output_file.exists() and not overwrite:
        log.warning(f"Found file {output_file!s}. NOT generating a new grid.")
        data_df = polars.read_parquet(output_file)
        data = [PointingData(**row) for row in data_df.to_dicts()]
    else:
        if npoints is None:
            raise ValueError("npoints is required if output_file does not exist.")

        log.info("Generating new pointing model grid.")
        altaz = get_random_sample(npoints, alt_range=alt_range, az_range=az_range)
        data: list[PointingData] = []
        for alt, az in altaz:
            data.append(PointingData(alt=alt, az=az, rot=0.0))

        write_data(data, output_file, write_csv=write_csv)

    log.info("Creating connection to Tron and waiting for keys.")
    tron = TronConnection(
        "APO.pointing_model",
        "sdss5-hub.apo.nmsu.edu",
        models=["tcc", "mcp"],
    )
    await tron.start()
    await asyncio.sleep(5)

    # Copy of the HAL TCC helper for slewing.
    tcc = TCCHelper(tron)

    for irow, pdata in enumerate(data):
        if pdata.done or pdata.failed:
            log.warning(f"Skipping row {irow!r} as it is already done or failed.")
            continue

        alt = pdata.alt
        az = pdata.az

        ra, dec = to_icrs(alt, az)
        pdata.ra = ra
        pdata.dec = dec

        log.warning(f"Slewing to target #{irow+1}: RA={ra:.4f}, Dec={dec:.4f}.")

        if not tcc.axes_are_clear():
            raise RuntimeError("Axes are not clear. Aborting.")

        track_command = f"track {ra},{dec} icrs /rottype=obj /rotang=0"
        await tcc.do_slew(  # This raises internally.
            track_command=track_command,
            callback=log_reply,
        )

        # Make sure the guider knows the position of the field.
        log.debug("Setting a fake jaeger field.")
        await tron.send_command(
            "jaeger",
            f"configuration fake-field {ra} {dec} 0.0",
        )

        log.info("Astrometrically solving the field with cherno.")
        exp_time: float = 5

        while True:
            cmd_time = await tron.send_command("tcc", "show time")
            tai0 = cmd_time.replies.get("tai")[0]

            cmd_acq = await tron.send_command(
                "cherno",
                f"acquire -t {exp_time} --mode astrometrynet --no-continuous --no-apply",
                callback=log_reply,
            )

            did_fail = cmd_acq.status.did_fail
            acquisition_valid = cmd_acq.replies.get("acquisition_valid")[0]

            try:
                astrometry_fit = cmd_acq.replies.get("astrometry_fit")
                if float(astrometry_fit[2]) == -999.0:
                    log.warning(
                        "Field was not solved with fit_SP or failed to acquire."
                    )
                    astrometry_fit = None
            except KeyError:
                log.warning("The astrometry_fit keyword was not output.")
                astrometry_fit = None

            if did_fail or astrometry_fit is None or not acquisition_valid:
                exp_time += 10
                if exp_time <= 25:
                    log.error(
                        "Failed to solve field. Retrying "
                        f"with exp_time={exp_time} seconds."
                    )
                    continue

                log.error("Failed to solve field. Skipping this pointing.")
                pdata.failed = True
                break

            pdata.gimg = int(astrometry_fit[0])
            pdata.n_cameras_solved = int(astrometry_fit[1])
            pdata.ra_bore = float(astrometry_fit[2])
            pdata.dec_bore = float(astrometry_fit[3])
            pdata.offset_ra = float(astrometry_fit[7])
            pdata.offset_dec = float(astrometry_fit[8])
            pdata.offset_rot = float(astrometry_fit[9])
            pdata.tai_ref = float(tai0 + exp_time / 2)

            log.info("Retrieving pdata.")

            ptcorr_string = f"ptcorr {pdata.ra_bore}, {pdata.dec_bore} icrs 0,0,0,0,{pdata.tai_ref} instrument"
            log.debug(f"Running command tcc {ptcorr_string!r}.")

            pdata_cmd = await tron.send_command(
                "tcc",
                ptcorr_string,
                callback=log_reply,
            )

            ptcorr = pdata_cmd.replies.get("ptcorr")
            ptdata = pdata_cmd.replies.get("ptdata")
            pdata.ptcorr_azcorr = ptcorr[0]
            pdata.ptcorr_altcorr = ptcorr[1]
            pdata.ptcorr_xpos = ptcorr[2]
            pdata.ptcorr_ypos = ptcorr[3]
            pdata.ptdata_azphys = ptdata[0]
            pdata.ptdata_altphys = ptdata[1]
            pdata.ptdata_azmount = ptdata[2]
            pdata.ptdata_altmount = ptdata[3]
            pdata.ptdata_rotphys = ptdata[4]

            pdata.done = True

            log.info(f"Field solved. Data: {pdata.model_dump_json(indent=2)}")

            break

        write_data(data, output_file, write_csv=write_csv)


def write_data(
    data: list[PointingData],
    output_file: pathlib.Path,
    write_csv: bool = True,
):
    """Writes the data to disk.

    Parameters
    ----------
    data
        The list of pointing data.
    output_file
        The output file. Expected to be the path to a Parquet file.
    write_csv
        If True, writes the data to a CSV file as well. The file will have the same
        name as ``output_file`` but with the extension replaced to ``.csv``.

    """

    data_df = polars.DataFrame([d.model_dump() for d in data], schema=SCHEMA)
    data_df.write_parquet(output_file)

    if write_csv:
        csv_file = output_file.with_suffix(".csv")
        data_df.write_csv(csv_file)


def log_reply(reply: Reply):
    """Logs a reply to a command."""

    log.debug(reply.string)
