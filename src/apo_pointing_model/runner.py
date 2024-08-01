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
    npoints: int,
    output_file: str | pathlib.Path,
    overwrite: bool = False,
):
    """Collects pointing model data."""

    output_file = pathlib.Path(output_file)

    if output_file.exists() and not overwrite:
        log.info(f"Found file {output_file!s}.")
        data_df = polars.read_parquet(output_file)
        data = [PointingData(**row) for row in data_df.to_dicts()]
    else:
        log.info("Generating new pointing model data.")
        altaz = get_random_sample(npoints)
        data: list[PointingData] = []
        for alt, az in altaz:
            data.append(PointingData(alt=alt, az=az, rot=0.0))

        write_to_parquet(data, output_file)

    log.info("Creating connection to Tron and waiting for keys.")
    tron = TronConnection(
        "APO.pointing_model",
        "sdss5-hub.apo.nmsu.edu",
        models=["tcc", "mcp"],
    )
    await tron.start()
    await asyncio.sleep(10)

    tcc = TCCHelper(tron)

    for irow, pdata in enumerate(data):
        if pdata.done:
            log.warning(f"Skipping row {irow!r} as it is already done.")
            continue

        alt = pdata.alt
        az = pdata.az
        ra, dec = to_icrs(alt, az)

        log.info(f"Slewing to target #{irow+1}: RA={ra:.4f}, Dec={dec:.4f}.")

        if not tcc.axes_are_clear():
            raise RuntimeError("Axes are not clear. Aborting.")

        track_command = rf"track {ra},{dec} icrs \rottype=obj \rotang=0"
        # This raises internally.
        await tcc.do_slew(track_command=track_command, callback=log_reply)

        # Make sure the guider knows the position of the field.
        await tron.send_command(
            "jaeger",
            f"configuration fake-field {ra} {dec} 0.0",
        )

        log.debug("Solving field with cherno.")
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
            if not did_fail:
                valid = cmd_acq.replies.get("acquisition_valid")[0]
            else:
                valid = False

            if did_fail or not valid:
                exp_time += 5
                if exp_time <= 20:
                    log.error(
                        "Failed to solve field. Retrying "
                        f"with exp_time={exp_time} seconds."
                    )
                    continue

                log.error("Failed to solve field. Skipping this pointing.")
                pdata.failed = True

                break

            astrometry_fit = cmd_acq.replies.get("astrometry_fit")

            if astrometry_fit[2] == -999:
                log.error("Field was not solved with fit_SP. Skipping this pointing.")
                pdata.failed = True
                break

            pdata.gimg = astrometry_fit[0]
            pdata.ra_bore = astrometry_fit[2]
            pdata.dec_bore = astrometry_fit[3]
            pdata.offset_ra = astrometry_fit[7]
            pdata.offset_dec = astrometry_fit[8]
            pdata.offset_rot = astrometry_fit[9]
            pdata.tai_ref = tai0 + exp_time / 2

            log.debug("Retrieving pdata.")

            pdata_cmd = await tron.send_command(
                "tcc",
                f"ptcorr {pdata.ra_bore}, {pdata.dec_bore} icrs "
                f"0,0,0,0,{pdata.tai_ref} instrument”",
                callback=log_reply,
            )

            ptcorr = pdata_cmd.replies.get("ptcorr")
            pdata = pdata_cmd.replies.get("pdata")
            pdata.ptcorr_azcorr = ptcorr[0]
            pdata.ptcorr_altcorr = ptcorr[1]
            pdata.ptcorr_xpos = ptcorr[2]
            pdata.ptcorr_ypos = ptcorr[3]
            pdata.ptdata_azphys = pdata[0]
            pdata.ptdata_altphys = pdata[1]
            pdata.ptdata_azmount = pdata[2]
            pdata.ptdata_altmount = pdata[3]
            pdata.ptdata_rotphys = pdata[4]

            pdata.done = True

            log.info(f"Field solved. Data: {data[0].model_dump_json(indent=2)}")

            break

    write_to_parquet(data, output_file)


def write_to_parquet(data: list[PointingData], output_file: str | pathlib.Path):
    data_df = polars.DataFrame(
        [d.model_dump() for d in data],
        schema={
            "alt": polars.Float64,
            "az": polars.Float64,
            "rot": polars.Float64,
            "ra": polars.Float64,
            "dec": polars.Float64,
            "done": polars.Boolean,
            "mjd": polars.Int32,
            "gimg": polars.Int32,
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
        },
    )
    data_df.write_parquet(output_file)


def log_reply(reply: Reply):
    """Logs a reply to a command."""

    log.debug(reply.string)
