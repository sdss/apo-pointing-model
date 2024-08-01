#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-07-31
# @Filename: __main__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import logging
import pathlib

import click

from sdsstools.daemonizer import cli_coro


@click.group()
@click.option("-v", "--verbose", is_flag=True, help="Shows verbose output.")
def apo_pointing_model(verbose: bool = False):
    """Command-line interface to the APO pointing model."""

    from apo_pointing_model import log

    if verbose:
        log.sh.setLevel(logging.DEBUG)
    else:
        log.sh.setLevel(logging.INFO)


@apo_pointing_model.command()
@click.argument("OUTPUT_FILE", type=click.Path(dir_okay=True))
@click.argument("N_POINTS", type=int, required=False)
@click.option(
    "-a",
    "--alt-range",
    nargs=2,
    type=float,
    default=(28, 85),
    show_default=True,
    help="Altitude range.",
)
@click.option(
    "-z",
    "--az-range",
    nargs=2,
    type=float,
    default=(0, 359),
    show_default=True,
    help="Azimuth range.",
)
@cli_coro()
async def run(
    n_points: int | None,
    output_file: str,
    alt_range: tuple[float, float],
    az_range: tuple[float, float],
):
    """Collects APO pointing model data."""

    from apo_pointing_model.runner import get_pointing_data

    path_ = pathlib.Path(output_file)
    if not path_.exists() and n_points is None:
        raise click.BadArgumentUsage(
            "N_POINTS is required if OUTPUT_FILE does not exist."
        )

    await get_pointing_data(
        n_points,
        output_file,
        alt_range=alt_range,
        az_range=az_range,
    )
