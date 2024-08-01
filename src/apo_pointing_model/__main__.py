#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-07-31
# @Filename: __main__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import logging

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
@click.argument("N_POINTS", type=int)
@click.argument("OUTPUT_FILE", type=click.Path(dir_okay=True))
@cli_coro()
async def run(n_points: int, output_file: str):
    """Collects APO pointing model data."""

    from apo_pointing_model.runner import get_pointing_data

    await get_pointing_data(n_points, output_file)
