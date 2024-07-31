#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-07-31
# @Filename: __main__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import click


@click.group()
def apo_pointing_model():
    """Command-line interface to the APO pointing model."""

    pass
