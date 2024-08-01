#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-07-31
# @Filename: __init__.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import warnings

from clu.exceptions import CluWarning

from sdsstools import get_logger, get_package_version


warnings.simplefilter("ignore", category=CluWarning)


NAME = "sdss-apo-pointing-model"

log = get_logger(NAME)

__version__ = get_package_version(path=__file__, package_name=NAME)
