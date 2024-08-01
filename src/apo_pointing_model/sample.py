#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-07-31
# @Filename: sample.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import numpy
from astropy import units as uu
from astropy.coordinates import (
    AltAz,
    EarthLocation,
    SkyCoord,
    uniform_spherical_random_surface,
)
from astropy.time import Time


__all__ = ["get_random_sample", "to_icrs"]


def get_random_sample(
    n_points: int,
    alt_range: tuple[float, float] | None = None,
    az_range: tuple[float, float] | None = None,
    min_alt: float = 30,
):
    """Provides a random sample of RA/Dec points on the surface of a sphere.

    Parameters
    ----------
    n_points
        Number of points to return. This number is ensured even if ``alt_range``
        or ``az_range`` are provided.
    alt_range
        The range of altitude to which to limit the sample, in degrees.
    az_range
        The range of azimuth to which to limit the sample, in degrees.
    min_alt
        The minimum altitude in degrees to consider a point.

    Returns
    -------
    coordinates
        A 2D array with the Alt/Az coordinates of the points on the sphere.

    """

    points = numpy.zeros((0, 2), dtype=numpy.float64)

    lco = EarthLocation.of_site("Las Campanas Observatory")
    now = Time.now()

    while True:
        sph_points = uniform_spherical_random_surface(n_points)
        altaz = AltAz(
            alt=sph_points.lat.deg * uu.deg,
            az=sph_points.lon.deg * uu.deg,
            location=lco,
            obstime=now,
        )

        if alt_range is not None:
            alt = altaz.alt.deg
            altaz = altaz[(alt > alt_range[0]) & (alt < alt_range[1])]
        if az_range is not None:
            az = altaz.az.deg
            altaz = altaz[(az > az_range[0]) & (az < az_range[1])]

        altaz = SkyCoord(altaz[altaz.alt.deg > min_alt])

        altaz_array = numpy.array([altaz.alt.deg, altaz.az.deg], dtype=numpy.float64).T

        points = numpy.vstack((points, altaz_array))

        if points.shape[0] >= n_points:
            return points[0:n_points, :]


def to_icrs(alta: float, az: float) -> SkyCoord:
    """Converts Alt/Az to ICRS coordinates."""

    apo = EarthLocation.of_site("Apache Point Observatory")
    now = Time.now()

    altaz = SkyCoord(
        AltAz(
            alt=alta * uu.deg,
            az=az * uu.deg,
            location=apo,
            obstime=now,
        )
    )

    icrs = altaz.transform_to("icrs")

    return icrs.ra.deg, icrs.dec.deg
