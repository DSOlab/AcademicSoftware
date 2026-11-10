"""Convenience wrapper for the GPT3/GMF tropospheric delay model.

This computation is deliberately stateless: cached receiver coordinates on a
function object can give wrong answers for a different station, or when the
function is called concurrently by two notebooks.
"""
import math

import numpy as np

from dsoclasses.geodesy import transformations
from dsoclasses.troposphere import gmf, gpt3


def tropo_delay(rsta, t, el, gpt3_grid):
    """Slant tropospheric delay [m] for receiver XYZ [m], elevation [rad]."""
    station = np.asarray(rsta, dtype=float)
    if station.shape != (3,) or not np.all(np.isfinite(station)):
        raise ValueError('rsta must be a finite 3-vector in ECEF metres')
    if not math.isfinite(float(el)) or not (0 < el <= math.pi / 2):
        raise ValueError('Tropospheric delay requires 0 < elevation <= 90 degrees')
    lat, lon, hgt = transformations.car2ell(*station)
    # gpt3 expects longitude first; GMF expects latitude first.
    meteo = gpt3.gpt3(t, lon, lat, hgt, gpt3_grid)
    zhd = gpt3.saastamoinen_zhd(lat, hgt, meteo['p'])
    zwd = gpt3.askne_zwd(meteo['e'], meteo['Tm'], meteo['la'])
    gmfh, gmfw = gmf.gmf(t, lat, lon, hgt, math.pi / 2 - el)
    return float(zhd * gmfh + zwd * gmfw)
