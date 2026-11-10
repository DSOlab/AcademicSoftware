#! /usr/bin/python

import math
import numpy as np


# Cartesian coordinates to (r, Azimouth, Elevation) from rsta to rsat in ([m], [rad], [rad])
def azele(rsat, rsta):
    # topocentric matrix
    lat, lon, hgt = car2ell(*rsta)
    R = geodetic2lvlh(lat, lon).transpose()
    # Require valid ECEF vectors, including for callers passing Python lists.
    rsat = np.asarray(rsat, dtype=float)
    rsta = np.asarray(rsta, dtype=float)
    if rsat.shape != (3,) or rsta.shape != (3,) or not np.all(np.isfinite(rsat)) or not np.all(np.isfinite(rsta)):
        raise ValueError("rsat and rsta must be finite ECEF XYZ vectors")
    enu = R @ (rsat - rsta)
    r = np.linalg.norm(enu)
    if r == 0:
        raise ValueError("Satellite and station positions must differ")
    az = np.arctan2(enu[0], enu[1])
    el = np.arcsin(np.clip(enu[2] / r, -1.0, 1.0))
    return r, az, el


""" The resulting matrix from this function, R, can be used in the sense:
    |δX|        |e|
    |δY|  = R * |n|
    |δZ|        |u|
"""


def geodetic2lvlh(lat, lon):
    enu = np.identity(3)
    sphi = math.sin(lat)
    cphi = math.cos(lat)
    slmb = math.sin(lon)
    clmb = math.cos(lon)
    # unit vector e, i.e. east
    enu[0][0] = -slmb
    enu[1][0] = clmb
    enu[2][0] = 0e0
    # unit vector n, i.e. north
    enu[0][1] = -sphi * clmb
    enu[1][1] = -sphi * slmb
    enu[2][1] = cphi
    # unit vector u, i.e. up
    enu[0][2] = cphi * clmb
    enu[1][2] = cphi * slmb
    enu[2][2] = sphi
    return enu


def ell2car(lat, lon, h, a=6378137e0, f=1e0 / 298.257223563e0):
    """Convert geodetic latitude/longitude [rad] and height [m] to ECEF XYZ."""
    if not np.all(np.isfinite([lat, lon, h, a, f])):
        raise ValueError("Geodetic coordinates and ellipsoid must be finite")
    if not (-np.pi/2 <= lat <= np.pi/2) or a <= 0 or not (0 < f < 1):
        raise ValueError("Invalid geodetic latitude or ellipsoid")
    e2 = (2e0 - f) * f  # eccentricity squared
    sf = math.sin(lat)
    cf = math.cos(lat)
    sl = math.sin(lon)
    cl = math.cos(lon)
    Rn = a / math.sqrt(1e0 - e2 * sf * sf)  # normal radius of curvature
    x = (Rn + h) * cf * cl
    y = (Rn + h) * cf * sl
    z = ((1e0 - e2) * Rn + h) * sf
    return x, y, z


def car2ell(x, y, z, a=6378137e0, f=1e0 / 298.257223563e0):
    """Convert ECEF XYZ metres to geodetic latitude, longitude, height.

    The geocentre has no unique geodetic latitude/longitude; reject it
    explicitly rather than returning an arbitrary polar coordinate.
    """
    if not np.all(np.isfinite([x, y, z, a, f])) or a <= 0 or not (0 < f < 1):
        raise ValueError("ECEF coordinates and ellipsoid must be finite and valid")
    if x == y == z == 0:
        raise ValueError("Geodetic coordinates are undefined at the Earth centre")
    aeps2 = a * a * 1e-32
    e2 = (2e0 - f) * f  # eccentricity squared
    e4t = e2 * e2 * 1.5e0
    ep2 = 1.0e0 - e2
    ep = math.sqrt(ep2)
    aep = a * ep

    # Compute distance from polar axis squared
    p2 = x * x + y * y

    # Compute longitude.
    lon = math.atan2(y, x) if p2 != 0e0 else 0e0

    # Ensure that Z-coordinate is unsigned.
    absz = abs(z)

    if p2 > aeps2:
        # Continue unless at the poles
        # Compute distance from polar axis.
        p = math.sqrt(p2)
        # Normalize.
        s0 = absz / a
        pn = p / a
        zp = ep * s0
        # Prepare Newton correction factors.
        c0 = ep * pn
        c02 = c0 * c0
        c03 = c02 * c0
        s02 = s0 * s0
        s03 = s02 * s0
        a02 = c02 + s02
        a0 = math.sqrt(a02)
        a03 = a02 * a0
        d0 = zp * a03 + e2 * s03
        f0 = pn * a03 - e2 * c03
        # Prepare Halley correction factor.
        b0 = e4t * s02 * c02 * pn * (a0 - ep)
        s1 = d0 * f0 - b0 * s0
        cp = ep * (f0 * f0 - b0 * c0)
        # Evaluate latitude and height.
        lat = math.atan(s1 / cp)
        s12 = s1 * s1
        cp2 = cp * cp
        hgt = (p * cp + absz * s1 - a * math.sqrt(ep2 * s12 + cp2)) / math.sqrt(
            s12 + cp2
        )
    else:
        # Special case: pole.
        lat = math.pi / 2e0
        hgt = absz - aep

    # Restore sign of latitude.
    if z < 0.0e0:
        lat = -lat
    return lat, lon, hgt
