#! /usr/bin/python

import math
import numpy as np

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
    e2 = (2e0-f)*f # eccentricity squared
    sf = sin(lat)
    cf = cos(lat)
    sl = sin(lon)
    cl = cos(lon)
    Rn = a / math.sqrt(1e0 - e2 * sf * sf) # normal radius of curvature
    x = (Rn + h) * cf * cl
    y = (Rn + h) * cf * sl
    z = ((1e0 - e2) * Rn + h) * sf
    return x, y, z

def car2ell(x, y, z, a=6378137e0, f=1e0 / 298.257223563e0):
    aeps2 = a * a * 1e-32
    e2 = (2e0-f)*f # eccentricity squared
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

    if (p2 > aeps2):
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
        hgt = (p * cp + absz * s1 - a * math.sqrt(ep2 * s12 + cp2)) /\
          math.sqrt(s12 + cp2)
    else:
        # Special case: pole.
        lat = math.pi / 2e0
        hgt = absz - aep

# Restore sign of latitude.
    if z < 0.e0: lat = -lat
    return lat, lon, hgt
