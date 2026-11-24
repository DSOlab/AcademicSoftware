"""Small Astropy/pyTMD interfaces for the GNSS and DORIS lessons.

All Cartesian input/output coordinates are metres.  Positions are specified
in the terrestrial system of the user's orbit/reference products; this module
does not transform between distinct ITRF realizations.  A station solid-Earth
tide is a *displacement*, not a new reference-frame realization.
"""

from __future__ import annotations

import math

import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation, get_body_barycentric_posvel, solar_system_ephemeris

from dsoclasses.geodesy.transformations import car2ell, geodetic2lvlh
from dsoclasses.time.astropy_time import to_astropy_time


def _xyz(value):
    xyz = np.asarray(value, dtype=float)
    if xyz.shape != (3,) or not np.all(np.isfinite(xyz)):
        raise ValueError('Expected a finite terrestrial XYZ vector in metres')
    return xyz


def itrs_to_gcrs_xyz(station_xyz, epoch, time_system='GPS'):
    """Convert one fixed terrestrial XYZ vector into GCRS [metres].

    Astropy uses IERS Earth-orientation parameters (UT1/polar motion).
    Station tidal displacements, if needed, must be applied beforehand.
    """
    x, y, z = _xyz(station_xyz)
    site = EarthLocation.from_geocentric(x, y, z, unit=u.m)
    return np.asarray(site.get_gcrs(to_astropy_time(epoch, time_system))
                      .cartesian.xyz.to_value(u.m), dtype=float)


def sun_moon_geocentric_xyz(epoch, time_system='GPS', ephemeris='builtin'):
    """Geocentric Sun and Moon vectors [metres], in ICRS-aligned axes.

    These are geometric barycentric ephemeris differences, suitable for
    teaching; they are not light-time-corrected apparent sky directions.
    """
    time = to_astropy_time(epoch, time_system)
    with solar_system_ephemeris.set(ephemeris):
        earth, _ = get_body_barycentric_posvel('earth', time)
        sun, _ = get_body_barycentric_posvel('sun', time)
        moon, _ = get_body_barycentric_posvel('moon', time)
    return ((sun - earth).xyz.to_value(u.m),
            (moon - earth).xyz.to_value(u.m))


def solid_earth_tide_xyz(station_xyz, epoch, time_system='GPS'):
    """IERS-style solid Earth tide displacement, ECEF XYZ [m], via pyTMD.

    Uses pyTMD's ephemeris-based model and a tide-free convention.  The
    returned local 'R' component is geocentric radial; mapping it onto the
    geodetic Up direction has a submillimetre-to-millimetre approximation
    for ordinary terrestrial stations, sufficient for a code-only lesson.
    """
    import pyTMD

    position = _xyz(station_xyz)
    lat, lon, _ = car2ell(*position)
    # pyTMD's 'datetime' input convention means UTC numpy datetime64.
    utc_time = to_astropy_time(epoch, time_system).utc.datetime
    displacement = pyTMD.compute.SET_displacements(
        np.array([math.degrees(lon)]),
        np.array([math.degrees(lat)]),
        np.array([np.datetime64(utc_time)]),
        method='ephemerides', standard='datetime', type='trajectory',
        tide_system='tide_free', variable=['N', 'E', 'R'],
    )
    enu = np.array([
        float(np.asarray(displacement['E']).reshape(-1)[0]),
        float(np.asarray(displacement['N']).reshape(-1)[0]),
        float(np.asarray(displacement['R']).reshape(-1)[0]),
    ])
    if not np.all(np.isfinite(enu)):
        raise ValueError('pyTMD returned non-finite solid Earth tides')
    return geodetic2lvlh(lat, lon) @ enu
