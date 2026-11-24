"""Explicit GNSS epoch-label to Astropy Time conversion.

RINEX/SP3 epochs are labels in the product's declared time system, not
necessarily UTC.  The existing attosecond helpers remain the precision-safe
representation for SP3 interpolation; Astropy handles time *scales* whenever
we request celestial coordinates, Earth orientation, or tidal corrections.
"""

from __future__ import annotations

from datetime import timedelta

from astropy.time import Time, TimeDelta

from dsoclasses.time.pyattotime import Epoch, at2pt, to_attoseconds


def to_astropy_time(epoch, time_system: str = 'GPS') -> Time:
    """Interpret an epoch label on its declared scale (not as implicit UTC).

    Supported labels are GPS, UTC, TAI, TT, BDT and GLO.  The GLONASS
    conversion approximates UTC(SU) as UTC, suitable for the teaching data;
    other unspecified scale offsets are deliberately not guessed.

    An attodatetime's sub-microsecond remainder is added separately so
    conversion never first rounds the scientific label to milliseconds.
    """
    if not isinstance(time_system, str):
        raise TypeError('time_system must be a declared time-scale string')
    system = time_system.strip().upper()
    total_as = to_attoseconds(epoch)
    us_as = (total_as // 10**12) * 10**12
    base = at2pt(Epoch(us_as))
    remainder = total_as - us_as

    if system in ('GPS', 'QZS'):
        # GPS time = TAI - 19 seconds.  Astropy's Time has no 'gps' scale;
        # its GPS *format* is a different API for seconds since 1980-01-06.
        result = Time(base, scale='tai') + TimeDelta(19, format='sec')
    elif system == 'BDT':
        # BeiDou Time = TAI - 33 seconds.
        result = Time(base, scale='tai') + TimeDelta(33, format='sec')
    elif system in ('UTC', 'TAI', 'TT'):
        result = Time(base, scale=system.lower())
    elif system in ('GLO', 'GLONASS'):
        # GLONASS time is conventionally UTC(SU) + 3 hours.  The small
        # UTC(SU)-UTC difference is beyond this introductory scope.
        result = Time(base - timedelta(hours=3), scale='utc')
    else:
        raise ValueError(f'Unsupported time scale {time_system!r}; do not assume UTC')

    if remainder:
        result += TimeDelta(remainder / 10**18, format='sec')
    return result
