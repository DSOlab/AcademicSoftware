import datetime
from math import floor

month_day = [
    [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365],
    [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366],
]

mtab = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def cal2mjd(t):
    """Return (integer MJD, float seconds-of-day) for any supported GNSS epoch.

    Epoch normalization is shared with RINEX, SP3 and signal propagation.
    Doing the day split in integer attoseconds correctly handles both dates
    before 1970 and attodatetime's fractional nanoseconds. Only the final
    seconds-of-day is converted to float because the historical API returns
    floats, rather than claiming attosecond precision for an MJD float.

    This is calendar conversion, NOT conversion between GNSS time scales.
    """
    from dsoclasses.time.pyattotime import to_attoseconds

    asec_per_day = 86400 * 10**18
    days, remainder = divmod(to_attoseconds(t), asec_per_day)
    return 40587 + days, remainder / 10**18


def cal2fmjd(t):
    """Transform a calendar date to MJD and fraction of day.

    The input parameter t, can be either
    * a native python datetime instance, or
    * an attodatetime instance

    Return:
        mjd -> the (integral) MJD + the fraction of day
        The fraction of day is always in range [0, 1)
    """
    mjd, sec = cal2mjd(t)
    return float(mjd) + sec / 86400.0


def mjd2cal(mjd: float) -> datetime.datetime:
    epoch = datetime.datetime(1858, 11, 17)  # MJD 0
    days = floor(mjd)
    frac = mjd - days  # in [0,1)
    return epoch + datetime.timedelta(days=days, seconds=frac * 86400.0)


def to_utc(dt: datetime) -> datetime:
    """Return dt in UTC. Naive datetimes are assumed to already be UTC."""
    if dt.tzinfo is None:
        return dt.replace(tzinfo=datetime.timezone.utc)
    return dt.astimezone(datetime.timezone.utc)
