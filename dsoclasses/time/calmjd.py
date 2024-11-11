import datetime

month_day = [
    [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365],
    [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
]

mtab = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

def cal2mjd(t):
    """ Transform a calendar date to MJD and seconds of day.

        The input parameter t, can be either
        * a native python datetime instance, or
        * an attodatetime instance

        Return:
            mjd    -> the (integral) MJD and
            secday -> the (fractional) seconds of day
    """
    try:
        secday = t.nanosecond * 1e-9
    except:
        secday = 0.
    secday += (float(t.hour*3600) + (float(t.minute*60) + (float(t.second) + float(t.microsecond) * 1e-6)))

    yr = int(t.year)
    mn = int(t.month)
    dm = int(t.day)
    
    assert mn >= 1 and mn <= 12
    leap = int((mn == 2) and not (yr % 4) and (yr % 100 or not (yr % 400)))
    assert not ((dm < 1) or (dm > (mtab[mn - 1] + leap)))

# took a while to debug this ... negative integer division is different from C
    my = int( (mn - 14) / 12 )
    iypmy = yr + my

    mjd = (1461 * (iypmy + 4800)) // 4 + (367 * (mn - 2 - 12 * my)) // 12 - \
        (3 * ((iypmy + 4900) // 100)) // 4 + dm - 2432076
    
    return mjd, secday

def cal2fmjd(t):
    """ Transform a calendar date to MJD and fraction of day.

        The input parameter t, can be either
        * a native python datetime instance, or
        * an attodatetime instance

        Return:
            mjd -> the (integral) MJD + the fraction of day
            The fraction of day is always in range [0, 1)
    """
    mjd, sec = cal2mjd(t)
    return float(mjd) + sec / 86400.
