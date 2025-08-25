import attotime
import datetime
import numpy as np
from typing import Any


def at2pt(at):
    """Convert an instance of attotime to native datetime.datetime

    Warning! This will cuse loss of precision.
    Translate an attotime instance to a native python datetime instance.
    """
    return datetime.datetime(
        at.year, at.month, at.day, at.hour, at.minute, at.second, at.microsecond
    )


_ATTO_PER_SEC = 10**18
_ATTO_PER_NS = 10**9
_NS_PER_SEC   = 10**9

def datetime_to_attoseconds(dt: datetime.datetime) -> int:
    ns = np.datetime64(dt, "ns").astype("int64")  # exact integer ns since epoch
    sec = int(ns // _NS_PER_SEC)
    rem_ns = int(ns - sec * _NS_PER_SEC)
    return sec * _ATTO_PER_SEC + rem_ns * _ATTO_PER_NS


def to_attoseconds(t: Any) -> int:
    # 1) direct attotime-style methods/props
    for attr in ("to_attoseconds", "attoseconds", "as_attoseconds", "to_asec"):
        if hasattr(t, attr):
            v = getattr(t, attr)
            try:
                return int(v() if callable(v) else v)
            except Exception:
                pass

    # 2) common field patterns
    if hasattr(t, "sec") and hasattr(t, "asec"):
        return int(t.sec) * _ATTO_PER_SEC + int(t.asec)

    # 3) attodatetime-style conversions -> try getting a regular datetime
    for conv in ("to_datetime", "as_datetime", "datetime"):
        if hasattr(t, conv):
            dt = getattr(t, conv)
            dt = dt() if callable(dt) else dt
            if isinstance(dt, datetime.datetime):
                return datetime_to_attoseconds(dt)

    # 4) attodatetime-style components (year..second + subsecond fields)
    if all(hasattr(t, a) for a in ("year","month","day","hour","minute","second")):
        base = datetime.datetime(int(t.year), int(t.month), int(t.day),
                            int(t.hour), int(t.minute), int(t.second))
        asec = datetime_to_attoseconds(base)
        # add microseconds / nanoseconds / attoseconds if present
        if hasattr(t, "microsecond"):   asec += int(getattr(t, "microsecond")) * (10**12)
        if hasattr(t, "microseconds"):  asec += int(getattr(t, "microseconds")) * (10**12)
        if hasattr(t, "nanosecond"):    asec += int(getattr(t, "nanosecond"))   * (10**9)
        if hasattr(t, "nanoseconds"):   asec += int(getattr(t, "nanoseconds"))  * (10**9)
        if hasattr(t, "attosecond"):    asec += int(getattr(t, "attosecond"))
        if hasattr(t, "attoseconds"):   asec += int(getattr(t, "attoseconds"))
        return asec

    # 5) numpy.datetime64
    if isinstance(t, np.datetime64):
        ns = np.datetime64(t, "ns").astype("int64")
        return int(ns) * _ATTO_PER_NS

    # 6) python datetime
    if isinstance(t, datetime.datetime):
        return datetime_to_attoseconds(t)

    # 7) numeric seconds
    if isinstance(t, (int, float, np.integer, np.floating)):
        return int(round(float(t) * _ATTO_PER_SEC))

    raise TypeError(f"Unsupported time type: {t.__class__.__module__}.{t.__class__.__qualname__}  repr={t!r}")


def fsec2asec(fsec):
    isec = int(fsec)  # integral seconds
    imsec = int((fsec - isec) * 1e6)  # integral microseconds
    fnsec = float(fsec * 1e9 - imsec * 1e3)  # fractional nanoseconds
    assert (
        abs(
            float(
                attotime.attotimedelta(
                    seconds=isec, microseconds=imsec, nanoseconds=fnsec
                ).total_nanoseconds()
            )
            - fsec * 1e9
        )
        < 1e-1
    )
    return attotime.attotimedelta(seconds=isec, microseconds=imsec, nanoseconds=fnsec)
