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


def datetime_to_attoseconds(dt: datetime) -> int:
    if dt.tzinfo is not None:
        dt = dt.astimezone(datetime.timezone.utc).replace(tzinfo=None)
    ns = np.datetime64(dt, "ns").astype("int64")
    sec = int(ns // 1_000_000_000)
    rem_ns = int(ns - sec * 1_000_000_000)
    return sec * _ATTO_PER_SEC + rem_ns * _ATTO_PER_NS


def to_attoseconds(t) -> int:
    """
    Normalize supported time-like inputs to integer attoseconds:
      - attotime-like objects with .to_attoseconds() / .attoseconds / .as_attoseconds()
      - objects with (sec, asec) attributes
      - Python datetime.datetime
      - numpy.datetime64
    """
    for attr in ("to_attoseconds", "attoseconds", "as_attoseconds", "to_asec"):
        if hasattr(t, attr):
            v = getattr(t, attr)
            return int(v() if callable(v) else v)
    if hasattr(t, "sec") and hasattr(t, "asec"):
        return int(t.sec) * _ATTO_PER_SEC + int(t.asec)
    if isinstance(t, datetime):
        return datetime_to_attoseconds(t)
    if isinstance(t, np.datetime64):
        ns = np.datetime64(t, "ns").astype("int64")
        return int(ns) * _ATTO_PER_NS
    raise TypeError("Unsupported time type")


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
