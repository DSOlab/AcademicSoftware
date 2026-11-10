"""Centralized, high-resolution handling of GNSS epoch *labels*.

All GNSS modules should go through this module when comparing, interpolating,
or shifting an epoch. A timestamp is normalized to **a Python integer number
of attoseconds** (10**-18 seconds) since the proleptic Gregorian
1970-01-01 00:00:00. Python integers do not overflow for ordinary dates.

The input label's time scale is deliberately left unchanged: GPS, UTC, TAI,
and GLONASS time labels must NOT be mixed without an explicit conversion.
Timezone-aware datetimes are rejected for this reason. In particular, the
"Unix epoch" above is only a convenient calendar origin, not an implicit UTC
or Unix/GPS-timescale conversion.

Accepted public inputs:
    * datetime.datetime (whole microseconds);
    * numpy.datetime64 (up to nanoseconds);
    * attotime.attodatetime (Decimal fractional nanoseconds, up to attoseconds);
    * documented absolute-attosecond objects with to_attoseconds();
    * Epoch (the internal wrapper for an already-normalized attosecond integer).

Use at2pt() only for plotting or display: a Python datetime cannot represent
fractions of a microsecond. Never use it in scientific interpolation.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
from decimal import Decimal, InvalidOperation, ROUND_HALF_EVEN
from typing import Any
import warnings

import numpy as np

try:
    import attotime as _attotime
except ImportError:  # The basic teaching notebooks also run without attotime.
    _attotime = None

_ATTO_PER_SECOND = 10**18
_ATTO_PER_MICROSECOND = 10**12
_ATTO_PER_NANOSECOND = 10**9
_EPOCH = datetime(1970, 1, 1)


@dataclass(frozen=True)
class Epoch:
    """A normalized time label with integer attosecond precision.

    This wrapper distinguishes *absolute attoseconds* from a plain integer,
    historically accepted by to_attoseconds() as Unix-style seconds. It is
    useful when passing time values between internal modules without creating
    a new third-party datetime instance for every sample.
    """

    attoseconds: int

    def __post_init__(self) -> None:
        if isinstance(self.attoseconds, bool) or not isinstance(self.attoseconds, int):
            raise TypeError('Epoch.attoseconds must be an integer')


def _decimal(value: Any, *, what: str) -> Decimal:
    """Convert a short numerical quantity without introducing binary rounding.

    Decimal(str(float)) represents the intended user-facing decimal (unlike
    Decimal(float)); an existing Decimal is used directly. Reject NaN/inf.
    """
    try:
        result = value if isinstance(value, Decimal) else Decimal(str(value))
    except (InvalidOperation, ValueError, TypeError) as exc:
        raise TypeError(f'{what} must be a finite number') from exc
    if not result.is_finite():
        raise ValueError(f'{what} must be finite')
    return result


def seconds_to_attoseconds(value: Any) -> int:
    """Convert seconds/duration to integer attoseconds, rounding at 1 as.

    Never multiply a floating-point absolute timestamp by 10**18. A *short*
    duration may be supplied as float, str, int or Decimal; converting through
    a decimal string retains nanosecond-scale values accurately.
    """
    attoseconds = _decimal(value, what='Seconds') * _ATTO_PER_SECOND
    return int(attoseconds.to_integral_value(rounding=ROUND_HALF_EVEN))


def _calendar_to_attoseconds(value: Any) -> int:
    """Read calendar fields before any lossy to_datetime()/timetuple() path.

    attotime's `nanosecond` is a Decimal **fraction of a microsecond** in
    [0, 1000), *not* a full nanosecond field or an absolute timestamp. Hence
    `microsecond * 1e12 + nanosecond * 1e9` is the correct subsecond total.
    """
    tz = getattr(value, 'tzinfo', None)
    if tz is not None and getattr(value, 'utcoffset', lambda: None)() is not None:
        raise ValueError('Use a timezone-naive epoch in the declared GNSS time system')
    base = datetime(int(value.year), int(value.month), int(value.day),
                    int(value.hour), int(value.minute), int(value.second))
    days = (base - _EPOCH).days
    second_of_day = base.hour * 3600 + base.minute * 60 + base.second
    microseconds = int(getattr(value, 'microsecond', 0))
    if not 0 <= microseconds < 1_000_000:
        raise ValueError('Invalid microsecond field')
    # The actual attotime library exposes a Decimal here, e.g. 123.456789012.
    nanoseconds = _decimal(getattr(value, 'nanosecond', 0), what='nanosecond field')
    if not 0 <= nanoseconds < 1000:
        raise ValueError('nanosecond field must be in [0, 1000)')
    remainder = nanoseconds * _ATTO_PER_NANOSECOND
    if remainder != remainder.to_integral_value():
        raise ValueError('Fractional nanoseconds exceed attosecond precision')
    return ((days * 86400 + second_of_day) * _ATTO_PER_SECOND
            + microseconds * _ATTO_PER_MICROSECOND + int(remainder))


def to_attoseconds(value: Any) -> int:
    """Convert a supported epoch label to absolute integer attoseconds.

    No float Unix timestamps, guessed time-scale conversions, or implicit
    truncation to milliseconds/microseconds/nanoseconds are performed.
    """
    if isinstance(value, Epoch):
        return value.attoseconds
    if isinstance(value, datetime):
        return _calendar_to_attoseconds(value)
    if isinstance(value, np.datetime64):
        if np.isnat(value):
            raise ValueError('NaT cannot be used as a GNSS epoch')
        # datetime64 uses signed integer ticks. Check that nanosecond casting
        # roundtrips, rather than allowing numpy's silent date-range overflow.
        nanoseconds = value.astype('datetime64[ns]')
        if nanoseconds.astype(value.dtype) != value:
            raise ValueError('numpy datetime64 epoch is outside nanosecond range')
        return int(nanoseconds.astype('int64')) * _ATTO_PER_NANOSECOND
    if all(hasattr(value, name) for name in
           ('year', 'month', 'day', 'hour', 'minute', 'second')):
        return _calendar_to_attoseconds(value)
    # Some high-resolution timestamp types expose a *total* since 1970;
    # do not read `attoseconds` as a property, since it may mean only the
    # fractional part of one second in other APIs.
    for name in ('to_attoseconds', 'as_attoseconds', 'to_asec'):
        method = getattr(value, name, None)
        if method is not None:
            result = method() if callable(method) else method
            if isinstance(result, bool) or not isinstance(result, (int, np.integer)):
                raise TypeError(f'{name}() must return integral absolute attoseconds')
            return int(result)
    if hasattr(value, 'sec') and hasattr(value, 'asec'):
        # Known representation: absolute seconds plus *subsecond* attoseconds.
        return int(value.sec) * _ATTO_PER_SECOND + int(value.asec)
    if isinstance(value, (int, float, Decimal, np.integer, np.floating)) and not isinstance(value, bool):
        # Historical behavior: bare numeric inputs are *seconds* since 1970.
        return seconds_to_attoseconds(value)
    raise TypeError(f'Unsupported GNSS epoch type: {type(value).__name__}')


def datetime_to_attoseconds(value: datetime) -> int:
    """Compatibility entry point for native datetime values."""
    if not isinstance(value, datetime):
        raise TypeError('Expected datetime.datetime')
    return to_attoseconds(value)


def from_attoseconds(value: int):
    """Build a datetime-like object from an integer attosecond epoch.

    With attotime installed, returns an attodatetime and retains every
    attosecond. Without it, falls back to native datetime with rounding to
    microseconds; callers needing sub-microsecond arithmetic must install
    attotime (or work directly with Epoch and integer attoseconds).
    """
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError('Absolute attoseconds must be an integer')
    whole, subsecond = divmod(value, _ATTO_PER_SECOND)
    base = _EPOCH + timedelta(seconds=whole)
    microsecond, remaining_as = divmod(subsecond, _ATTO_PER_MICROSECOND)
    if _attotime is not None:
        nanosecond = Decimal(remaining_as) / _ATTO_PER_NANOSECOND
        return _attotime.attodatetime(
            base.year, base.month, base.day, base.hour, base.minute,
            base.second, microsecond, nanosecond,
        )
    # Round to the nearest microsecond, including carries into next second.
    rounded = int((Decimal(subsecond) / _ATTO_PER_MICROSECOND).to_integral_value(
        rounding=ROUND_HALF_EVEN))
    if remaining_as:
        warnings.warn('attotime is not installed: returning datetime rounded to microseconds',
                      RuntimeWarning, stacklevel=2)
    return base + timedelta(microseconds=rounded)


def shift_epoch_seconds(epoch: Any, seconds: Any):
    """Return epoch + seconds without intermediate floating Unix timestamps.

    This is the common path for signal-transmission-time and numerical
    differentiation in the GNSS algorithms. With attotime installed the
    resulting epoch keeps attosecond precision even if the input was native
    datetime or numpy.datetime64.
    """
    return from_attoseconds(to_attoseconds(epoch) + seconds_to_attoseconds(seconds))


def at2pt(value: Any) -> datetime:
    """Convert an epoch to native Python datetime for plotting **only**.

    Fractions below a microsecond are necessarily lost. It is never a safe
    normalization step before SP3 interpolation or light-time iteration.
    """
    if isinstance(value, datetime):
        return value
    # Round down to a whole microsecond for display, not model arithmetic.
    total = to_attoseconds(value)
    microseconds = total // _ATTO_PER_MICROSECOND
    return _EPOCH + timedelta(microseconds=microseconds)


def fsec2asec(seconds: Any):
    """Build an attotimedelta if possible, native timedelta otherwise.

    Kept for compatibility with the original notebooks. Prefer
    shift_epoch_seconds() in new code, so the same conversion is used for
    datetime, datetime64 and attodatetime input epochs.
    """
    asec = seconds_to_attoseconds(seconds)
    if _attotime is not None:
        return _attotime.attotimedelta(nanoseconds=Decimal(asec) / _ATTO_PER_NANOSECOND)
    microseconds = int((Decimal(asec) / _ATTO_PER_MICROSECOND).to_integral_value(
        rounding=ROUND_HALF_EVEN))
    return timedelta(microseconds=microseconds)
