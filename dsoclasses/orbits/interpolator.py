"""Conservative, local interpolation of precise SP3 positions and clocks.

The reader and numerical interpolator intentionally have separate jobs:

* :class:`dsoclasses.orbits.sp3c.Sp3` parses SP3-c/SP3-d files, units and flags.
* :class:`OrbitInterpolator` interpolates ONE satellite's sampled time series.
* :class:`Sp3Interpolator` builds one orbit interpolator per satellite.

There is no global polynomial fit across an entire day and no silent
extrapolation.  Missing clock samples are independent of available positions
unless ``exclude_missing_clock_values=True`` explicitly requests otherwise.

SP3 timestamps are labels in the file's declared time system (e.g. GPS or
TAI).  This module does NOT convert between time systems or apply leap seconds.
"""

from __future__ import annotations

from bisect import bisect_left
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.interpolate import BarycentricInterpolator, CubicSpline, PchipInterpolator

from dsoclasses.time.pyattotime import Epoch, to_attoseconds, seconds_to_attoseconds


_ATTO_PER_SECOND = 10**18
_METHODS = frozenset({"Linear", "Polynomial", "Barycentric", "CubicSpline", "PchipInterpolator"})


def _as_nanoseconds(value: Any) -> int:
    """Legacy helper; use to_attoseconds() for lossless internal arithmetic.

    Its nanosecond return value is retained only for compatibility with
    older importers of this private function. It is not used by interpolation.
    """
    return to_attoseconds(value) // 10**9


def flag_is_on(flags: str, selected: Sequence[str]) -> bool:
    """Test whether any requested SP3 event code is present in a record."""
    return bool(set(str(flags).upper()) & {code.upper() for code in selected})


class OrbitInterpolator:
    """Interpolate one satellite's position and (optionally) clock samples.

    Parameters
    ----------
    satid:
        SP3 identifier, such as ``'L40'`` or ``'G02'``.
    tarray, x, y, z, clk:
        Sample epochs and SI-valued series.  Non-finite rows are unavailable.
        ``clk=None`` means there are no clock observations.
    interval_in_sec:
        Half-width of the local data window, not a full interval width.
    min_data_pts:
        Minimum count of *valid* samples required for interpolation between
        tabulated epochs.  An exact available tabulated position is returned
        even when the window contains fewer than this number of samples.
    itype:
        ``'Linear'``, ``'CubicSpline'``, ``'PchipInterpolator'``, or a local
        ``'Polynomial'`` / ``'Barycentric'`` fit.  Polynomial now means a real
        polynomial, rather than the former (misleading) alias for linear.
    stencil_pts:
        Number of nearby nodes used in polynomial/barycentric interpolation.
        A small local stencil avoids fitting high-degree polynomials to a
        whole orbital arc.  Must be at least ``min_data_pts``.
    max_gap_in_sec:
        Optional maximum gap between the two samples bracketing a query.
        By default, use 3 times the satellite's median sample spacing.  A
        missing single epoch is usually safe; a long outage is not.  Set to
        a positive explicit value to override the inferred threshold.

    Notes
    -----
    Interpolation uses no cache.  This favors predictable, testable behavior
    over speed and eliminates stale-model reuse when the window changes.
    """

    def __init__(
        self,
        satid: str,
        *,
        tarray: Sequence,
        x: Sequence[float],
        y: Sequence[float],
        z: Sequence[float],
        clk: Sequence[float] | None = None,
        interval_in_sec: float = 1800.0,
        min_data_pts: int = 4,
        itype: str = "Linear",
        stencil_pts: int = 9,
        max_gap_in_sec: float | None = None,
    ) -> None:
        if itype not in _METHODS:
            raise ValueError(f"Unknown interpolation method {itype!r}; choose from {sorted(_METHODS)}")
        if not np.isfinite(interval_in_sec) or interval_in_sec <= 0:
            raise ValueError("interval_in_sec must be positive and finite")
        if isinstance(min_data_pts, bool) or int(min_data_pts) != min_data_pts or min_data_pts < 2:
            raise ValueError("min_data_pts must be an integer >= 2")
        if isinstance(stencil_pts, bool) or int(stencil_pts) != stencil_pts or stencil_pts < 2:
            raise ValueError("stencil_pts must be an integer >= 2")
        if max_gap_in_sec is not None and (not np.isfinite(max_gap_in_sec) or max_gap_in_sec <= 0):
            raise ValueError("max_gap_in_sec must be positive and finite")

        self._satellite = satid.upper()
        self._itype = itype
        self._minpts = int(min_data_pts)
        self._stencil_pts = max(int(stencil_pts), self._minpts)
        self._half_window_as = seconds_to_attoseconds(interval_in_sec)
        self._explicit_gap_as = (
            None if max_gap_in_sec is None else seconds_to_attoseconds(max_gap_in_sec)
        )
        self._t = tuple(to_attoseconds(v) for v in tarray)
        self._x = np.asarray(x, dtype=np.float64)
        self._y = np.asarray(y, dtype=np.float64)
        self._z = np.asarray(z, dtype=np.float64)
        self._c = np.full(len(self._t), np.nan) if clk is None else np.asarray(clk, dtype=np.float64)

        n = len(self._t)
        if n == 0 or any(len(series) != n for series in (self._x, self._y, self._z, self._c)):
            raise ValueError("Non-empty timestamps and equally sized x, y, z, clock arrays are required")
        if any(later <= earlier for earlier, later in zip(self._t, self._t[1:])):
            raise ValueError("Satellite epochs must be unique and strictly increasing")

        # The position series uses a shared mask, so X/Y/Z are always fitted
        # at the same epochs.  The clock mask is independent by design.
        xyz = np.column_stack((self._x, self._y, self._z))
        pos_valid = np.all(np.isfinite(xyz), axis=1)
        self._pos_times = tuple(t for t, valid in zip(self._t, pos_valid) if valid)
        self._positions = xyz[pos_valid]
        clock_valid = np.isfinite(self._c)
        self._clk_times = tuple(t for t, valid in zip(self._t, clock_valid) if valid)
        self._clocks = self._c[clock_valid]

        self._position_gap_as = self._gap_limit(self._pos_times)
        # The XYZ cadence is usually the best indication of the SP3 product
        # sampling rate.  Do not mistake two isolated clock corrections for a
        # valid, long-period clock time series.  Users can override the limit.
        self._clock_gap_as = (
            self._position_gap_as if self._position_gap_as is not None
            else self._gap_limit(self._clk_times)
        )

    def _gap_limit(self, epochs: tuple[int, ...]) -> int | None:
        """Determine a conservative maximum bracketing gap for this series."""
        if self._explicit_gap_as is not None:
            return self._explicit_gap_as
        if len(epochs) < 2:
            return None
        gaps = sorted(later - earlier for earlier, later in zip(epochs, epochs[1:]))
        # Use Python integer arithmetic: absolute attoseconds overflow int64,
        # and converting epoch values to float would erase fine resolution.
        median = gaps[len(gaps) // 2] if len(gaps) % 2 else (
            gaps[len(gaps) // 2 - 1] + gaps[len(gaps) // 2]) // 2
        # This automatic threshold is a *heuristic*, not a physical model.
        return 3 * median

    def _evaluate(
        self,
        query_as: int,
        epochs: tuple[int, ...],
        values: np.ndarray,
        *,
        gap_limit: int | None,
    ) -> np.ndarray:
        """Evaluate a scalar or vector time series at one epoch.

        Reject extrapolation, demand enough nearby valid samples, and reject
        a query which lies inside an unacceptably large observation gap.
        """
        if not epochs or query_as < epochs[0] or query_as > epochs[-1]:
            raise ValueError(f"No observations bracket the query for {self._satellite}")
        insertion = bisect_left(epochs, query_as)
        if insertion < len(epochs) and epochs[insertion] == query_as:
            # Exact available data should be returned exactly, without
            # spline roundoff and without unnecessary neighbor requirements.
            return np.asarray(values[insertion])
        if insertion == 0 or insertion == len(epochs):
            raise ValueError(f"Interpolation would extrapolate for {self._satellite}")
        if gap_limit is not None and int(epochs[insertion] - epochs[insertion - 1]) > gap_limit:
            raise ValueError(f"Observation gap too large for {self._satellite}")

        # Filter the window using integral attoseconds; never cast absolute
        # 2024-era epochs to float or numpy int64.
        indices = [i for i, epoch in enumerate(epochs)
                   if abs(epoch - query_as) <= self._half_window_as]
        ts = [epochs[i] for i in indices]
        ys = values[indices]
        if len(ts) < self._minpts:
            raise ValueError(
                f"Insufficient nearby samples for {self._satellite}: "
                f"found {len(ts)}, need {self._minpts}"
            )
        if ts[0] > query_as or ts[-1] < query_as:
            raise ValueError(f"Interpolation window does not bracket query for {self._satellite}")

        if self._itype in ("Polynomial", "Barycentric"):
            # Select *nearest* samples, then restore chronological order.
            # Using a bounded local stencil is safer numerically than fitting
            # 50+ nodes across a full day (Runge oscillations / roundoff).
            count = min(self._stencil_pts, len(ts))
            bracket = bisect_left(ts, query_as)
            # Always include the two bracketing observations.  Otherwise a
            # locally clustered set of nearest nodes might all fall on one
            # side of an irregularly sampled query, causing extrapolation.
            obligatory = {bracket - 1, bracket}
            remaining = [
                i for i in sorted(range(len(ts)), key=lambda i: abs(ts[i] - query_as))
                if i not in obligatory
            ]
            nearest = np.array(sorted(obligatory | set(remaining[:count - 2])), dtype=int)
            ts, ys = [ts[i] for i in nearest], ys[nearest]

        # Subtract Python integer attoseconds *before* converting to float.
        # This preserves tiny time offsets; using absolute Unix seconds as
        # floats would lose microseconds or worse in modern epochs.
        rel_seconds = np.array([(t - ts[0]) / _ATTO_PER_SECOND for t in ts], dtype=float)
        q_seconds = (query_as - ts[0]) / _ATTO_PER_SECOND
        if self._itype == "Linear":
            if ys.ndim == 1:
                return np.asarray(np.interp(q_seconds, rel_seconds, ys))
            return np.array([np.interp(q_seconds, rel_seconds, ys[:, col]) for col in range(ys.shape[1])])
        if self._itype == "CubicSpline":
            model = CubicSpline(rel_seconds, ys, axis=0)
        elif self._itype == "PchipInterpolator":
            model = PchipInterpolator(rel_seconds, ys, axis=0)
        else:
            model = BarycentricInterpolator(rel_seconds, ys, axis=0)
        return np.asarray(model(q_seconds))

    def sat_at(self, t: Any) -> tuple[float, float, float, float]:
        """Return ``(x, y, z, clock)`` in metres and seconds.

        An unusable position raises ``ValueError``.  An unavailable clock
        returns ``NaN`` without discarding a perfectly valid position.
        """
        query_as = to_attoseconds(t)
        pos = self._evaluate(
            query_as, self._pos_times, self._positions, gap_limit=self._position_gap_as
        )
        try:
            clock = float(
                self._evaluate(query_as, self._clk_times, self._clocks, gap_limit=self._clock_gap_as)
            )
        except ValueError:
            # No clock, gap, or too few valid clock samples does not imply
            # anything about the quality of the independent XYZ orbit.
            clock = float("nan")
        return float(pos[0]), float(pos[1]), float(pos[2]), clock


class Sp3Interpolator:
    """Provide per-satellite interpolation from a shared SP3 epoch mapping.

    Input ``data`` retains the original project convention:
    ``{epoch: {satellite_id: [x_m, y_m, z_m, clock_s, flags]}}``.
    """

    def __init__(
        self,
        data: Mapping,
        *,
        interval_in_sec: float = 1800.0,
        min_data_pts: int = 4,
        itype: str = "Linear",
        stencil_pts: int = 9,
        max_gap_in_sec: float | None = None,
        exclude_missing_clock_values: bool = False,
        exclude_flag_events: Sequence[str] = (),
        time_system: str = "UNKNOWN",
    ) -> None:
        if not data:
            raise ValueError("SP3 contains no epochs")
        self._time_system = time_system
        # Preserve the common time base, even where only one of the satellites
        # has an observation.  This is important for mixed GNSS SP3 products.
        self._t = tuple(sorted({to_attoseconds(t) for t in data}))
        index_by_time = {epoch: idx for idx, epoch in enumerate(self._t)}
        satellites = {str(sat).upper() for records in data.values() for sat in records}
        if not satellites:
            raise ValueError("No satellite position records for the requested GNSS systems")
        self._interpolators: dict[str, OrbitInterpolator] = {}
        for sat in sorted(satellites):
            n = len(self._t)
            xyz = np.full((n, 3), np.nan)
            clocks = np.full(n, np.nan)
            for epoch, records in data.items():
                # ``records`` can use upper or lower-case satellite IDs.
                entry = next((value for key, value in records.items() if key.upper() == sat), None)
                if entry is None:
                    continue
                if len(entry) != 5:
                    raise ValueError(f"Invalid record for {sat}: expected [x, y, z, clock, flags]")
                x, y, z, clk, flags = entry
                if flag_is_on(str(flags), exclude_flag_events):
                    continue
                if exclude_missing_clock_values and (clk is None or not np.isfinite(float(clk))):
                    # This option intentionally drops XYZ alongside the clock.
                    continue
                i = index_by_time[to_attoseconds(epoch)]
                xyz[i] = (float(x), float(y), float(z))
                if clk is not None:
                    clocks[i] = float(clk)
            self._interpolators[sat] = OrbitInterpolator(
                satid=sat,
                tarray=[Epoch(t) for t in self._t],
                x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], clk=clocks,
                interval_in_sec=interval_in_sec,
                min_data_pts=min_data_pts,
                itype=itype,
                stencil_pts=stencil_pts,
                max_gap_in_sec=max_gap_in_sec,
            )
        self._sat_ids = sorted(self._interpolators)

    def sat_at(self, satid: str, t: Any) -> tuple[float, float, float, float]:
        """Interpolate a satellite identified by its SP3 three-character ID."""
        satid = satid.upper()
        if satid not in self._interpolators:
            raise KeyError(f"Satellite {satid!r} not found; available: {self._sat_ids}")
        return self._interpolators[satid].sat_at(t)

    @property
    def time_base(self) -> np.ndarray:
        """Shared, sorted SP3 epochs in ``numpy.datetime64[ns]`` format."""
        # The historical property exposes numpy datetime64[ns]. Users needing
        # sub-nanosecond epochs should use time_base_attoseconds instead.
        return np.array([np.datetime64(t // 10**9, 'ns') for t in self._t])

    @property
    def time_base_attoseconds(self) -> tuple[int, ...]:
        """Return the full-precision, normalized epoch labels (Python ints)."""
        return self._t

    @property
    def satellites(self) -> list[str]:
        """Satellite identifiers in the selected constellations."""
        return self._sat_ids.copy()

    def sat_list(self) -> list[str]:
        """Compatibility helper used in the original tutorial."""
        return self.satellites

    def time_system(self) -> str:
        """Return the SP3 header time system; no time-scale conversion occurs."""
        return self._time_system

    @classmethod
    def from_sp3(
        cls,
        filename,
        satsys: Sequence[str] | str,
        *,
        strict_epoch_count: bool = False,
        **kwargs,
    ) -> "Sp3Interpolator":
        """Read SP3-c/d and interpolate the satellite systems requested.

        ``strict_epoch_count=False`` preserves the original notebook API and
        tolerates a stale declared header count with an explicit warning.
        Choose True when validating unedited scientific orbit products.  This
        option belongs to the SP3 reader, not the interpolation algorithm.
        """
        from dsoclasses.orbits.sp3c import Sp3

        source = Sp3(filename, strict_epoch_count=strict_epoch_count)
        if "time_system" in kwargs:
            raise TypeError("from_sp3 determines time_system from the SP3 header")
        data = source.get_system_pos(satsys, toSI=True)
        return cls(data, time_system=source.time_sys, **kwargs)
