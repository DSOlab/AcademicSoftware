# interpolator.py
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator, BarycentricInterpolator

# ---------------------------------------------------------------------------
# External time utilities
#   - to_attoseconds(t): normalize attotime/datetime/np.datetime64 -> int attoseconds
#   - ATTO_PER_SEC: int, 10**18
#   - ATTO_PER_NS:  int, 10**9  (attoseconds per nanosecond)
# ---------------------------------------------------------------------------
from dsoclasses.time.pyattotime import (
    to_attoseconds,
    _ATTO_PER_SEC,
    _ATTO_PER_NS,
)

# =============================================================================
# Helpers (non-time)
# =============================================================================


def flag_is_on(flag_str: str, selected: Sequence[str]) -> bool:
    """Return True if any of the selected one-letter SP3 flags are present."""
    if not flag_str or not selected:
        return False
    s = set(str(flag_str).replace(" ", ""))
    return any(ch in s for ch in selected)


def _as_np_datetime64_ns(seq_dt) -> np.ndarray:
    """Convert an iterable of datetimes to a sorted, unique numpy datetime64[ns] array."""
    uniq = sorted(set(seq_dt))
    return np.array(uniq, dtype="datetime64[ns]")


def _chebyshev_lobatto_indices(m: int, n: int) -> np.ndarray:
    """
    Return n integer indices in [0, m-1] approximating Chebyshev-Lobatto nodes
    clustered at both ends of the interval.
    """
    if n <= 1:
        return np.array([0], dtype=int)
    if n >= m:
        return np.arange(m, dtype=int)

    k = np.arange(n, dtype=float)
    x = np.cos(np.pi * k / (n - 1))  # [-1, 1]
    u = (x + 1.0) * 0.5  # [0, 1]
    idxf = u * (m - 1)  # [0, m-1]
    idx = np.rint(idxf).astype(int)

    # ensure strictly increasing, unique indices
    seen = set()
    out: List[int] = []
    for i in idx:
        if i not in seen:
            out.append(i)
            seen.add(i)
    if len(out) < n:
        # fill with nearest unused endpoints first, then inner points
        candidates = list(range(m))
        dist_to_end = np.minimum(candidates, np.abs(np.array(candidates) - (m - 1)))
        order = np.argsort(dist_to_end)
        for c in np.array(candidates)[order]:
            if c not in seen:
                out.append(c)
                seen.add(c)
                if len(out) == n:
                    break
    return np.array(sorted(out), dtype=int)


# =============================================================================
# OrbitInterpolator: per-satellite on shared time base
# =============================================================================


class OrbitInterpolator:
    """
    Interpolates a single satellite's position (and optionally clock) on a
    shared time base (numpy datetime64[ns]), with NaNs where data is missing.

    Supported interpolation types (per-axis):
      - "Linear" (alias: "Polynomial"): np.interp over a masked window
      - "CubicSpline": scipy.interpolate.CubicSpline
      - "PchipInterpolator": scipy.interpolate.PchipInterpolator
      - "Barycentric": scipy.interpolate.BarycentricInterpolator
         * nodes in the window are chosen with Chebyshev-Lobatto–like
           index selection (clustered near the endpoints).
    """

    def __init__(
        self,
        satid: str,
        *,
        tarray: np.ndarray,  # datetime64[ns] shared time base
        x: Sequence[float],
        y: Sequence[float],
        z: Sequence[float],
        clk: Optional[Sequence[float]] = None,
        interval_in_sec: float = 1800.0,  # half-width of window used to build local model
        min_data_pts: int = 4,
        itype: str = "Linear",
        stencil_pts: int = 9,  # used by Barycentric; desired nodes in window
    ) -> None:
        self._satellite = satid
        self._t: np.ndarray = np.asarray(tarray, dtype="datetime64[ns]")
        self._x = np.asarray(x, dtype=float)
        self._y = np.asarray(y, dtype=float)
        self._z = np.asarray(z, dtype=float)
        self._c = None if clk is None else np.asarray(clk, dtype=float)

        self._itype = "Linear" if itype == "Polynomial" else itype
        self._dsec = float(interval_in_sec)
        self._minpts = int(min_data_pts)
        self._stencil_pts = int(stencil_pts)

        # caches
        self._last_key: Optional[Tuple[int, int, str, Tuple[int, ...]]] = None
        self._xspl = None
        self._yspl = None
        self._zspl = None
        self._cspl = None
        self._grid_base_ns: Optional[int] = None

        # sanity
        n = len(self._t)
        if not (len(self._x) == len(self._y) == len(self._z) == n):
            raise ValueError("x, y, z arrays must match the shared time base length")

    # ---------------------------------------------------------------------

    def _window_indices(self, t: Any) -> Tuple[int, int]:
        """Return indices [start, stop) for epochs within ±self._dsec of t."""
        # Convert query to attoseconds, then to nanoseconds to match the time base
        tq_ns = to_attoseconds(t) // _ATTO_PER_NS

        # Nanoseconds since epoch for the time base
        ns = self._t.astype("datetime64[ns]").astype("int64")

        # Window half-width in nanoseconds
        d_ns = int(self._dsec * 1_000_000_000)

        # find first index where (tq - ti) <= d
        start = None
        for i, ti in enumerate(ns):
            if (tq_ns - ti) <= d_ns:
                start = i
                break
        if start is None:
            return -1, -1

        # find stop index where (ti - tq) > d
        stop = None
        for j, ti in enumerate(ns[start:]):
            if (ti - tq_ns) > d_ns:
                stop = start + j
                break
        if stop is None:
            stop = len(self._t)
        return start, stop

    def _slice_seconds(self, start: int, stop: int) -> np.ndarray:
        """
        Convert self._t[start:stop] to float seconds relative to the first node.
        Work purely in nanoseconds to avoid overflow, then scale.
        """
        ts = (
            self._t[start:stop].astype("datetime64[ns]").astype("int64")
        )  # ns since epoch
        base_ns = int(ts[0])
        rel_ns = (ts - base_ns).astype(np.float64)
        return rel_ns / 1_000_000_000.0

    def _masked_window(
        self, start: int, stop: int
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
        """
        Build masked per-axis arrays for the [start, stop) window, dropping rows
        with NaNs in any of X/Y/Z. Returns (tx_sec, xx, yy, zz, cc_or_None).
        """
        tx_sec = self._slice_seconds(start, stop)
        xw = self._x[start:stop]
        yw = self._y[start:stop]
        zw = self._z[start:stop]
        mask = ~(np.isnan(xw) | np.isnan(yw) | np.isnan(zw))
        tx = tx_sec[mask]
        xx = xw[mask]
        yy = yw[mask]
        zz = zw[mask]
        if self._c is not None:
            cw = self._c[start:stop]
            cc = cw[mask]
        else:
            cc = None
        return tx, xx, yy, zz, cc

    def _barycentric_subset(
        self, start: int, stop: int, desired: int
    ) -> Tuple[
        np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray], int
    ]:
        """
        Choose ~Chebyshev-Lobatto indices inside [start, stop), mask NaNs, and return
        the sub-selected arrays plus the base_ns used for evaluation.
        """
        m = max(0, stop - start)
        if m < 2:
            return np.array([]), np.array([]), np.array([]), np.array([]), None, 0

        # initial indices clustered near ends
        rel_idx = _chebyshev_lobatto_indices(m, min(desired, m))
        idx = start + rel_idx

        # mask NaNs
        xw = self._x[idx]
        yw = self._y[idx]
        zw = self._z[idx]
        valid = ~(np.isnan(xw) | np.isnan(yw) | np.isnan(zw))
        idx = idx[valid]
        if idx.size < self._minpts:
            # fallback: expand to full masked window
            tx, xx, yy, zz, cc = self._masked_window(start, stop)
            base_ns = int(self._t[start].astype("datetime64[ns]").astype("int64"))
            return tx, xx, yy, zz, cc, base_ns

        # seconds relative to first selected node
        ts = self._t[idx].astype("datetime64[ns]").astype("int64")
        base_ns = int(ts[0])
        rel_ns = (ts - base_ns).astype(np.float64)
        tx = rel_ns / 1_000_000_000.0

        xx = self._x[idx]
        yy = self._y[idx]
        zz = self._z[idx]
        cc = self._c[idx] if self._c is not None else None
        return tx, xx, yy, zz, cc, base_ns

    def _build_model(self, start: int, stop: int):
        """
        Build/cache the local interpolation model depending on self._itype.
        For Barycentric, subselect nodes clustered at endpoints.
        """
        if self._itype == "Linear":
            self._xspl = self._yspl = self._zspl = self._cspl = None
            self._grid_base_ns = int(
                self._t[start].astype("datetime64[ns]").astype("int64")
            )
            key = (start, stop, "Linear", ())
            self._last_key = key
            return

        if self._itype == "Barycentric":
            tx, xx, yy, zz, cc, base_ns = self._barycentric_subset(
                start, stop, self._stencil_pts
            )
            if tx.size < self._minpts:
                raise RuntimeError(
                    f"Not enough valid data for Barycentric interpolation: have {tx.size}, need ≥ {self._minpts}"
                )
            self._xspl = BarycentricInterpolator(tx, xx)
            self._yspl = BarycentricInterpolator(tx, yy)
            self._zspl = BarycentricInterpolator(tx, zz)
            self._cspl = (
                BarycentricInterpolator(tx, cc)
                if cc is not None and np.any(~np.isnan(cc))
                else None
            )
            self._grid_base_ns = base_ns
            key = (start, stop, "Barycentric", tuple(tx.shape))
            self._last_key = key
            return

        # CubicSpline / Pchip: use masked full window
        tx, xx, yy, zz, cc = self._masked_window(start, stop)
        if tx.size < self._minpts:
            raise RuntimeError(
                f"Not enough valid data in stencil: have {tx.size}, need ≥ {self._minpts}"
            )
        if self._itype == "CubicSpline":
            self._xspl = CubicSpline(tx, xx)
            self._yspl = CubicSpline(tx, yy)
            self._zspl = CubicSpline(tx, zz)
            self._cspl = (
                CubicSpline(tx, cc)
                if cc is not None and np.any(~np.isnan(cc))
                else None
            )
        elif self._itype == "PchipInterpolator":
            self._xspl = PchipInterpolator(tx, xx)
            self._yspl = PchipInterpolator(tx, yy)
            self._zspl = PchipInterpolator(tx, zz)
            self._cspl = (
                PchipInterpolator(tx, cc)
                if cc is not None and np.any(~np.isnan(cc))
                else None
            )
        else:
            raise ValueError(f"Unknown interpolation type '{self._itype}'")
        self._grid_base_ns = int(
            self._t[start].astype("datetime64[ns]").astype("int64")
        )
        key = (start, stop, self._itype, ())
        self._last_key = key

    # ---------------------------------------------------------------------

    def sat_at(self, t: Any) -> Tuple[float, float, float, float]:
        """Interpolate satellite at epoch 't' (attotime-like or datetime)."""
        start, stop = self._window_indices(t)
        if start < 0 or stop < 0 or (stop - start) < self._minpts:
            raise RuntimeError(
                f"ERROR Cannot find suitable interval for interpolating satellite orbit at {t}"
            )

        key = (start, stop, self._itype, self._last_key[3] if self._last_key else ())
        if self._last_key != key:
            self._build_model(start, stop)

        # evaluate at query time (seconds relative to base_ns)
        tq_ns = to_attoseconds(t) // _ATTO_PER_NS
        tq_sec = (tq_ns - self._grid_base_ns) / 1_000_000_000.0

        if self._itype == "Linear":
            tx, xx, yy, zz, cc = self._masked_window(start, stop)
            if tx.size < 2:
                raise RuntimeError(
                    f"Not enough valid data for linear interpolation near {t}"
                )
            x = np.interp(tq_sec, tx, xx)
            y = np.interp(tq_sec, tx, yy)
            z = np.interp(tq_sec, tx, zz)
            c = (
                np.interp(tq_sec, tx, cc)
                if cc is not None and np.any(~np.isnan(cc))
                else np.nan
            )
            return float(x), float(y), float(z), float(c)

        x = float(self._xspl(tq_sec))
        y = float(self._yspl(tq_sec))
        z = float(self._zspl(tq_sec))
        c = float(self._cspl(tq_sec)) if self._cspl is not None else np.nan
        return x, y, z, c


# =============================================================================
# Sp3Interpolator: orchestrates per-satellite interpolators over shared time base
# =============================================================================


class Sp3Interpolator:
    """
    Build per-satellite OrbitInterpolators on a SINGLE shared datetime64[ns] time base.

    Expected input 'data' format (as produced by Sp3.get_system_pos(...)):
      data: Dict[datetime, Dict[str, List]] where
        data[epoch][satid] = [x, y, z, clk, flag]
    """

    def __init__(
        self,
        data: Dict,  # Dict[datetime, Dict[satid, [x,y,z,clk,flag]]]
        *,
        interval_in_sec: float = 1800.0,
        min_data_pts: int = 4,
        itype: str = "Linear",
        stencil_pts: int = 9,  # for Barycentric
        exclude_missing_clock_values: bool = False,
        exclude_flag_events: Sequence[str] = (),
    ) -> None:
        # Build shared time base
        self._t: np.ndarray = _as_np_datetime64_ns(data.keys())

        # Satellite list
        sat_ids = set()
        for d in data.values():
            sat_ids.update(d.keys())
        self._sat_ids = sorted(sat_ids)

        # Map time -> index (FIXED: enumerate to bind both idx and t)
        self._gkeys = {
            int(np.datetime64(t, "ns").astype("int64")): idx
            for idx, t in enumerate(self._t)
        }
        n = len(self._t)

        self._interpolators: Dict[str, OrbitInterpolator] = {}

        for sat in self._sat_ids:
            X = np.full(n, np.nan)
            Y = np.full(n, np.nan)
            Z = np.full(n, np.nan)
            C = np.full(n, np.nan)

            for k_dt, sats in data.items():
                if sat not in sats:
                    continue
                idx = self._gkeys[int(np.datetime64(k_dt, "ns").astype("int64"))]
                x, y, z, clk, flag = sats[sat]
                if exclude_flag_events and flag_is_on(str(flag), exclude_flag_events):
                    continue
                X[idx] = float(x)
                Y[idx] = float(y)
                Z[idx] = float(z)
                if exclude_missing_clock_values and (clk is None or np.isnan(clk)):
                    pass
                else:
                    C[idx] = float(clk) if clk is not None else np.nan

            self._interpolators[sat] = OrbitInterpolator(
                satid=sat,
                tarray=self._t,
                x=X,
                y=Y,
                z=Z,
                clk=C,
                interval_in_sec=interval_in_sec,
                min_data_pts=min_data_pts,
                itype=itype,
                stencil_pts=stencil_pts,
            )

    def sat_at(self, satid: str, t: Any) -> Tuple[float, float, float, float]:
        """Interpolate 'satid' at epoch 't' (attotime-like or datetime)."""
        if satid not in self._interpolators:
            raise KeyError(f"Satellite '{satid}' not found in SP3 dataset.")
        return self._interpolators[satid].sat_at(t)

    @property
    def time_base(self) -> np.ndarray:
        return self._t

    @property
    def satellites(self) -> List[str]:
        return list(self._sat_ids)

    # -------- convenience constructors --------
    @classmethod
    def from_sp3(
        cls,
        filename: str,
        satsys: Sequence[str],
        **kwargs,
    ) -> "Sp3Interpolator":
        """
        Convenience: build directly from an SP3 filename.
        Equivalent to: Sp3(filename).get_system_pos(satsys)
        """
        from dsoclasses.orbits.sp3c import Sp3  # lazy import

        sp3 = Sp3(filename)
        data = sp3.get_system_pos(list(satsys), toSI=True)
        return cls(data, **kwargs)
