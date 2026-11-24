"""Small, shared building blocks for the two cycle-slip teaching algorithms.

Code ranges are metres, carrier phases are RINEX cycles, frequencies are Hz,
and times are elapsed seconds.  No orbit, atmospheric model or station clock
is needed: geometry and non-dispersive terms cancel in the combinations.

This module deliberately does not silently 'fix' questionable observations.
Repairs are accepted only when two independent combinations identify one
integer pair; otherwise the ambiguity arc is broken and an event is reported.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Sequence

import numpy as np

from dsoclasses.gnss.systems import C, GPS_L1_FREQ, GPS_L2_FREQ
from dsoclasses.time.pyattotime import to_attoseconds


@dataclass(frozen=True)
class DualFrequencySeries:
    """One uninterrupted satellite's dual-frequency observations.

    ``t_s`` is elapsed time, ``p*_m`` are code pseudoranges, and ``l*_cycles``
    are carrier phase measurements. The input is never mutated by a detector.
    """

    t_s: Sequence[float]
    p1_m: Sequence[float]
    p2_m: Sequence[float]
    l1_cycles: Sequence[float]
    l2_cycles: Sequence[float]
    f1_hz: float
    f2_hz: float
    satellite: str = "G00"
    signals: tuple[str, str, str, str] = ()

    def __post_init__(self):
        for key in ("t_s", "p1_m", "p2_m", "l1_cycles", "l2_cycles"):
            value = np.asarray(getattr(self, key), dtype=float)
            if value.ndim != 1 or not np.all(np.isfinite(value)):
                raise ValueError(f"{key} must be a finite one-dimensional array")
            value = value.copy()
            value.setflags(write=False)
            object.__setattr__(self, key, value)
        n = len(self.t_s)
        if n < 3 or any(len(getattr(self, key)) != n for key in
                        ("p1_m", "p2_m", "l1_cycles", "l2_cycles")):
            raise ValueError("At least three aligned dual-frequency epochs are required")
        if np.any(np.diff(self.t_s) <= 0):
            raise ValueError("Epochs must increase strictly")
        if not (np.isfinite(self.f1_hz) and np.isfinite(self.f2_hz)
                and self.f1_hz > self.f2_hz > 1e8):
            raise ValueError("Frequencies must be distinct GNSS values in Hz "
                             "(f1 > f2 > 1e8); do not pass MHz")

    @property
    def wavelength1_m(self) -> float:
        return C / self.f1_hz

    @property
    def wavelength2_m(self) -> float:
        return C / self.f2_hz

    @property
    def wavelength_wide_m(self) -> float:
        return C / (self.f1_hz - self.f2_hz)

    @property
    def gf_metres_per_tecu(self) -> float:
        """Conversion factor of the geometry-free phase observable [m/TECU]."""
        return (40.3e16 / self.f1_hz**2 *
                ((self.f1_hz / self.f2_hz)**2 - 1))

    def with_changes(self, *, slips=None, code_outliers=None, phase_outliers=None):
        """Teaching helper: add known *cumulative* cycle slips or isolated errors.

        ``slips={index: (delta_L1_cycles, delta_L2_cycles)}``; a slip changes
        every subsequent phase until repaired. A code/phase outlier affects
        only its named index. No passed-in data is overwritten.
        """
        l1 = self.l1_cycles.copy()
        l2 = self.l2_cycles.copy()
        p1 = self.p1_m.copy()
        p2 = self.p2_m.copy()
        for i, pair in (slips or {}).items():
            if not 0 <= i < len(l1) or len(pair) != 2 or any(int(v) != v for v in pair):
                raise ValueError("Invalid slip index or non-integer cycle pair")
            l1[i:] += pair[0]
            l2[i:] += pair[1]
        for i, pair in (code_outliers or {}).items():
            if not 0 <= i < len(l1) or len(pair) != 2:
                raise ValueError("Invalid code-outlier index or pair")
            p1[i] += pair[0]
            p2[i] += pair[1]
        for i, pair in (phase_outliers or {}).items():
            if not 0 <= i < len(l1) or len(pair) != 2:
                raise ValueError("Invalid phase-outlier index or pair")
            l1[i] += pair[0]
            l2[i] += pair[1]
        return replace(self, p1_m=p1, p2_m=p2, l1_cycles=l1, l2_cycles=l2)


@dataclass(frozen=True)
class CycleSlipEvent:
    index: int
    time_s: float
    kind: str  # 'repaired', 'code_outlier', 'phase_outlier', 'unresolved', 'gap'
    slip_l1: int | None = None
    slip_l2: int | None = None
    mw_step_cycles: float = float("nan")
    gf_step_m: float = float("nan")


@dataclass
class CycleSlipResult:
    method: str
    original: DualFrequencySeries
    mw_cycles: np.ndarray
    gf_m: np.ndarray
    mw_prediction: np.ndarray
    gf_prediction_m: np.ndarray
    mw_innovation: np.ndarray
    gf_innovation_m: np.ndarray
    repaired_l1_cycles: np.ndarray
    repaired_l2_cycles: np.ndarray
    status: list[str]
    arc_id: np.ndarray
    events: list[CycleSlipEvent] = field(default_factory=list)
    tecr_observed: np.ndarray | None = None
    tecr_predicted: np.ndarray | None = None

    @property
    def slips(self) -> list[CycleSlipEvent]:
        return [e for e in self.events if e.kind == "repaired"]


def combinations(data: DualFrequencySeries) -> tuple[np.ndarray, np.ndarray]:
    """Return Melbourne-Wuebbena [wide-lane cycles], geometry-free [metres].

    For positive phase jumps (dN1,dN2), the MW jump is dN1-dN2 and the GF
    jump is lambda1*dN1-lambda2*dN2. No satellite/receiver biases are needed
    for *changes* within a fixed-signal continuous arc.
    """
    wide_code = ((data.f1_hz * data.p1_m + data.f2_hz * data.p2_m)
                 / (data.f1_hz + data.f2_hz))
    mw = data.l1_cycles - data.l2_cycles - wide_code / data.wavelength_wide_m
    gf = data.wavelength1_m * data.l1_cycles - data.wavelength2_m * data.l2_cycles
    return mw, gf


def _local_stats(values, *, sigma_floor):
    """Median and robust scale; use a physical floor while history is short."""
    a = np.asarray(values, dtype=float)
    if len(a) == 0:
        raise ValueError("Cannot estimate a level without observations")
    centre = float(np.median(a))
    sigma = (1.4826 * float(np.median(np.abs(a - centre)))
             if len(a) >= 4 else 0.)
    return centre, max(sigma_floor, sigma)


def _integer_pair(mw_step, gf_step, data, *, sigma_mw, sigma_gf,
                  max_cycles=100):
    """Search integer pairs satisfying *both* independent jump equations.

    Reject a plausible-looking but non-unique integer choice. This matters
    especially with noisy code data or isolated pseudorange outliers.
    """
    if not all(np.isfinite((mw_step, gf_step, sigma_mw, sigma_gf))):
        return None
    l1, l2 = data.wavelength1_m, data.wavelength2_m
    candidates = []
    for d in range(round(mw_step) - 2, round(mw_step) + 3):
        n2_float = (gf_step - l1 * d) / (l1 - l2)
        for n2 in range(round(n2_float) - 3, round(n2_float) + 4):
            n1 = n2 + d
            if max(abs(n1), abs(n2)) > max_cycles or (n1 == n2 == 0):
                continue
            r_mw = mw_step - d
            r_gf = gf_step - (l1 * n1 - l2 * n2)
            score = (r_mw / sigma_mw)**2 + (r_gf / sigma_gf)**2
            candidates.append((score, n1, n2, r_mw, r_gf))
    candidates.sort(key=lambda row: row[0])
    if not candidates:
        return None
    score, n1, n2, residual_mw, residual_gf = candidates[0]
    # A minimum alone does NOT prove an integer ambiguity is identifiable.
    if (abs(residual_mw) > max(0.35, 3 * sigma_mw)
            or abs(residual_gf) > max(0.018, 3 * sigma_gf)):
        return None
    if len(candidates) > 1 and candidates[1][0] - score < 9.0:
        return None
    return n1, n2


def _persistent_jump(i, data, mw, gf, *, mw_prediction, gf_prediction,
                     mw_step, gf_step, mw_limit, gf_limit, cumulative1,
                     cumulative2, max_gap_s):
    """Offline one-epoch check: distinguish an isolated spike from a step.

    Returns None if no next epoch exists; never claim a last-epoch repair.
    """
    j = i + 1
    if j >= len(data.t_s) or data.t_s[j] - data.t_s[i] > max_gap_s:
        return None
    next_mw = mw[j] - (cumulative1 - cumulative2)
    next_gf = (gf[j] - data.wavelength1_m * cumulative1
               + data.wavelength2_m * cumulative2)
    mw_persists = (abs((next_mw - mw_prediction) - mw_step)
                   < max(0.25, mw_limit * .65))
    # The ionosphere can drift: only require that a GF *step* is still there.
    gf_persists = (abs((next_gf - gf_prediction) - gf_step)
                   < max(0.018, gf_limit * .65))
    return bool((mw_persists if abs(mw_step) > mw_limit else True) and (
        gf_persists if abs(gf_step) > gf_limit else True))


def teaching_gps_series(n=420, *, seed=2024, code_sigma_m=.055,
                        phase_sigma_m=.002) -> DualFrequencySeries:
    """Deterministic clean 1 Hz GPS arc; students can inject known errors.

    Includes realistic dual-frequency ionospheric dispersion, independent
    code/carrier noise, kinematic range and fixed carrier ambiguities.
    """
    if not isinstance(n, int) or n < 90:
        raise ValueError("Use at least 90 epochs for the teaching example")
    random = np.random.default_rng(seed)
    f1, f2 = GPS_L1_FREQ * 1e6, GPS_L2_FREQ * 1e6
    t = np.arange(n, dtype=float)
    distance = 23_000_000. + 85. * t + .014 * t**2
    tec = 16. + .005 * t + .14 * np.sin(t / 110.)
    i1 = 40.3e16 * tec / f1**2
    i2 = 40.3e16 * tec / f2**2
    p1 = distance + i1 + random.normal(0, code_sigma_m, n)
    p2 = distance + i2 + random.normal(0, code_sigma_m, n)
    l1 = ((distance - i1 + random.normal(0, phase_sigma_m, n))
          / (C / f1)) + 140
    l2 = ((distance - i2 + random.normal(0, phase_sigma_m, n))
          / (C / f2)) + 93
    return DualFrequencySeries(t, p1, p2, l1, l2, f1, f2,
                               satellite="G21", signals=("C1W", "C2W", "L1W", "L2W"))


_GPS_L1 = (("C1W", "L1W"), ("C1C", "L1C"), ("C1X", "L1X"), ("C1L", "L1L"))
_GPS_L2 = (("C2W", "L2W"), ("C2P", "L2P"), ("C2X", "L2X"),
           ("C2L", "L2L"), ("C2S", "L2S"), ("C2C", "L2C"), ("C2D", "L2D"))


def gps_series_from_rinex(filename, satellite="G21", *, max_samples=400) -> DualFrequencySeries:
    """Read a single GPS satellite with one *fixed* code/phase signal pair.

    A tracking-code switch can look exactly like a slip. Choose the pair
    observed most often, keep it fixed, and leave actual missing epochs as
    gaps. This is an exploratory source for the 30 s classroom RINEX file.
    """
    from dsoclasses.rinex.gnss.rinex import GnssRinex

    if not satellite.startswith("G") or len(satellite) != 3:
        raise ValueError("This introductory RINEX helper supports GPS PRNs")
    reader = GnssRinex(filename)
    rows = [(block.t(), block.dct[satellite]) for block in reader
            if satellite in block.dct]
    if not rows:
        raise LookupError(f"No RINEX observations for {satellite}")
    options = [(a, b) for a in _GPS_L1 for b in _GPS_L2]
    def complete(record, pair):
        for code in (*pair[0], *pair[1]):
            item = record.get(code)
            if not isinstance(item, dict):
                return False
            try:
                if not np.isfinite(float(item['value'])):
                    return False
            except (TypeError, ValueError, KeyError):
                return False
        return True
    chosen = max(options, key=lambda pair: sum(complete(rec, pair) for _, rec in rows))
    kept = [(t, rec) for t, rec in rows if complete(rec, chosen)]
    if max_samples is not None:
        kept = kept[:int(max_samples)]
    if len(kept) < 3:
        raise LookupError(f"Not enough complete dual-frequency observations for {satellite}")
    origin = to_attoseconds(kept[0][0])
    c1, l1 = chosen[0]
    c2, l2 = chosen[1]
    return DualFrequencySeries(
        t_s=[(to_attoseconds(t) - origin) / 1e18 for t, _ in kept],
        p1_m=[rec[c1]['value'] for _, rec in kept],
        p2_m=[rec[c2]['value'] for _, rec in kept],
        l1_cycles=[rec[l1]['value'] for _, rec in kept],
        l2_cycles=[rec[l2]['value'] for _, rec in kept],
        f1_hz=GPS_L1_FREQ * 1e6, f2_hz=GPS_L2_FREQ * 1e6,
        satellite=satellite, signals=(c1, c2, l1, l2),
    )
