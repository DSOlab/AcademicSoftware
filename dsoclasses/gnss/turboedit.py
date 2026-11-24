"""Teaching-scale, offline TurboEdit-style dual-frequency cycle-slip editing.

Blewitt (1990), doi:10.1029/GL017i003p00199. This is a transparent
implementation of the Melbourne--Wuebbena / geometry-free principle, not a
line-for-line reproduction of the GIPSY software. The one-epoch look-ahead
used to distinguish isolated bad data is an explicitly *offline* addition.
"""

from __future__ import annotations

import numpy as np

from .cycle_slip_common import (
    CycleSlipEvent, CycleSlipResult, DualFrequencySeries, _integer_pair,
    _local_stats, _persistent_jump, combinations,
)


def turboedit(data: DualFrequencySeries, *, history: int = 20,
              max_gap_s: float = 60.0) -> CycleSlipResult:
    """Detect and repair confident two-frequency integer phase jumps.

    Input: time [s], pseudoranges [m], phases [cycles], frequencies [Hz].
    A valid repair is a *pair* of integer jumps verified by MW and GF.
    Isolated outliers are marked, unresolved jumps break the phase arc.
    No correction is applied at a gap, or without sufficient clean history.

    Thresholds are practical teaching choices, not GIPSY tuning constants.
    """
    if history < 6 or max_gap_s <= 0:
        raise ValueError("history must be >= 6 and max_gap_s positive")
    mw, gf = combinations(data)
    n = len(mw)
    pred_mw = np.full(n, np.nan)
    pred_gf = np.full(n, np.nan)
    delta_mw = np.full(n, np.nan)
    delta_gf = np.full(n, np.nan)
    phase1 = np.full(n, np.nan)
    phase2 = np.full(n, np.nan)
    statuses = ["warmup"] * n
    arcs = np.zeros(n, dtype=int)
    events = []
    good = []  # (epoch, corrected MW, corrected GF), never outliers
    cumulative1 = cumulative2 = 0
    arc = 0

    for i, t in enumerate(data.t_s):
        if i and t - data.t_s[i-1] > max_gap_s:
            events.append(CycleSlipEvent(i, float(t), "gap"))
            arc += 1
            good.clear()
            cumulative1 = cumulative2 = 0
        arcs[i] = arc
        current_mw = mw[i] - (cumulative1 - cumulative2)
        current_gf = (gf[i] - data.wavelength1_m * cumulative1
                      + data.wavelength2_m * cumulative2)

        if len(good) < 6:
            phase1[i] = data.l1_cycles[i] - cumulative1
            phase2[i] = data.l2_cycles[i] - cumulative2
            good.append((t, current_mw, current_gf))
            statuses[i] = "warmup"
            continue

        recent = good[-history:]
        centre_mw, sigma_mw = _local_stats([v[1] for v in recent],
                                           sigma_floor=.055)
        # Linear extrapolation of the local ionospheric trend; centering
        # time avoids poorly conditioned fits at large epoch labels.
        fit = good[-min(len(good), 12):]
        dt = np.asarray([row[0]-fit[-1][0] for row in fit])
        g = np.asarray([row[2] for row in fit])
        slope = float(np.polyfit(dt, g, 1)[0])
        expected_gf = float(g[-1] + slope * (t - fit[-1][0]))
        # Physical phase-noise floor plus robust prior detrended scatter.
        scatter = g - (slope * dt + float(np.median(g-slope*dt)))
        _, sigma_gf = _local_stats(scatter, sigma_floor=.003)
        mw_limit = max(.45, 5 * sigma_mw)
        gf_limit = max(.035, 5 * sigma_gf)

        pred_mw[i], pred_gf[i] = centre_mw, expected_gf
        delta_mw[i] = current_mw - centre_mw
        delta_gf[i] = current_gf - expected_gf
        suspect = abs(delta_mw[i]) > mw_limit or abs(delta_gf[i]) > gf_limit
        if not suspect:
            phase1[i] = data.l1_cycles[i] - cumulative1
            phase2[i] = data.l2_cycles[i] - cumulative2
            good.append((t, current_mw, current_gf))
            statuses[i] = "ok"
            continue

        persists = _persistent_jump(
            i, data, mw, gf, mw_prediction=centre_mw,
            gf_prediction=expected_gf, mw_step=delta_mw[i],
            gf_step=delta_gf[i], mw_limit=mw_limit, gf_limit=gf_limit,
            cumulative1=cumulative1, cumulative2=cumulative2,
            max_gap_s=max_gap_s,
        )
        if persists is False:
            kind = ("code_outlier" if abs(delta_gf[i]) <= gf_limit
                    else "phase_outlier")
        elif persists is True:
            pair = _integer_pair(delta_mw[i], delta_gf[i], data,
                                 sigma_mw=max(.07, sigma_mw),
                                 sigma_gf=max(.003, sigma_gf))
            if pair is not None:
                step1, step2 = pair
                cumulative1 += step1
                cumulative2 += step2
                phase1[i] = data.l1_cycles[i] - cumulative1
                phase2[i] = data.l2_cycles[i] - cumulative2
                good.append((t, mw[i]-(cumulative1-cumulative2),
                             gf[i]-data.wavelength1_m*cumulative1
                             +data.wavelength2_m*cumulative2))
                kind = "repaired"
                events.append(CycleSlipEvent(i, float(t), kind, step1, step2,
                                             float(delta_mw[i]), float(delta_gf[i])))
                statuses[i] = kind
                continue
            kind = "unresolved"
        else:
            kind = "unresolved"

        events.append(CycleSlipEvent(i, float(t), kind,
                                     mw_step_cycles=float(delta_mw[i]),
                                     gf_step_m=float(delta_gf[i])))
        statuses[i] = kind
        if kind == "unresolved":
            # Never silently continue an uncertain carrier-phase ambiguity.
            arc += 1
            good.clear()
            cumulative1 = cumulative2 = 0

    return CycleSlipResult("TurboEdit-style", data, mw, gf, pred_mw, pred_gf,
                           delta_mw, delta_gf, phase1, phase2,
                           statuses, arcs, events)
