"""Educational TEC-rate + Melbourne--Wuebbena cycle-slip detection.

Based on Z. Liu, published online in 2010, Journal of Geodesy 85 (2011),
171--183, doi:10.1007/s00190-010-0426-y. We predict the ionospheric TEC
rate (with a linear acceleration term) from preceding clean measurements and
jointly use the MW jump to solve both integer slips. This is a simplified,
offline teaching implementation, not a literal reproduction of every paper
parameter. The next-epoch persistence check is an offline teaching addition.

The published 0.15 TECU/s threshold was evaluated at ~1 Hz; this function
deliberately refuses sparse 30-second observations by default.
"""

from __future__ import annotations

import numpy as np

from .cycle_slip_common import (
    CycleSlipEvent, CycleSlipResult, DualFrequencySeries, _integer_pair,
    _local_stats, _persistent_jump, combinations,
)


def liu2010(data: DualFrequencySeries, *, tecr_threshold_tecu_s: float = .15,
            rate_history: int = 30, max_gap_s: float = 2.5) -> CycleSlipResult:
    """Locate and confidently repair dual-frequency integer cycle slips.

    Retains ~30 preceding, clean TEC rates to predict the next rate and
    acceleration. Both MW and TEC-rate evidence must agree with one integer
    pair; otherwise a break/outlier is explicitly recorded. Intended for 1 Hz
    series, including the ``teaching_gps_series`` demonstration.
    """
    if rate_history < 8 or tecr_threshold_tecu_s <= 0 or max_gap_s <= 0:
        raise ValueError("rate_history >= 8 and positive thresholds are required")
    if np.median(np.diff(data.t_s)) > 2.0:
        raise ValueError("Liu's 0.15 TECU/s test requires approximately 1 Hz "
                         "data; use the real 30-second RINEX only to explore "
                         "MW/GF, not to claim 1-Hz repair performance")

    mw, gf = combinations(data)
    n = len(mw)
    pred_mw = np.full(n, np.nan)
    pred_gf = np.full(n, np.nan)
    delta_mw = np.full(n, np.nan)
    delta_gf = np.full(n, np.nan)
    obs_tecr = np.full(n, np.nan)
    predict_tecr = np.full(n, np.nan)
    phase1 = np.full(n, np.nan)
    phase2 = np.full(n, np.nan)
    statuses = ["warmup"] * n
    arcs = np.zeros(n, dtype=int)
    events = []
    good = []  # clean (time, MW, GF), excluding isolated outliers
    cumulative1 = cumulative2 = 0
    arc = 0
    scale = data.gf_metres_per_tecu

    for i, t in enumerate(data.t_s):
        if i and t - data.t_s[i-1] > max_gap_s:
            events.append(CycleSlipEvent(i, float(t), "gap"))
            good.clear()
            cumulative1 = cumulative2 = 0
            arc += 1
        arcs[i] = arc
        current_mw = mw[i] - (cumulative1 - cumulative2)
        current_gf = gf[i] - data.wavelength1_m*cumulative1 + data.wavelength2_m*cumulative2

        if good:
            obs_tecr[i] = (current_gf-good[-1][2]) / (scale*(t-good[-1][0]))
        if len(good) < rate_history+1:
            phase1[i] = data.l1_cycles[i] - cumulative1
            phase2[i] = data.l2_cycles[i] - cumulative2
            good.append((t, current_mw, current_gf))
            continue

        recent = good[-rate_history:]
        centre_mw, sigma_mw = _local_stats([v[1] for v in recent],
                                           sigma_floor=.055)
        # TECU/s over the preceding 30 clean intervals, then estimate
        # current rate and its acceleration by a linear regression.
        rate_data = good[-(rate_history+1):]
        rate_t = np.asarray([(rate_data[j][0]+rate_data[j-1][0])/2
                             for j in range(1, len(rate_data))])
        rates = np.asarray([(rate_data[j][2]-rate_data[j-1][2]) /
                            (scale*(rate_data[j][0]-rate_data[j-1][0]))
                            for j in range(1, len(rate_data))])
        relative_t = rate_t - rate_t[-1]
        slope, intercept = np.polyfit(relative_t, rates, 1)
        estimated_rate = float(intercept+slope*((t+good[-1][0])/2-rate_t[-1]))
        # A smooth model's TEC-rate residual must exceed the threshold.
        fit_errors = rates - (slope*relative_t+intercept)
        _, sigma_rate = _local_stats(fit_errors, sigma_floor=.018)
        tecr_limit = max(tecr_threshold_tecu_s, 5*sigma_rate)
        dt_s = t-good[-1][0]
        expected_gf = good[-1][2] + estimated_rate*scale*dt_s
        mw_limit = max(.45, 5*sigma_mw)
        gf_limit = tecr_limit*scale*dt_s

        pred_mw[i], pred_gf[i] = centre_mw, expected_gf
        predict_tecr[i] = estimated_rate
        delta_mw[i] = current_mw-centre_mw
        delta_gf[i] = current_gf-expected_gf
        suspect = abs(delta_mw[i]) > mw_limit or abs(delta_gf[i]) > gf_limit
        if not suspect:
            phase1[i] = data.l1_cycles[i]-cumulative1
            phase2[i] = data.l2_cycles[i]-cumulative2
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
                                 sigma_gf=max(.003, sigma_rate*scale*dt_s))
            if pair is not None:
                step1, step2 = pair
                cumulative1 += step1
                cumulative2 += step2
                phase1[i] = data.l1_cycles[i]-cumulative1
                phase2[i] = data.l2_cycles[i]-cumulative2
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
            arc += 1
            good.clear()
            cumulative1 = cumulative2 = 0

    return CycleSlipResult("Liu TEC-rate + MW", data, mw, gf, pred_mw,
                           pred_gf, delta_mw, delta_gf, phase1, phase2,
                           statuses, arcs, events, obs_tecr, predict_tecr)
