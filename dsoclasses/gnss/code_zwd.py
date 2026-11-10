"""Weighted code-only least squares for station XYZ, epoch clocks and ZWD.

This is a deliberately transparent batch solver for the ``tropo_zwd.ipynb``
academic example. Every iteration rebuilds the design matrix; the final
residuals are computed from the *updated* parameters, not from a previous
linearization. No interpolation of external products is done here.

Model (metres)::

    P3 + c * delta_t_sat = ||satellite_xyz - station_xyz||
                            + receiver_clock_range + ZHD*m_h + ZWD*m_w
                            + relativity_range + sagnac_range

The notebook prepares the code bias, antenna and atmospheric metadata. The
solver estimates a receiver-clock range for every kept epoch and a piecewise
linear ZWD through user-provided time knots.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from dsoclasses.gnss.systems import C


@dataclass
class CodeZwdSolution:
    """An auditable solution and its post-fit diagnostics."""
    parameters: np.ndarray
    residuals: np.ndarray  # observed minus *final* nonlinear prediction [m]
    weights: np.ndarray  # inverse variances [1/m^2], zero for rejected rows
    sigma_hat0: float
    iterations: int
    rank: int
    number_of_epochs: int
    number_of_zwd_knots: int


def solve_code_zwd(raw_obs, approx_xyz, n_epochs, *, sigma0=5.0,
                   max_iterations=8, max_rejections=3, rejection_sigma=3.0,
                   convergence_metres=1e-4):
    """Solve code-only XYZ, clocks and piecewise-linear zenith wet delays.

    A non-full-rank design, insufficient observations, non-finite input,
    nonconvergence, or a pathological geometry raises a specific error.
    Gross errors may be removed iteratively using the a priori per-row
    sigma ``sigma0 / sin(elevation)``. Rejected rows remain in the returned
    residual array but have zero statistical weight.
    """
    observations = list(raw_obs)
    if not observations:
        raise ValueError('No usable code observations; cannot solve for ZWD')
    if not isinstance(n_epochs, int) or n_epochs <= 0:
        raise ValueError('n_epochs must be positive')
    if not math.isfinite(float(sigma0)) or sigma0 <= 0:
        raise ValueError('sigma0 must be finite and positive [m]')
    if not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError('max_iterations must be positive')
    if not isinstance(max_rejections, int) or max_rejections < 0:
        raise ValueError('max_rejections must be non-negative')

    station0 = np.asarray(approx_xyz, dtype=float)
    if station0.shape != (3,) or not np.all(np.isfinite(station0)):
        raise ValueError('approx_xyz must be finite ECEF XYZ metres')

    # Index validation is performed before allocation; a missing knot cannot
    # otherwise be distinguished from a genuinely unobservable parameter.
    max_knot = -1
    for i, obs in enumerate(observations):
        if not 0 <= obs['nepoch'] < n_epochs:
            raise ValueError(f'Observation {i} references invalid clock epoch')
        if len(obs['zwd_basis']) != 2:
            raise ValueError(f'Observation {i} needs two ZWD basis weights')
        knots = obs['zwd_basis']
        if knots[1][0] != knots[0][0] + 1:
            raise ValueError(f'Observation {i} has non-adjacent ZWD knots')
        if knots[0][0] < 0:
            raise ValueError(f'Observation {i} has a negative ZWD knot index')
        if not math.isclose(knots[0][1] + knots[1][1], 1.0, abs_tol=1e-10):
            raise ValueError(f'Observation {i} has invalid ZWD basis normalization')
        max_knot = max(max_knot, knots[1][0])
    n_zwd = max_knot + 1
    n_params = 3 + n_epochs + n_zwd
    n_obs = len(observations)
    if n_obs <= n_params:
        raise ValueError(f'Not enough code observations: {n_obs} rows for {n_params} unknowns')
    x = np.concatenate([station0, np.zeros(n_epochs),
                        np.full(n_zwd, float(np.mean([obs['dT']['zwd']
                                                      for obs in observations])))])
    elevation = np.asarray([obs['el'] for obs in observations], dtype=float)
    if np.any(~np.isfinite(elevation)) or np.any(elevation <= 0) or np.any(elevation > math.pi / 2):
        raise ValueError('All elevations must be within (0, 90] degrees')
    sigma_i = sigma0 / np.sin(elevation)
    weights = 1.0 / sigma_i**2

    def linearize(parameters):
        """Build the nonlinear residual and model Jacobian at parameters."""
        residual = np.zeros(n_obs, dtype=float)
        jacobian = np.zeros((n_obs, n_params), dtype=float)
        for i, obs in enumerate(observations):
            sat = np.asarray(obs['rsat'], dtype=float)
            if sat.shape != (3,) or not np.all(np.isfinite(sat)):
                raise ValueError(f'Observation {i}: invalid satellite XYZ')
            diff = parameters[:3] - sat
            distance = float(np.linalg.norm(diff))
            if distance <= 0 or not math.isfinite(distance):
                raise ValueError(f'Observation {i}: invalid geometric range')
            dt_m = parameters[3 + obs['nepoch']]
            knot0, coeff0 = obs['zwd_basis'][0]
            knot1, coeff1 = obs['zwd_basis'][1]
            wet_delay = (coeff0 * parameters[3 + n_epochs + knot0]
                         + coeff1 * parameters[3 + n_epochs + knot1])
            dT = obs['dT']['zhd'] * obs['dT']['mfh'] + wet_delay * obs['dT']['mfw']
            computed = distance + dt_m + dT + obs['drel'] + obs['dsagnac']
            observed = obs['p3'] + C * obs['csat']
            residual[i] = observed - computed
            jacobian[i, :3] = diff / distance
            jacobian[i, 3 + obs['nepoch']] = 1
            jacobian[i, 3 + n_epochs + knot0] = obs['dT']['mfw'] * coeff0
            jacobian[i, 3 + n_epochs + knot1] = obs['dT']['mfw'] * coeff1
        if not np.all(np.isfinite(residual)) or not np.all(np.isfinite(jacobian)):
            raise ValueError('Model produced non-finite residuals or derivatives')
        return residual, jacobian

    iteration = 0
    rejected_rounds = 0
    for iteration in range(1, max_iterations + 1):
        residual, matrix = linearize(x)
        active = weights > 0
        if np.count_nonzero(active) <= n_params:
            raise ValueError('Too few observations remain after outlier rejection')
        sqrt_w = np.sqrt(weights[active])
        system = matrix[active] * sqrt_w[:, None]
        rhs = residual[active] * sqrt_w
        dx, _, rank, _ = np.linalg.lstsq(system, rhs, rcond=None)
        if rank != n_params:
            raise ValueError(f'Rank-deficient code/ZWD system: rank {rank}/{n_params}. '
                             'Check satellite geometry, clock epochs and ZWD knots.')
        x += dx
        residual, _ = linearize(x)  # Final-state residual; not stale dlw-Jw@dx.

        # Rejection is based on *a priori* per-observation standard deviations.
        # A very large noise scale should be corrected at the data-model level,
        # not silently used to inflate the empirical threshold.
        outliers = active & (np.abs(residual) / sigma_i > rejection_sigma)
        if np.any(outliers):
            if rejected_rounds >= max_rejections:
                raise RuntimeError('Too many rounds of gross-error rejection')
            weights[outliers] = 0.0
            rejected_rounds += 1
            continue  # Re-solve after removing gross errors.
        if np.linalg.norm(dx[:3]) < convergence_metres and np.max(np.abs(dx[3:])) < convergence_metres:
            break
    else:
        raise RuntimeError(f'Code/ZWD least squares did not converge in {max_iterations} iterations')

    inliers = weights > 0
    dof = int(np.count_nonzero(inliers) - n_params)
    if dof <= 0:
        raise ValueError('No positive degrees of freedom for variance estimation')
    sigma_hat = math.sqrt(float(np.dot(residual[inliers]**2, weights[inliers])) / dof)
    return CodeZwdSolution(parameters=x, residuals=residual, weights=weights,
                           sigma_hat0=sigma_hat, iterations=iteration, rank=n_params,
                           number_of_epochs=n_epochs, number_of_zwd_knots=n_zwd)
