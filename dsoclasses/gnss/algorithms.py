"""GNSS geometry and signal propagation utilities (SI units throughout)."""
from __future__ import annotations

import math

import numpy as np

from dsoclasses.gnss import systems as gs
from dsoclasses.time.pyattotime import shift_epoch_seconds


def _xyz(value, what: str) -> np.ndarray:
    """Validate a finite 3-vector before using it in geometric computations."""
    xyz = np.asarray(value, dtype=float)
    if xyz.shape != (3,) or not np.all(np.isfinite(xyz)):
        raise ValueError(f'{what} must be a finite ECEF XYZ 3-vector [m]')
    return xyz


def geometric_range(rsta, rsat) -> float:
    """Euclidean range [m] between station and satellite ECEF positions."""
    return float(np.linalg.norm(_xyz(rsat, 'Satellite') - _xyz(rsta, 'Station')))


def _subtract_seconds(epoch, seconds):
    """Subtract short durations through the canonical high-precision API.

    The interpolation query now accepts attodatetime, datetime and numpy
    datetime64 uniformly. With attotime installed, light-time shifts retain
    sub-microsecond precision even for a native datetime reception epoch.
    """
    return shift_epoch_seconds(epoch, -seconds)


def sat_at_emission_time(rsta, t_reception, sp3_interpolator, sat,
                         tolerance=(0.1, 0.1, 0.1), max_iterations=10):
    """Find satellite position and clock at the time of signal transmission.

    The observed reception epoch is iteratively moved backwards by the
    vacuum geometric light-travel time ``rho / c``. The algorithm *does not*
    apply receiver/satellite clock, Sagnac, or atmospheric corrections; those
    are separate terms in the observation equation.

    Parameters
    ----------
    rsta : sequence of three floats
        Approximate ECEF receiver position in metres.
    t_reception : datetime-like
        Receive epoch, on the same timescale as the SP3 interpolator.
    sp3_interpolator : object
        Exposes ``sat_at(satellite_id, epoch) -> (x, y, z, clock)``.
    sat : str
        Satellite identifier, e.g. G21.
    tolerance : tuple of three positive floats
        Maximum per-coordinate change [m] between consecutive iterations.
    max_iterations : int
        Hard limit preventing accidental infinite loops.

    Returns
    -------
    (x, y, z, clock, flight_time)
        Position in metres, SP3 clock in seconds, flight time in seconds.
    """
    station = _xyz(rsta, 'Station')
    tol = np.asarray(tolerance, dtype=float)
    if tol.shape != (3,) or np.any(~np.isfinite(tol)) or np.any(tol <= 0):
        raise ValueError('tolerance must contain three positive finite metres')
    if not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError('max_iterations must be a positive integer')

    position = _xyz(sp3_interpolator.sat_at(sat, t_reception)[:3], 'Satellite')
    for _ in range(max_iterations):
        flight_time = float(np.linalg.norm(position - station) / gs.C)
        if not math.isfinite(flight_time) or flight_time <= 0:
            raise ValueError('Satellite/receiver geometry gives invalid flight time')
        emission_epoch = _subtract_seconds(t_reception, flight_time)
        x, y, z, clock = sp3_interpolator.sat_at(sat, emission_epoch)
        new_position = _xyz((x, y, z), 'Satellite')
        if not math.isfinite(float(clock)):
            raise ValueError(f'{sat} has no finite satellite clock at emission')
        if np.all(np.abs(new_position - position) < tol):
            # The last iteration's coordinates may shift the flight time by
            # less than tol/c. Return a time consistent with the final XYZ.
            final_flight_time = float(np.linalg.norm(new_position - station) / gs.C)
            return float(x), float(y), float(z), float(clock), final_flight_time
        position = new_position
    raise RuntimeError(f'Emission-time iteration did not converge for {sat} '
                       f'after {max_iterations} iterations')


def satellite_velocity(epoch, interpolator, sat, step_seconds=0.5) -> tuple[np.ndarray, np.ndarray]:
    """Central-difference satellite ECEF velocity and position [m/s, m].

    A symmetric difference has lower truncation error than the notebook's
    original forward difference. Queries must have interpolator support on
    *both* sides of the target epoch; no extrapolation is performed.
    """
    if not math.isfinite(float(step_seconds)) or step_seconds <= 0:
        raise ValueError('step_seconds must be finite and positive')
    before = _subtract_seconds(epoch, step_seconds)
    # Calling with a negative interval adds seconds, for both datetime
    # and attotime-like epochs.
    after = _subtract_seconds(epoch, -step_seconds)
    left = _xyz(interpolator.sat_at(sat, before)[:3], 'Satellite')
    right = _xyz(interpolator.sat_at(sat, after)[:3], 'Satellite')
    position = _xyz(interpolator.sat_at(sat, epoch)[:3], 'Satellite')
    return (right - left) / (2.0 * step_seconds), position


def relativistic_range_term(epoch, interpolator, sat, step_seconds=0.5) -> float:
    """Relativistic eccentricity term in metres, 2*(r dot v)/c.

    **Sign convention:** This is the term moved to the *computed* side of
    ``P + c*SP3_clock = geometric_range + relativity + ...`` in the teaching
    notebook. The satellite-clock correction itself is ``-2*r.v/c**2``.

    WARNING: Before enabling this term with a particular precise clock
    product, verify whether its clock convention *already* includes the
    relativistic effect. Double-applying the correction biases pseudoranges.
    """
    velocity, position = satellite_velocity(epoch, interpolator, sat, step_seconds)
    return float(2.0 * np.dot(position, velocity) / gs.C)


def sagnac_range_term(rsat, rsta) -> float:
    """First-order Earth-rotation correction [m] in receive-time geometry.

    Sign convention: add to the *computed* range with ECEF r_sat at signal
    emission time and r_sta at reception time.
    """
    satellite = _xyz(rsat, 'Satellite')
    station = _xyz(rsta, 'Station')
    return float((gs.OmegaEarth / gs.C) *
                 (satellite[0] * station[1] - satellite[1] * station[0]))
