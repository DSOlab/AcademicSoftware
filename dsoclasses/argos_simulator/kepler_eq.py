import numpy as np


def kepler_eq(t, t0, mu, a, e, nloop=5):
    """
    Solve Kepler's equation by fixed-point iteration.

    This is a python translation of the KeplerEq.m function for argos simulator
    written by E. Schrama.

    Parameters
    ----------
    t : float or array_like
        Epoch(s) [s].
    t0 : float or array_like
        Time of last perigee passage [s].
    mu : float or array_like
        Gravitational parameter μ [m^3/s^2].
    a : float or array_like
        Semi-major axis [m].
    e : float or array_like
        Eccentricity (0 <= e < 1 for elliptic orbits).
    nloop : int, optional
        Number of fixed-point iterations for E = M + e sin(E).
        (Default: 5, matching the MATLAB style loop.)

    Returns
    -------
    true_anom : np.ndarray
        True anomaly ν(t) [rad].
    ecc_anom : np.ndarray
        Eccentric anomaly E(t) [rad].
    mean_anom : np.ndarray
        Mean anomaly M(t) [rad].
    check : np.ndarray
        Residual of Kepler’s equation: E - e sin(E) - M (should be ~0).

    Notes
    -----
    - This is a direct translation of the MATLAB routine that uses
      simple fixed-point iteration (not Newton). Increase `nloop`
      if you need tighter convergence for high eccentricities.
    - All inputs are broadcast to a common shape; outputs have that shape.
    """
    t = np.asarray(t, dtype=float)
    t0 = np.asarray(t0, dtype=float)
    mu = np.asarray(mu, dtype=float)
    a = np.asarray(a, dtype=float)
    e = np.asarray(e, dtype=float)

    # Broadcast everything to a common shape
    t, t0, mu, a, e = np.broadcast_arrays(t, t0, mu, a, e)

    n = np.sqrt(mu / a**3)  # mean motion [rad/s]
    mean_anom = n * (t - t0)  # M(t)
    ecc_anom = mean_anom.copy()  # initial guess E0 = M

    for _ in range(int(nloop)):
        ecc_anom = mean_anom + e * np.sin(ecc_anom)

    check = ecc_anom - e * np.sin(ecc_anom) - mean_anom
    true_anom = np.arctan2(np.sqrt(1.0 - e**2) * np.sin(ecc_anom), np.cos(ecc_anom) - e)

    return true_anom, ecc_anom, mean_anom, check
