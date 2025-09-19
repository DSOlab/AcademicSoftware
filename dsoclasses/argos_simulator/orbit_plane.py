import numpy as np


def orbit_plane(a, e, mu, true_anom):
    """
    Vectorized translation of the MATLAB OrbitPlane function.

    This is a python translation of the OrbitPlane.m function for argos simulator
    written by E. Schrama.

    Parameters
    ----------
    a : float or array_like
        Semi-major axis [m].
    e : float or array_like
        Eccentricity.
    mu : float or array_like
        Gravitational parameter GM [m^3/s^2].
    true_anom : float or array_like
        True anomaly ν [rad].

    Returns
    -------
    q1, q2, r, v, q3, q4 : ndarray
    q1, q2: in-plane coordinates (perigee axis and perpendicular)
    r:      radius
    v:      speed
    q3, q4: in-plane velocity components (along/perp. to perigee axis)
    """
    a = np.asarray(a, dtype=float)
    e = np.asarray(e, dtype=float)
    nu = np.asarray(true_anom, dtype=float)
    mu = np.asarray(mu, dtype=float)

    # Broadcast all to a common shape
    a, e, mu, nu = np.broadcast_arrays(a, e, mu, nu)

    p = a * (1.0 - e * e)
    ct = np.cos(nu)
    st = np.sin(nu)

    r = p / (1.0 + e * ct)
    v = np.sqrt(mu * (2.0 / r - 1.0 / a))

    h0 = np.sqrt(mu / p)
    h1 = h0 * e * st
    h2 = h0 * (1.0 + e * ct)

    q1 = r * ct
    q2 = r * st
    q3 = h1 * ct - h2 * st
    q4 = h1 * st + h2 * ct

    return q1, q2, r, v, q3, q4
