import numpy as np
from dsoclasses.argos_simulator.kepler_eq import kepler_eq
from dsoclasses.argos_simulator.orbit_plane import orbit_plane
from dsoclasses.argos_simulator.earth_coord import earth_coord


# assumes you already defined:
#   kepler_eq(t, t0, mu, a, e, nloop=100)
#   orbit_plane(a, e, mu, true_anom)
#   earth_coord(q1, q2, Ielem, Welem, Oelem)
# (the Python versions you translated earlier)


def kepler_orbit(major, eccen, incli, ranode, omega, t, t0, gm, nloop=100):
    """
    Keplerian elements -> inertial position & velocity.

    Parameters
    ----------
    major : float or array_like
        Semi-major axis a [m].
    eccen : float or array_like
        Eccentricity e (0 <= e < 1).
    incli : float or array_like
        Inclination i [rad].
    ranode : float or array_like
        RAAN Ω [rad].
    omega : float or array_like
        Argument of perigee ω [rad].
    t : float or array_like
        Epoch(s) [s].
    t0 : float or array_like
        Time of last perigee passage [s].
    gm : float or array_like
        Gravitational parameter μ [m^3/s^2].
    nloop : int, optional
        Iterations for Kepler solver (fixed-point). Default 100.

    Returns
    -------
    pos : np.ndarray, shape (..., 3)
        Inertial position [m].
    vel : np.ndarray, shape (..., 3)
        Inertial velocity [m/s].
    """
    # 1) True anomaly
    true_anom, _, _, _ = kepler_eq(t, t0, gm, major, eccen, nloop=nloop)

    # 2) Perifocal (in-plane) components
    q1, q2, _, _, q3, q4 = orbit_plane(major, eccen, gm, true_anom)

    # 3) Rotate perifocal -> inertial
    pos = earth_coord(q1, q2, incli, omega, ranode)  # (..., 3)
    vel = earth_coord(q3, q4, incli, omega, ranode)  # (..., 3)

    return pos, vel
