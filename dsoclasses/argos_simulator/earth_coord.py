import numpy as np


def earth_coord(q1, q2, Ielem, Welem, Oelem):
    """
    Convert in-plane (perifocal) components to inertial components.

    This is a python translation of the EarthCoord.m function for argos simulator
    written by E. Schrama.

    Parameters
    ----------
    q1, q2 : float or np.ndarray
        In-plane components (along-perigee and perpendicular).
    Ielem : float or np.ndarray
        Inclination [rad].
    Welem : float or np.ndarray
        Argument of perigee [rad].
    Oelem : float or np.ndarray
        RAAN [rad].

    Returns
    -------
    coords : np.ndarray
        Inertial coordinates stacked as [..., 3] = (x, y, z).
    """
    q1 = np.asarray(q1, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    I = np.asarray(Ielem, dtype=float)
    W = np.asarray(Welem, dtype=float)
    O = np.asarray(Oelem, dtype=float)

    cI, sI = np.cos(I), np.sin(I)
    cW, sW = np.cos(W), np.sin(W)
    cO, sO = np.cos(O), np.sin(O)

    c1 = -cI * sW * sO + cW * cO
    c2 = -cI * cW * sO - sW * cO
    c3 = cI * sW * cO + cW * sO
    c4 = cI * cW * cO - sW * sO
    c5 = sI * sW
    c6 = sI * cW

    xs = q1 * c1 + q2 * c2
    ys = q1 * c3 + q2 * c4
    zs = q1 * c5 + q2 * c6

    return np.stack([xs, ys, zs], axis=-1)
