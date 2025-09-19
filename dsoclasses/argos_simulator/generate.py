import numpy as np


def generate(t, p, v, b, f0, minelev, c=3e8):
    """
    Generate visible-beacon observations and Doppler shifts.

    This is a python translation of the generate.m function for argos simulator
    written by E. Schrama.

    Parameters
    ----------
    t : array_like, shape (N,)
        Epochs [s].
    p : array_like, shape (N, 3)
        Satellite inertial/ECEF positions [m].
    v : array_like, shape (N, 3)
        Satellite inertial/ECEF velocities [m/s].
    b : array_like, shape (3,)
        Beacon position [m].
    f0 : float
        Beacon carrier frequency [Hz].
    minelev : float
        Minimum elevation angle at the beacon [rad].
    c : float, optional
        Speed of light [m/s]. Default 3e8 (to match original MATLAB).

    Returns
    -------
    idx : np.ndarray, dtype=int, shape (K,)
        Indices (into inputs) where the satellite is visible.
    time : np.ndarray, shape (K,)
        Observation times at those indices [s].
    doppler : np.ndarray, shape (K,)
        Observed Doppler shift at those times [Hz].
    """
    t = np.asarray(t, dtype=float)
    p = np.asarray(p, dtype=float)
    v = np.asarray(v, dtype=float)
    b = np.asarray(b, dtype=float).reshape(
        3,
    )

    # Basic shape checks
    if p.ndim != 2 or p.shape[1] != 3:
        raise ValueError("p must be shape (N, 3).")
    if v.ndim != 2 or v.shape[1] != 3 or v.shape[0] != p.shape[0]:
        raise ValueError("v must be shape (N, 3) and match p in length.")
    if t.shape[0] != p.shape[0]:
        raise ValueError("t length must match p/v length.")
    if b.shape != (3,):
        raise ValueError("b must be a 3-vector (shape (3,)).")

    rad = np.linalg.norm(b)  # beacon distance from origin
    if rad == 0.0:
        raise ValueError("Beacon position has zero norm; zenith angle is undefined.")

    idx_list = []
    time_list = []
    doppler_list = []

    for i in range(t.shape[0]):
        los = p[i] - b  # satellite -> beacon vector
        dis = np.linalg.norm(los)
        if dis == 0.0:
            # Degenerate (satellite at beacon): skip
            continue

        # cos(angle between los and beacon->origin)
        inpro = float(los @ b) / (dis * rad)
        inpro = np.clip(inpro, -1.0, 1.0)
        zenith = np.arccos(inpro)

        # Visible if elevation > minelev  <=>  zenith < (pi/2 - minelev)
        if zenith < (np.pi / 2.0 - minelev):
            # LOS unit vector from satellite toward beacon (negative sat->beacon)
            n_hat = -los / dis
            vee = float(n_hat @ v[i])  # radial velocity toward beacon

            df = (c / (c - vee)) * f0 - f0  # Doppler shift

            idx_list.append(i)
            time_list.append(t[i])
            doppler_list.append(df)

    return (
        np.asarray(idx_list, dtype=int),
        np.asarray(time_list, dtype=float),
        np.asarray(doppler_list, dtype=float),
    )
