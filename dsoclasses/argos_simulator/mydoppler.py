import numpy as np


def mydoppler(idx, time, Dobs, t_all, pos_ecef, vel_ecef, beacon_xyz, f0):
    """
    Compute modeled Doppler at the same indices (non-vectorized loop to match MATLAB).

    This is a direct, loop-based translation of the original MATLAB `mydoppler`
    and keeps the same call signature (even though `time`, `Dobs`, and `t`
    are not used in the math). Indices are 0-based (Python-style).

    This is a python translation of the mydoppler.m function for argos simulator
    written by E. Schrama.
    """
    c = 3e8
    xb0, yb0, zb0 = beacon_xyz
    Dcom = np.zeros_like(time, dtype=float)

    for i in range(len(idx)):
        j = idx[i]
        xlos = pos_ecef[j, 0] - xb0
        ylos = pos_ecef[j, 1] - yb0
        zlos = pos_ecef[j, 2] - zb0
        dis = np.sqrt(xlos**2 + ylos**2 + zlos**2)

        xn = -xlos / dis
        yn = -ylos / dis
        zn = -zlos / dis

        vee = xn * vel_ecef[j, 0] + yn * vel_ecef[j, 1] + zn * vel_ecef[j, 2]
        Dcom[i] = (c / (c - vee)) * f0 - f0

    return Dcom
