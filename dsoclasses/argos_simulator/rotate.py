import numpy as np


def rotate(t, pos_eci, vel_eci, t0, a0):
    """
    Rotate inertial position/velocity into Earth-fixed using a z-rotation by GAST.
    This is a python translation of the rotate.m function for argos simulator
    written by E. Schrama.

    Rotate inertial (ECI) position/velocity to Earth-fixed (ECEF)
    using a z-rotation by angle a(t) = a0 + (t - t0)*dadt, with
    R3(-a) and dR3/dt = R3dot*(-dadt).

    Inputs:
      t        : (N,) times [s]
      pos_eci  : (N,3)
      vel_eci  : (N,3)

    Returns:
      pos_ecef, vel_ecef : (N,3)
    """
    t = np.asarray(t, dtype=float).ravel()
    px, py, pz = pos_eci.T
    vx, vy, vz = vel_eci.T

    dadt = 2.0 * np.pi / ((365.25 / 366.25) * 86400.0)  # rad/s
    a = a0 + (t - t0) * dadt
    ca = np.cos(-a)
    sa = np.sin(-a)

    # Rotate positions
    px1 = ca * px - sa * py
    py1 = sa * px + ca * py
    pz1 = pz

    # R3dot * p  (matrix from MATLAB code) * (-dadt)
    term = -dadt
    r3dot_p_x = ((-sa) * px + (-ca) * py) * term
    r3dot_p_y = ((ca) * px + (-sa) * py) * term
    r3dot_p_z = 0.0

    # R3 * v
    rvx = ca * vx - sa * vy
    rvy = sa * vx + ca * vy
    rvz = vz

    vx1 = r3dot_p_x + rvx
    vy1 = r3dot_p_y + rvy
    vz1 = r3dot_p_z + rvz  # r3dot_p_z = 0

    pos_ecef = np.column_stack((px1, py1, pz1))
    vel_ecef = np.column_stack((vx1, vy1, vz1))
    return pos_ecef, vel_ecef
