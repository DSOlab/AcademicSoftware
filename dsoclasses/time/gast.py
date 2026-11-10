from datetime import datetime
from scipy.spatial.transform import Rotation as R
import numpy as np


def approx_gast(t):
    # Use the same calendar normalization as the RINEX and orbit modules;
    # native datetime, numpy datetime64 and attodatetime are all accepted.
    from dsoclasses.time.calmjd import cal2mjd
    _, sec_day = cal2mjd(t)
    omega_earth = 7.292115e-5  ## [rad/s]
    return omega_earth * float(sec_day)


def R3(theta):
    r = R.from_euler("z", theta, degrees=False).as_matrix()
