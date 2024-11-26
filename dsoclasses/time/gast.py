#! /usr/bin/python

from datetime import datetime
from scipy.spatial.transform import Rotation as R
import numpy as np

def approx_gast(t):
    sec_day = (t - t.replace(hour=0, minute=0, second=0, microsecond=0)).\
                total_seconds()
    omega_earth = 7.292115e-5 ## [rad/s]
    return omega_earth * float(sec_day)

def R3(theta):
    r = R.from_euler('z', theta, degrees=False).as_matrix()
