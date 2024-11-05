#! /usr/bin/python

import sys
import datetime
from dsoclasses.troposphere import gmf
import numpy as np
import random

lat = 35.78
lon = 23.79
t = datetime.datetime(2024, 11, 4, 11)

f1, f2 = gmf.gmf(t, np.radians(lat), np.radians(lon), 111.11, np.radians(53.78))
print("{:.6f} {:.6f}".format(f1,f2))
