#! /usr/bin/python

import sys
import datetime
from dsoclasses.troposphere import gpt3
from dsoclasses.time import calmjd
import numpy as np
import random

lat = 35.78
lon = 23.79
print("input (lat, lon) = ({:.2f}, {:.2f})".format(lat, lon))
tl, tr, bl, br = gpt3.get_grid_nodes(np.radians(lon), np.radians(lat), sys.argv[1])
print(tl, tr, bl, br)

for i in range(1000):
    lon = random.uniform(-np.pi, np.pi)
    lat = random.uniform(-np.pi/2, np.pi/2)
    try:
        gpt3.get_grid_nodes(lon, lat, sys.argv[1])
        print("ok")
    except:
        print("Failed for (lat,lon)={:.2f}, {:.2f})".format(np.degrees(lat), np.degrees(lon)))


lat = np.pi/6.
lon = -np.pi/8.
hgt = 111.11
t = datetime.datetime(2024, 11, 4, 11)
print(gpt3.gpt3(t, lon, lat, hgt, sys.argv[1]))
print("gpt3_5({:.15e}, {:.15e}, {:.15e}, {:.10e}, {:d})".format(calmjd.cal2fmjd(t), lat, lon, hgt, 0))
