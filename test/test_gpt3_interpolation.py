#! /usr/bin/python

import sys
import datetime
from dsoclasses.troposphere import gpt3
import numpy as np
import random

lat = 35.78
lon = 23.79
print("input (lat, lon) = ({:.2f}, {:.2f})".format(lat, lon))
tl, tr, bl, br = gpt3.gpt3grid(np.radians(lon), np.radians(lat), sys.argv[1])
print(tl, tr, bl, br)

for i in range(1000):
    lon = random.uniform(-np.pi, np.pi)
    lat = random.uniform(-np.pi/2, np.pi/2)
    try:
        gpt3.gpt3grid(lon, lat, sys.argv[1])
        print("ok")
    except:
        print("Failed for (lat,lon)={:.2f}, {:.2f})".format(np.degrees(lat), np.degrees(lon)))
