#! /usr/bin/python

import numpy as np
import math

from dsoclasses.orbits import elements

# example
r = np.array([-6045e3, -3490e3, 2500e3])   # [m]
v = np.array([-3.457e3, 6.618e3, 2.533e3]) # [m/s]
# in [m], [m^2] ...
refres = {'specific angular momentum': 58310e6,
            'inclination': math.radians(153.2e0),
            'Omega': math.radians(255.3e0),
            'eccentricity': 0.1712,
            'omega': math.radians(20.07e0),
            'true anomaly': math.radians(28.45e0),
            'perigee radius': 7284e3,
            'apogee radius':  10290e3,
            'semimajor': 8788e3,
            'period':  2.278e0*3600e0}

res = elements.state2elements(r, v)
print("Equatorial cartesian to Elements:")
for k,d in res.items():
    print('{:20s} diff: {:+.2f}'.format(k, d-refres[k]))

re, ve = elements.elements2state(res)
print("Elements to cartesian:")
print("dx  = {:.0f} ({:.0f}-{:.0f}) m".format(    re[0]-r[0], re[0], r[0]))
print("dy  = {:.0f} ({:.0f}-{:.0f}) m".format(    re[1]-r[1], re[1], r[1]))
print("dz  = {:.0f} ({:.0f}-{:.0f}) m".format(    re[2]-r[2], re[2], r[2]))
print("dvx = {:.3f} ({:.3f}-{:.3f}) m/sec".format(ve[0]-v[0], ve[0], v[0]))
print("dvy = {:.3f} ({:.3f}-{:.3f}) m/sec".format(ve[1]-v[1], ve[1], v[1]))
print("dvz = {:.3f} ({:.3f}-{:.3f}) m/sec".format(ve[2]-v[2], ve[2], v[2]))
