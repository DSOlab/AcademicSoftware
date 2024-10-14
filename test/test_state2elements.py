#! /usr/bin/python

import numpy as np
import math

from dsoclasses.orbits import elements

# example
r = np.array([-6045e3, -3490e3, 2500e3])
v = np.array([-3.457e3, 6.618e3, 2.533e3])
# in [km], [km^2] ...
refres = {'specific angular momentum': 58310,
            'inclination': math.radians(153.2e0),
            'Omega': math.radians(255.3e0),
            'eccentricity': 0.1712,
            'omega': math.radians(20.07e0),
            'true anomaly': math.radians(28.45e0),
            'perigee radius': 7284e0,
            'apogee radius':  10290e0,
            'semimajor': 8788e0,
            'period':  2.278e0*3600e0}

res = elements.state2elements(r,v,398600e9,True)
print("Equatorial cartesian to Elements:")
for k,d in res.items():
    print('{:20s} diff: {:+.2f}'.format(k, d-refres[k]))

re, ve = elements.elements2state(res, 398600e9, False)
print("Elements to cartesian:")
print("dx  = {:.0f} km".format(re[0]-r[0]))
print("dy  = {:.0f} km".format(re[1]-r[1]))
print("dz  = {:.0f} km".format(re[2]-r[2]))
print("dvx = {:.3f} km/sec".format(ve[0]-v[0]))
print("dvy = {:.3f} km/sec".format(ve[1]-v[1]))
print("dvz = {:.3f} km/sec".format(ve[2]-v[2]))
