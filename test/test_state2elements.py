#! /usr/bin/python

import numpy as np
import math

from dsoclasses.orbits import elements

# example
r = np.array([-6045., -3490., 2500.])   # [km]
v = np.array([-3.457, 6.618, 2.533]) # [km/s]

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
print("Angular momentum (kmˆ2/s) = {:.1f}".format(res['specific angular momentum']))
print("Eccentricity              = {:.6f}".format(res['eccentricity']))
print("Right ascension (deg)     = {:.3f}".format(np.degrees(res['Omega'])))
print("Inclination (deg)         = {:.3f}".format(np.degrees(res['inclination'])))
print("Argument of perigee (deg) = {:.4f}".format(np.degrees(res['omega'])))
print("True anomaly (deg)        = {:.4f}".format(np.degrees(res['true anomaly'])))
print("Semimajor axis (km)       = {:.1f}".format(elements.Coe(res).semimajor()))
r2, v2 = elements.elements2state(res)
print("r (km)   = ({:.1f} {:.1f} {:.1f})".format(r2[0], r2[1], r2[2]))
print("v (km/s) = ({:.4f} {:.4f} {:.4f})".format(v2[0], v2[1], v2[2]))

r2, v2 = elements.elements2state({'specific angular momentum': 80000, 
    'eccentricity': 1.4, 
    'Omega': np.radians(40.), 
    'inclination': np.radians(30.), 
    'omega': np.radians(60.), 
    'true anomaly': np.radians(30.)})
print("r (km)   = ({:.1f} {:.1f} {:.1f})".format(r2[0], r2[1], r2[2]))
print("v (km/s) = ({:.4f} {:.4f} {:.4f})".format(v2[0], v2[1], v2[2]))
res = elements.state2elements(r2, v2)
print("Angular momentum (kmˆ2/s) = {:.1f}".format(res['specific angular momentum']))
print("Eccentricity              = {:.6f}".format(res['eccentricity']))
print("Right ascension (deg)     = {:.3f}".format(np.degrees(res['Omega'])))
print("Inclination (deg)         = {:.3f}".format(np.degrees(res['inclination'])))
print("Argument of perigee (deg) = {:.4f}".format(np.degrees(res['omega'])))
print("True anomaly (deg)        = {:.4f}".format(np.degrees(res['true anomaly'])))
print("Semimajor axis (km)       = {:.1f}".format(elements.Coe(res).semimajor()))

r_vec = np.array((2777.102781, 1606.637620, -7017.697233))
v_vec = np.array((-21929.090388, 65692.475315, 6358.213206))
res = elements.state2elements(r_vec, v_vec)
print("Equatorial cartesian to Elements:")
print("Angular momentum (kmˆ2/s) = {:.1f}".format(res['specific angular momentum']))
print("Eccentricity              = {:.6f}".format(res['eccentricity']))
print("Right ascension (deg)     = {:.3f}".format(np.degrees(res['Omega'])))
print("Inclination (deg)         = {:.3f}".format(np.degrees(res['inclination'])))
print("Argument of perigee (deg) = {:.4f}".format(np.degrees(res['omega'])))
print("True anomaly (deg)        = {:.4f}".format(np.degrees(res['true anomaly'])))
print("Semimajor axis (km)       = {:.1f}".format(elements.Coe(res).semimajor()))
r2, v2 = elements.elements2state(res)
print("r (km)   = ({:.1f} {:.1f} {:.1f})".format(r2[0], r2[1], r2[2]))
print("v (km/s) = ({:.4f} {:.4f} {:.4f})".format(v2[0], v2[1], v2[2]))

