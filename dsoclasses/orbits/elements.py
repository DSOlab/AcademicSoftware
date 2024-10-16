##
## Algorithm from:
## Orbital Mechanics for Engineering Students, 
## 

import sys
import numpy as np
from scipy.spatial.transform import Rotation as R

def kepler(mean_anomaly, eccentricity, epsilon=1e-15):
    """ Iteratively solve Kepler’s equation in for an object in an elliptical 
        orbit, i.e.
                              M + e * sinE = E
        Returns E, i.e. eccentric anomaly in [rad]
    """
    e = eccentricity
    M = mean_anomaly
# first guess
    Ei   = M
    Eip1 = -1e3
# perform iterations while |E(i) - E(i-1)| <= epsilon
# define MAX_ITERATIONS so that we avoid an infinte loop in case of 
# non-convergence
    MAX_ITERATIONS = 1000
    cit = 0
    while abs(Ei-Eip1) > epsilon and cit < MAX_ITERATIONS:
        Eip1 = M + e * np.sin(Ei)
        cit += 1

# check convergence
    if cit >= MAX_ITERATIONS:
        print('Error. Kepler\'s equation failed to converge! reached max_iterations, giving up!', file=sys.stderr)
        raise RuntimeError

    return Eip1

def true_anomaly(eccentric_anomaly, eccentricity):
    """ Compute true anomaly (v) given the eccentricity (e) and the
        eccentric anomaly (E), according to:
            v = E + 2 arctan( β sin E / (1 − β cos E) )
            with 
            β = e / (1 + √(1 − e²)).
    """
    e = eccentricity
    E = eccentric_anomaly
    beta = e / (1e0 + np.sqrt(1-e*e))
    v = E + 2 * np.arctan2(beta * np.sin(E), (1 - beta * np.cos(E)))
    return v

def eccentric_anomaly(true_anomaly, eccentricity):
    """ Compute eccentric anomaly (E) given the true anomaly (v).
        E = arctan( sinv*√(1 − e²) / (cosv + e) )
    """
    v = true_anomaly
    e = eccentricity
    return np.arctan2(np.sin(v) * np.sqrt(1e0-e*e), (np.cos(v) + e))

def mean_anomaly(eccentric_anomaly, eccentricity):
    """ Use Kepler's equation to compute mean anomaly (M) from eccentric 
        anomaly (E) and eccentricity (e)
        M + e * sinE = E
    """ 
    E = eccentric_anomaly
    e = eccentricity
    return E - np.sin(E) * e

class Coe:
    def __init__(self, dct):
        if 'mu' not in dct:
            self.mu=398600.435507
        else:
            self.mu = dct['mu']
        self.theta = dct['true anomaly']
        self.h     = dct['specific angular momentum']
        self.e     = dct['eccentricity']
        self.raan  = dct['Omega']
        self.w     = dct['omega']
        self.inc   = dct['inclination']

    def semimajor(self):
        return self.h*self.h / self.mu / (1. - self.e*self.e)

    def period(self):
        return 2.*np.pi / np.sqrt(self.mu) * self.semimajor()**(3./2.)

def elements2state(elements_dict, mu=398600.435507):
# get elements
    theta = elements_dict['true anomaly']
    h = elements_dict['specific angular momentum']
    e = elements_dict['eccentricity']
    raan = elements_dict['Omega']
    w = elements_dict['omega']
    inc = elements_dict['inclination']
# position and velocity in the orbital frame (z-axis perpendicular to
# orbital plane, x-axis pointing to periapsis of the orbit)
    rp = (h*h/mu) * (1./(1.+e*np.cos(theta))) * (np.cos(theta)*np.array([1.,0.,0.]) + np.sin(theta) * np.array([0.,1.,0.]))
    vp = (mu/h) * (-np.sin(theta)*np.array([1.,0.,0.]) + (e + np.cos(theta)) * np.array([0.,1.,0.]))
# transformation matrix R, from perifocal to geocentric equatorial coordinates
    Q = R.from_euler('ZXZ', [raan, inc, w]).as_matrix()
# transform to equatorial
    return Q @ rp, Q @ vp

def state2elements(pos, vel, mu=398600.435507):
    """ Computes the classical orbital elements (coe) from the state vector 
        (r,v).
        pos : position vector in the geocentric equatorial frame [km]
        vel : velocity vector in the geocentric equatorial frame [km/sec]
    """
    R = pos
    V = vel
    r = np.linalg.norm(R)
    v = np.linalg.norm(V)
# radial velocity
# Note that if vr > 0, the satellite is flying away from perigee. 
# If vr < 0, it is flying towards perigee.
    vr = np.dot(R,V)/r
# specific angular momentum and its magnitude [km^2/sec]
    H = np.cross(R,V)
    h = np.linalg.norm(H)
# Inclination (i)
# Recall that i must lie between 0 and 180 [deg], so there is no quadrant 
# ambiguity. If 90 < i ≤ 180, the orbit is retrograde.
    i = np.arccos(H[2]/h)
# vector defining the nodal line (and its magnitude)
    N = np.cross(np.array([0e0,0e0,1e0]),H)
    n = np.linalg.norm(N)
# Right Ascencion of the Ascending Node (Ω)
    Omega = np.arccos(N[0]/n)
    if N[1] < 0: Omega = 2*np.pi - Omega
# eccentricity vector and magnitude
    E = 1./mu * ((v*v - mu/r)*R - r*vr*V)
    e = np.linalg.norm(E)
# argument of perigee
    omega = np.arccos(np.dot(N,E)/n/e)
    if E[2] < 0: omega = 2*np.pi - omega
# true anomaly
    theta = np.arccos(np.dot(E,R)/e/r)
    if vr < 0: theta = 2*np.pi - theta 
# semi-major
    a = h*h/mu/(1.-e*e)

    return {'specific angular momentum': h,
            'inclination': i,
            'Omega': Omega,
            'eccentricity': e,
            'omega': omega,
            'true anomaly': theta}
