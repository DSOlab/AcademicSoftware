#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
import datetime
import argparse
import os, sys

from dsoclasses.orbits import sp3c, interpolator
from dsoclasses.geodesy import transformations
from dsoclasses.gnss import systems as gs
from dsoclasses.rinex.gnss.rinex import GnssRinex
from dsoclasses.time import gast

parser = argparse.ArgumentParser()
parser.add_argument("--sp3", 
    metavar='SP3_FILE',
    dest='sp3',
    required=True,
    help='The sp3[cd] file to extract satellite position from.')
parser.add_argument("--data-rinex", 
    metavar='DATA_RINEX',
    dest='rnx',
    required=True,
    help='The RINEX file to extract observation data from. Currently, RINEX v3.x are supported.')
parser.add_argument("--sat-sys", 
    metavar='SAT_SYS',
    dest='satsys',
    required=False,
    default = 'G',
    help='Satellite system(s) to be used for the Point Processing.')
parser.add_argument("--interpolation", 
    metavar='INTERP_ALGORITHM',
    dest='interp_type',
    required=False,
    default = 'CubicSpline',
    choices=['Polynomial', 'CubicSpline', 'PchipInterpolator'],
    help='Choose interpolation algorithm to be used.')

def fetch(dct, *args):
    """ Given a dictionary containing e.g.
        R05 : {'C1C': {'value': 23539032.631, 'lli': None, 'ssi': 6}, 'L1C': {'value': 125829717.51, 'lli': 0, 'ssi': 6}, 'D1C': {'value': -4149.772, 'lli': None, 'ssi': 6}, 'S1C': {'value': 41.719, 'lli': None, 'ssi': None}, 'C1P': {'value': 23539032.74, 'lli': None, 'ssi': 6}, 'L1P': {'value': 125829714.502, 'lli': 0, 'ssi': 6}, 'D1P': {'value': -4149.698, 'lli': None, 'ssi': 6}, 'S1P': {'value': 41.062, 'lli': None, 'ssi': None}, 'C2P': {'value': 23539038.067, 'lli': None, 'ssi': 6}, 'L2P': {'value': 97867622.91, 'lli': 0, 'ssi': 6}, 'D2P': {'value': -3227.451, 'lli': None, 'ssi': 6}, 'S2P': {'value': 38.531, 'lli': None, 'ssi': None}, 'C2C': {'value': 23539037.837, 'lli': None, 'ssi': 6}, 'L2C': {'value': 97867623.908, 'lli': 0, 'ssi': 6}, 'D2C': {'value': -3227.359, 'lli': None, 'ssi': 6}, 'S2C': {'value': 38.531, 'lli': None, 'ssi': 6}}
        return the observation (dictionary) first encountered, matched by *args.
        E.g. if the above dictionary is stored in dct, 
        fetch(dct, 'C1P', 'C1C', 'C2P') 
        will return {'value': 23539032.74, 'lli': None, 'ssi': 6}
    """
    for arg in args:
        if arg in dct:
            return dct[arg]
    return None

def pseudorange(rrec, rsat):
    """ Computes:
        r:    distance between satellite and receiver,
        drdx: partial of distance w.r.t. x (receiver) 
        drdy: partial of distance w.r.t. y (receiver) 
        drdz: partial of distance w.r.t. z (receiver)
    """
    r = np.linalg.norm(rsat-rrec)
    drdx = -(rsat[0] - rrec[0]) / r
    drdy = -(rsat[1] - rrec[1]) / r
    drdz = -(rsat[2] - rrec[2]) / r
    return r, drdx, drdy, drdz

def azele(rrec, rsat, R):
    """ compute azimouth and elevation 
    """
    r = np.linalg.norm(rsat-rrec)
    enu = R @ (rsat-rrec)
    az = np.arctan2(enu[0], enu[1])
    el = np.arcsin(enu[2] / r)
    return az, el

def main() -> int:

    #try:
        args = parser.parse_args()

        if not os.path.isfile(args.sp3):
            print('Error. Failed to locate sp3 file {:}'.format(args.sp3), file=sys.stderr)
            return 1

# construct an Interpolator
        intrp = interpolator.Sp3Interpolator(args.sp3, [args.satsys], 1800, 4, args.interp_type)

# Construct a RINEX instance to extract observations from
        rnx = GnssRinex(args.rnx)
        if not rnx.time_sys == intrp.time_sys:
            print("Error. RINEX time system ({:}) is not the same as the SP3 {:}".format(rnx.time_sys, intrp.time_sys), file=sys.stderr)
            assert rnx.time_sys == intrp.time_sys

# get approximate coordinates from RINEX header
        site_xyz = rnx.approx_cartesian()
        lat, lon, hgt = transformations.car2ell(site_xyz[0], site_xyz[1], site_xyz[2])
        Rt = transformations.geodetic2lvlh(lat, lon)
        R = Rt.transpose()

# LS matrices and info
        l = []; J = []; w = []; x0 = site_xyz + [0.]
        num_obs = 0; num_obs_used = 0;
        sats_used = []

# Loop throufh data block, and compute observations
        for block in rnx:
            #try:
                t = block.t()
# consider only GPS satellite, loop through each GPS satellite
                for sat, obs in block.filter_satellite_system("gps", False):
                    num_obs += 1
                    p1 = fetch(obs, 'C1P', 'C1L', 'C1C')
                    p2 = fetch(obs, 'C2L', 'C2W')
# get satellite position in ITRF
                    try:
                        x,y,z,clk = intrp.sat_at(sat, t)
                    except:
                        x = None
                    if p1 is not None and p2 is not None and x is not None:
                        p1 = p1['value']
                        p2 = p2['value']
                        p3 = (gs.GPS_L1_FREQ * gs.GPS_L1_FREQ * p1 - gs.GPS_L2_FREQ * gs.GPS_L2_FREQ * p2) / (gs.GPS_L1_FREQ*gs.GPS_L1_FREQ - gs.GPS_L2_FREQ * gs.GPS_L2_FREQ)
# compute azimouth and elevation
                        az, el = azele(np.array(site_xyz), np.array((x,y,z)), R)
                        if el > np.radians(10):
                            num_obs_used += 1
# make L3 linear combination
                            r, drdx, drdy, drdz = pseudorange(np.array(site_xyz), np.array((x,y,z)))
                            l.append(p3 - gs.C * clk - r + gs.C * x0[3])
                            J.append([drdx, drdy, drdz, 1e0])
                            w.append(1. / np.cos(el))
# add some info
                            if not sat in sats_used: sats_used.append(sat)
                    else:
                        print("Failed to find observable for sat {:} at {:}".format(sat, t))
            #except:
            #    pass

# Least Squares solution
        J = np.array(J); l = np.array(l);
        print(J.shape)
        print(np.diag(w).shape)
        print(l.shape)
        dx = np.linalg.inv(J.transpose() @ np.diag(w) @ J) @ J.transpose() @ np.diag(w) @ l
        print(dx)

# List of satellites used in the analysis
        print("Satellites used: #{:} ->".format(len(sats_used)), end='')
        for s in sats_used: print("{:} ".format(s), end='')
        print()

    #except Exception as err:
    #    print("Error. Exception caught:", err)
    #    return 100

if __name__ == "__main__":
    sys.exit(main())
