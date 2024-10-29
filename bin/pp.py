#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
import datetime
import argparse
import os, sys

from dsoclasses.orbits import sp3c, interpolator
from dsoclasses.geodesy import transformations
from dsoclasses.gnss import systems
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
parser.add_argument("--sat-id", 
    metavar='SAT_ID',
    dest='satid',
    required=False,
    default = None,
    help='Satellite id. Will be used to extract its position from the given sp3 file, hence, the id should match the id referenced therein. If not give, the first satellite encountered in the sp3 file will be used.')
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
    for arg in args: if args in dct: return dct[arg]
    return None

def main() -> int:

    try:
        args = parser.parse_args()

        if not os.path.isfile(args.sp3):
            print('Error. Failed to locate sp3 file {:}'.format(args.sp3), file=sys.stderr)
            return 1

# create an Sp3 instance,
        sp3 = sp3c.Sp3(args.sp3)
# set the id of the satellite we need,
        satid = args.satid if args.satid is not None else sp3.sat_ids[0]
# and extract its data
        data = sp3.get_satellite(satid, True)
        if data == {}:
            print('Error. Satellite {:} not found in Sp3 file.'.format(satid), file=sys.stderr)
            return 99
# construct an Interpolator
        intrp = interpolator.OrbitInterpolator(satid, data, 1800, 4, 'CubicSpline')

# Construct a RINEX instance to extract observations from
    rnx = GnssRinex(args.rnx)

# get approximate coordinates from RINEX header
    site_xyz = rnx.approx_cartesian()

# Loop throufh data block, and compute observations
    for block in rnx:
        try:
            t = block[t]
# consider only GPS satellite, loop through each GPS satellite
            for sat, obs in block.filter_satellite_system("gps", False)
                p1 = fetch(obs, 'C1P', 'C1L', 'C1C')['value']
                p2 = fetch(obs, 'C2L', 'C2W')['value']
                p3 = (GPS_L1_FREQ * GPS_L1_FREQ * p1 - GPS_L2_FREQ * GPS_L2_FREQ * p2) / (GPS_L1_FREQ*GPS_L1_FREQ - GPS_L2_FREQ * GPS_L2_FREQ)
# make L3 linear combination
        except:
            pass 
