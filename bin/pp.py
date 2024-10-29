#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
import datetime
import argparse
import os, sys

from dsoclasses.orbits import sp3c, interpolator
from dsoclasses.geodesy import transformations
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
# consider only GPS satellite
            gps = block.filter_satellite_system("gps", False)
# for every satellite observation block
            for obs in gps:
                fetch(obs, 'C1L', 'C1C', 'C1L')
        except:
            pass 
