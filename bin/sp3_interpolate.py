#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
import datetime
import argparse
import os, sys

from dsoclasses.orbits import sp3c, interpolator
from dsoclasses.geodesy import transformations
from dsoclasses.time import gast

## TODO: add clock interpolation (as cmd option)
##

parser = argparse.ArgumentParser()
parser.add_argument("--sp3", 
    metavar='SP3_FILE',
    dest='sp3',
    required=True,
    help='The sp3[cd] file to extract satellite position from.')
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
        intrp = interpolator.OrbitInterpolator(satid, data, 1800, 4, args.interp_type)

# interpolated values
        xi=[]; yi=[]; zi=[]; ti=[];
        t = sp3.t_start
        while t <= (sp3.t_start + datetime.timedelta(hours=24)):
            try:
                xc, yc, zc = intrp.at(t)
                ti.append(t)
                xi.append(xc); yi.append(yc); zi.append(zc);
            except:
                print("Warning. Failed interpolating at {:}".format(t), file=sys.stderr)
            t += datetime.timedelta(seconds=5)
        if len(ti) < 2: 
            print("Error. Too few orbital states to plot!", file=sys.stderr)
            return

# get the sp3 data for plotting
        sx=[]; sy=[]; sz=[]; st=[];
        for k,v in data.items():
            st.append(k); sx.append(v['x']); sy.append(v['y']); sz.append(v['z']);
        
# plot results and data points
        fig, axs = plt.subplots(3, 1, sharex=True)
        axs[0].scatter(st, sx, c='tab:red', s=.9, alpha=.5, edgecolors='none');#axs[0].set_title(r"Component $X$");
        axs[0].scatter(ti, xi, c='tab:green', s=.9, alpha=.5, edgecolors='none');#axs[0].set_title(r"Component $X$");
        axs[1].scatter(st, sy, c='tab:red', s=.9, alpha=.5, edgecolors='none');#axs[1].set_title(r"Component $Y$");
        axs[1].scatter(ti, yi, c='tab:green', s=.9, alpha=.5, edgecolors='none');#axs[1].set_title(r"Component $Y$");
        axs[2].scatter(st, sz, c='tab:red', s=.9, alpha=.5, edgecolors='none');#axs[2].set_title(r"Component $Z$");
        axs[2].scatter(ti, zi, c='tab:green', s=.9, alpha=.5, edgecolors='none');#axs[2].set_title(r"Component $Z$");
        fig.tight_layout()
        plt.show()

    except Exception as err:
        print("Error. Exception caught:", err)
        return 100

if __name__ == "__main__":
    sys.exit(main())