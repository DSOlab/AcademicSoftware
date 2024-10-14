#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import argparse
import math
import os

from dsoclasses.orbits import sp3c
from dsoclasses.geodesy import transformations

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
parser.add_argument("--base-sta", 
    dest='site',
    required=False,
    default = (4595220.016, 2039434.081, 3912626.007),
    nargs=3,
    type=float,
    metavar=('X', 'Y', 'Z'),
    help='')

def main():

    args = parser.parse_args()

    if not os.path.isfile(args.sp3):
        print('Error. Failed to locate sp3 file {:}'.format(args.sp3))
        sys.exit(1)

# create an Sp3 instance,
    sp3 = sp3c.Sp3(args.sp3)
# set the id of the satellite we need, 
    satid = args.satid if args.satid is not None else sp3.sat_ids[0]
# and extract its data
    data = sp3.get_satellite(satid, True)

# get the ellipsoidal coordinates of the site
    lat, lon, hgt = transformations.car2ell(args.site[0], args.site[1], args.site[2])

# compute the (e,n,u) to (x, y, z) rotation matrix
    Rt = transformations.geodetic2lvlh(lat, lon)
# anf get it's inverse
    R = Rt.transpose()

# iterate through the satellite positions and compute azimouth and zenith 
# angle (i.e, 90-elevation)
# Also collect satellite "entry points", i.e. points on the orbital arc where 
# the satellite starts to be above the (local) horizon
    azs = []; els = [];
    entries_az = []; entries_el = []; entries_t = [];
    last_elevation_was = -math.pi
    for k,v in data.items():
# we are only plotting one day of orbit
        if (k-sp3.t_start).total_seconds() <= 86400e0:
# cartesina difference vector
            dx = v['x'] - args.site[0]
            dy = v['y'] - args.site[1]
            dz = v['z'] - args.site[2]
# topocentric vector
            enu = R @ [dx, dy, dz]
# compute azimouth and zenith angle
            r = np.linalg.norm(enu)
            az = math.atan2(enu[0], enu[1])
            el = math.asin(enu[2] / r) 
            if el >= 0:
                azs.append(az)
                els.append(90e0 - math.degrees(el))
# store if entry point
            if last_elevation_was <= 0 and el > 0:
                entries_az.append(az)
                entries_el.append(90e0 - math.degrees(el))
                entries_t.append(k)
            last_elevation_was = el

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(azs, els, 'bo', markersize=1)
    ax.plot(entries_az, entries_el, 'rx')
    for i in range(len(entries_t)):
        ax.text(entries_az[i], entries_el[i], entries_t[i].strftime("%y:%m:%dT%H:%M"),
        fontsize=7, ha='right', va='bottom')
    ax.set_rmax(90)
    ax.set_rticks([30, 60])  # Less radial ticks
    ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.grid(True)

# Rotate the angle axis by -90 degrees (make 0 degrees point up)
    ax.set_theta_zero_location('N')  # Set 0° to the top (North)
    ax.set_theta_direction(-1)  # Make the angles increase clockwise

    ax.set_title("Skyplot of Satellite {:}".format(satid), va='bottom')
    plt.show()

if __name__ == "__main__":
    main()
