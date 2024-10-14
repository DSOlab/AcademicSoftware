#! /usr/bin/python

from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
import argparse
import os

from dsoclasses.orbits import sp3c, elements
from dsoclasses.geodesy import transformations
from dsoclasses.time import gast

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

    t=[];
    a=[]; e=[]; i=[]; raan=[]; omega=[]; u=[]; h=[]; T=[];
    for k,v in data.items():
# we are only plotting one day of orbit
        if (k-sp3.t_start).total_seconds() <= 86400e0:
            t.append(k)
# ECEF to ECI (approximate)
            R = Rotation.from_euler('z', gast.approx_gast(k), degrees=False)
            recef = np.array((v['x'], v['y'], v['z']))
            vecef = np.array((v['vx'], v['vy'], v['vz']))
            reci =  R.as_matrix() @ recef
            veci =  R.as_matrix() @ vecef +\
                    np.cross(np.array([0, 0, 7.292115e-5]), recef)
# cartesian (equatorial) to elements
            d = elements.state2elements(reci, veci)
# store elements to plot
            a.append(d['semimajor'] * 1e-6)
            e.append(d['eccentricity'])
            i.append(np.degrees(d['inclination']))
            raan.append(np.degrees(d['Omega']))
            omega.append(np.degrees(d['omega']))
            u.append(np.degrees(d['true anomaly']))
            h.append(d['specific angular momentum'] * 1e-6) ## to km^2
            T.append(d['period'] / 3600e0 ) ## to hours

    fig, axs = plt.subplots(4, 2)
    axs[0, 0].plot(t, i);    axs[0, 0].set_title(r"Inclination $i$ [$\degree$]");
    axs[1, 0].plot(t, raan); axs[1, 0].set_title(r"RAAN $\Omega$ [$\degree$]");
    axs[2, 0].plot(t, omega);axs[2, 0].set_title(r"Argument of Periapsis $\omega$ [$\degree$]"); 
    axs[3, 0].plot(t, u);    axs[3, 0].set_title(r"True Anomaly $\theta$ [$\degree$]");
    axs[0, 1].plot(t, a);    axs[0, 1].set_title(r"Semimajor Axis $a$ [km]");
    axs[1, 1].plot(t, e);    axs[1, 1].set_title(r"Eccentricity $e$ [-]");
    axs[2, 1].plot(t, h);    axs[2, 1].set_title(r"Specific Angular Momentum $h$ [$km^2/s$]");
    axs[3, 1].plot(t, T);    axs[3, 1].set_title(r"Period $T$ [hours]");
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
