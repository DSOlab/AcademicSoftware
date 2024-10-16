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
parser.add_argument("--check-transformation", 
    action='store_true',
    dest='check_transformation',
    help='If set to true, then the program will check and plot the state-to-elements transformation, i.e. it will transform back, from the computed elements to ECI coordinates and produce a plot with respective differences. The results should be nearly zero (both for position and velocity).')

def main() -> int:

    try:
        args = parser.parse_args()

        if not os.path.isfile(args.sp3):
            print('Error. Failed to locate sp3 file {:}'.format(args.sp3), file=sys.stderr)
            return 1

# create an Sp3 instance,
        sp3 = sp3c.Sp3(args.sp3)
# make sure the sp3 has both position and velocity
        if sp3.pos_vel != 'V':
            print('Sp3 file {:} is marked as \'P\', i.e. does not contain velocity estimates. Giving up!'.format(sp3.fn))
            return 2

# set the id of the satellite we need, 
        satid = args.satid if args.satid is not None else sp3.sat_ids[0]
# and extract its data
        data = sp3.get_satellite(satid, False)

        t=[];
        a=[]; e=[]; i=[]; raan=[]; omega=[]; u=[]; h=[]; T=[];
        dr_coc = []; dv_coc=[]; dr_ioi = []; dv_ioi=[];
        for k,v in data.items():
# we are only plotting one day of orbit
            if (k-sp3.t_start).total_seconds() <= 86400e0:
                t.append(k)
# ECEF to ECI (approximate)
                R = Rotation.from_euler('z', gast.approx_gast(k), degrees=False)
                recef = np.array((v['x'], v['y'], v['z']))
                vecef = np.array((v['vx']*1e-4, v['vy']*1e-4, v['vz']*1e-4)) # dm/sec to km/sec
                reci =  R.as_matrix() @ recef
                veci =  R.as_matrix() @ vecef +\
                        np.cross(np.array([0, 0, 7.292115e-5]), recef)
                reci=recef
                veci=vecef
# cartesian (equatorial) to elements
                d = elements.state2elements(reci, veci)

# store elements to plot
                a.append(elements.Coe(d).semimajor())
                e.append(d['eccentricity'])
                i.append(np.degrees(d['inclination']))
                raan.append(np.degrees(d['Omega']))
                omega.append(np.degrees(d['omega']))
                u.append(np.degrees(d['true anomaly']))
                h.append(d['specific angular momentum'])
                T.append(elements.Coe(d).period() / 3600e0 ) ## to hours

# statistics of conversion
                if args.check_transformation:
                    r2, v2 = elements.elements2state(d)
                    dr_coc.append(1e3*(r2 - reci))
                    dv_coc.append(1e3*(v2 - veci))
                    recef2 = R.inv().as_matrix() @ r2
                    vecef2 = R.inv().as_matrix() @ v2 - np.cross(np.array([0, 0, 7.292115e-5]), r2)
                    dr_ioi.append(recef2-recef)
                    dv_ioi.append(vecef2-vecef)

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

        if args.check_transformation:
            fig, axs = plt.subplots(2)
            fig.suptitle('ECI -> Elements -> ECI [m]')
            axs[0].scatter(t, [x[0] for x in dr_coc], alpha=.3, edgecolors='none', label=r"$\delta X$")
            axs[0].scatter(t, [x[1] for x in dr_coc], alpha=.3, edgecolors='none', label=r"$\delta Y$")
            axs[0].scatter(t, [x[2] for x in dr_coc], alpha=.3, edgecolors='none', label=r"$\delta Z$")
            axs[1].scatter(t, [x[0] for x in dv_coc], alpha=.3, edgecolors='none', label=r"$\delta v_{x}$")
            axs[1].scatter(t, [x[1] for x in dv_coc], alpha=.3, edgecolors='none', label=r"$\delta v_{y}$")
            axs[1].scatter(t, [x[2] for x in dv_coc], alpha=.3, edgecolors='none', label=r"$\delta v_{z}$")
            plt.legend()
            plt.show()

    except Exception as err:
        print("Error. Exception caught:", err)
        return 100

if __name__ == "__main__":
    sys.exit(main())
