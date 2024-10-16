#! /usr/bin/python

from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
import argparse
import os

from dsoclasses.orbits import sp3c
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
parser.add_argument("--ref-frame", 
    metavar='REF_FRAME',
    dest='ref_frame',
    required=False,
    default = 'ecef',
    choices=['ecef', 'eci'],
    help='Reference frame to plot the orbit at. Choose between \'ecef\' or \'eci\'')
parser.add_argument("--2d", 
    dest='plot2d',
    action='store_true',
    help='If set to true, the program will also produce 2D plots (i.e. time-vs-param) for position and velocity components.')

def main():

    args = parser.parse_args()
    ref_frame = args.ref_frame.lower()

    if not os.path.isfile(args.sp3):
        print('Error. Failed to locate sp3 file {:}'.format(args.sp3))
        sys.exit(1)

# create an Sp3 instance,
    sp3 = sp3c.Sp3(args.sp3)
# set the id of the satellite we need, 
    satid = args.satid if args.satid is not None else sp3.sat_ids[0]
# and extract its data
    data = sp3.get_satellite(satid, True)

# collect cartesian coordinates and epochs (t).
# if we are plotting eci, we need to transform the coordinates.
    t=[];x=[];y=[];z=[];vx=[];vy=[];vz=[];
    for k,v in data.items():
# we are only plotting one day of orbit
        if (k-sp3.t_start).total_seconds() <= 86400e0:
            t.append(k)
            if ref_frame == 'ecef':
                x.append(v['x']); y.append(v['y']); z.append(v['z']); 
                if args.plot2d:
                    vx.append(v['vx']); vy.append(v['vy']); vz.append(v['vz']); 
            else:
                R = Rotation.from_euler('z', gast.approx_gast(k), degrees=False)
                eci = np.array((v['x'], v['y'], v['z'])) @ R.as_matrix()
                x.append(eci[0]); y.append(eci[1]); z.append(eci[2]);
                if args.plot2d:
                    veci =  R.as_matrix() @ np.array((v['vx'], v['vy'], v['vz'])) + \
                            np.cross(np.array([0, 0, 7.292115e-5]), np.array((v['x'], v['y'], v['z'])))
                    vx.append(veci[0]); vy.append(veci[1]); vz.append(veci[2]);
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z, 'gray')
    ax.scatter3D(x, y, z, c=z, cmap='Greens')
    ax.set_title('Orbit 3D Plot of Satellite {:} (frame: {:})'.format(satid, ref_frame.upper()))
    ax.set_xlabel('X', labelpad=10)
    ax.set_ylabel('Y', labelpad=10)
    ax.set_zlabel('Z', labelpad=10)
    plt.show()
        
    if args.plot2d:
        fig, axs = plt.subplots(3,2)
        fig.suptitle(' [m]')
        axs[0,0].scatter(t,  x, alpha=.3, edgecolors='none', label=r"$X$")
        axs[1,0].scatter(t,  y, alpha=.3, edgecolors='none', label=r"$Y$")
        axs[2,0].scatter(t,  z, alpha=.3, edgecolors='none', label=r"$Z$")
        axs[0,1].scatter(t, vx, alpha=.3, edgecolors='none', label=r"$v_{x}$")
        axs[1,1].scatter(t, vy, alpha=.3, edgecolors='none', label=r"$v_{y}$")
        axs[2,1].scatter(t, vz, alpha=.3, edgecolors='none', label=r"$v_{z}$")
        plt.legend()
        plt.show()

if __name__ == "__main__":
    main()
