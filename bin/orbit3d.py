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
    t=[];x=[];y=[];z=[];
    for k,v in data.items():
# we are only plotting one day of orbit
        if (k-sp3.t_start).total_seconds() <= 86400e0:
            t.append(k)
            if ref_frame == 'ecef':
                x.append(v['x']); y.append(v['y']); z.append(v['z']); 
            else:
                R = Rotation.from_euler('z', gast.approx_gast(k), degrees=False)
                eci = np.array((v['x'], v['y'], v['z'])) @ R.as_matrix()
                x.append(eci[0]); y.append(eci[1]); z.append(eci[2]);
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z, 'gray')
    ax.scatter3D(x, y, z, c=z, cmap='Greens')
    ax.set_title('Orbit 3D Plot of Satellite {:} (frame: {:})'.format(satid, ref_frame.upper()))
    ax.set_xlabel('X', labelpad=10)
    ax.set_ylabel('Y', labelpad=10)
    ax.set_zlabel('Z', labelpad=10)
    plt.show()

if __name__ == "__main__":
    main()
