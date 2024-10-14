#! /usr/bin/python

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import argparse
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

# concatenate satellite position vectors, and transform geocentric cartesian 
# crd to ellipsoidal
    lons = []
    lats = []
    for k,v in data.items():
        # we are only plotting one day of orbit
        if (k-sp3.t_start).total_seconds() <= 86400e0:
            phi, lamda, _ = transformations.car2ell(v['x'], v['y'], v['z'])
            lons.append(lamda)
            lats.append(phi)

# llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
# are the lat/lon values of the lower left and upper right corners
# of the map.
# lat_ts is the latitude of true scale.
# resolution = 'c' means use crude resolution coastlines.
    m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
                llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m.drawcoastlines()
    m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='aqua')

# convert to map projection coords.
    xpt,ypt = m(np.degrees(lons),np.degrees(lats))
    m.plot(xpt,ypt,'go',linewidth=.1,markersize=1)

    plt.title("Groundtrack of Satellite {:} (Mercator Projection)".format(satid))
    plt.show()

if __name__ == "__main__":
    main()
