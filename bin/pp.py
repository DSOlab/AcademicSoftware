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
parser.add_argument("--cut-off-angle", 
    metavar='CUTOFF_ANGLE',
    dest='cutoff',
    required=False,
    default = 10.,
    help='Cut-off angle to use in [deg].')
parser.add_argument("--exclude-sats", 
    metavar='EXCLUDE_SATS',
    dest='sat_exclude',
    required=False,
    default = [],
    help='Satellites to be excluded from the analysis; they should be passed in as recorded in the RINEX file, using a whitespace sperated list (e.g. \'--exclude-sats G13 R21 E09\'.')
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

def infoplot(xaxis, yaxis, dct, title='', satlist=[]):
    fig, ax = plt.subplots()
    available_sats = []
    for t, sats in dct.items():
        for sat in sats:
            if sat['sat'] not in available_sats:
                available_sats.append(sat['sat'])
    for plot_sat in available_sats:
        x=[]; y=[];
        if satlist == [] or (plot_sat in satlist):
            for t, tobs in dct.items():
                for sat in tobs:
                    if sat['sat'] == plot_sat:
                        if xaxis == 't': 
                            x.append(t)
                        else:
                            x.append(sat[xaxis])
                        y.append(sat[yaxis])
        ax.scatter(x, y, label=plot_sat)
    ax.legend()
    ax.grid(True)
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    plt.show()

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
        intrp = interpolator.Sp3Interpolator(args.sp3, [args.satsys], 1800, 4, args.interp_type, True, ['M', 'E'])

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
        dl = []; J = []; w = []; x0 = site_xyz + [0.]
        num_obs = 0; num_obs_used = 0;
        rawobs = {}
        sats_used = []
        residuals = {}

# Loop throufh data block, and compute observations
        for block in rnx:
            #try:
                t = block.t()
# consider only GPS satellite, loop through each GPS satellite
                for sat, obs in block.filter_satellite_system("gps", False):
                    if sat not in args.sat_exclude:
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
                            if el > np.radians(args.cutoff):
# make L3 linear combination
                                r, drdx, drdy, drdz = pseudorange(np.array(site_xyz), np.array((x,y,z)))
                                p3res = p3 + gs.C * clk - r - x0[3]
                                if abs(p3res) < 8e3:
                                    num_obs_used += 1
                                    dl.append(p3 + gs.C * clk - (r + x0[3]))
                                    J.append([drdx, drdy, drdz, 1e0])
                                    w.append(1. / np.cos(el))
                                    
                                    if t not in rawobs: rawobs[t] = []
                                    rawobs[t].append({'sat':sat, 'xyzsat': np.array((x,y,z)), 'p3': p3, 'el': np.degrees(el), 'clksat': clk, 'mark': 'used', 'res': dl[-1]})
# add some info
                                    if not sat in sats_used: sats_used.append(sat)
                                else:
                                    print("Residual too high! Observation skipped, sat={:}@{:}, res={:}km".format(sat, t, p3res*1e-3))
                        else:
                            # print("Failed to find observable for sat {:} at {:}".format(sat, t))
                            pass
                #except:
                #    pass

# Plots
        infoplot('t', 'res', rawobs, 'Raw Observation Residuals', [])
        infoplot('el', 'res', rawobs, 'Raw Observation Residuals', [])
        print("Used {:} out of {:} observations ~{:.1f}%".format(num_obs_used, num_obs, num_obs_used * 100. / num_obs))

# iterative Least Squares solution
        x0 = np.array(x0)
        for i in range(0,5):
            J  = np.array(J)
            dl = np.array(dl)
# var-covar matrix
            Q = np.linalg.inv(J.transpose() @ np.diag(w) @ J)
# estimate corrections
            dx = Q @ J.transpose() @ np.diag(w) @ dl
# compute new estimates
            x0 = x0 + dx
            lat, lon, hgt = transformations.car2ell(x0[0], x0[1], x0[2])
            Rt = transformations.geodetic2lvlh(lat, lon)
            R = Rt.transpose()
# compute residuals
            resi = []
            for t,v in rawobs.items():
                for satobs in v:
                    if satobs['mark'] == 'used':
                        xyzrec = np.array((x0[0], x0[1], x0[2]))
                        r, _, _, _ = pseudorange(xyzrec, satobs['xyzsat'])
                        res = satobs['p3'] + gs.C * satobs['clksat'] - (r + x0[3])
                        resi.append(res)
                        satobs['res'] = res
# variance of unit weight
            numobsi = dl.shape[0]
            numpars = 4
            u = np.array((resi))
            sigma = np.sqrt(u @ np.diag(w) @ u.transpose() / (numobsi - numpars))
            print("LS iteration {:}".format(i))
            print("Num of observations: {:}".format(numobsi))
            print("Sigma0 = {:}".format(np.sqrt(sigma)))

# remove observations with large residuals and prepare matrices for new iteration
            J = []; dl = []; w = [];
            obs_flagged = 0;
            for t,v in rawobs.items():
                for satobs in v:
                    if satobs['mark'] == 'used':
                        xyzrec = np.array((x0[0], x0[1], x0[2]))
                        r, drdx, drdy, drdz = pseudorange(xyzrec, satobs['xyzsat'])
                        az, el = azele(xyzrec, satobs['xyzsat'], R)
                        satobs['el'] = el
                        if (abs(satobs['res']) <= 3e0 * sigma) or (i<2):
                            dl.append(satobs['res'])
                            J.append([drdx, drdy, drdz, 1e0])
                            w.append(1. / np.cos(satobs['el']))
                        else:
                            satobs['mark'] = 'outlier'
                            obs_flagged += 1
            print("Observations marked: {:}".format(obs_flagged))
# plot residuals for this iteration
            print(dx)
            infoplot('t', 'res', rawobs,  'Observation Residuals iteration {:}'.format(i+1), [])
            infoplot('el', 'res', rawobs, 'Observation Residuals iteration {:}'.format(i+1), [])


    #except Exception as err:
    #    print("Error. Exception caught:", err)
    #    return 100

if __name__ == "__main__":
    sys.exit(main())
