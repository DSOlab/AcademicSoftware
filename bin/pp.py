#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
import datetime, attotime
import argparse
import os, sys

from dsoclasses.orbits import sp3c, interpolator
from dsoclasses.geodesy import transformations
from dsoclasses.gnss import systems as gs
from dsoclasses.rinex.gnss.rinex import GnssRinex
from dsoclasses.time import gast
from dsoclasses.troposphere import gmf, gpt3

def athash(self): 
    s = "{:}{:}{:}{:}{:}{:}{:}{:}".format(self.year, self.month, self.day, self.hour, self.minute, self.second, self.microsecond, self.nanosecond)
    return int(s)
attotime.attodatetime.__hash__ = athash

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
    type=float,
    default = 10.,
    help='Cut-off angle to use in [deg].')
parser.add_argument("--sigma0", 
    metavar='SIGMA0_APRIORI',
    dest='sigma0',
    required=False,
    type=float,
    default = 10.,
    help='A-priori std. deviation of unit weight [m].')
parser.add_argument("--rcvr-clock-poly", 
    metavar='DELTA_T_RECEIVER_ORDER',
    dest='rcvr_clk_order',
    required=False,
    type=int,
    default = 3,
    help='Order of polynomial to model the receiver clock offset.')
parser.add_argument("--gpt3-grid", 
    metavar='GPT35_GRID_FILE',
    dest='gpt3',
    required=False,
    default = None,
    help='The GPT3 (5x5 degree) data file to compute tropospheric delay. Can be downloaded from TU Vienna at https://vmf.geo.tuwien.ac.at/codes/gmf.m')
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

def infoplot(xaxis, yaxis, dct, title='', satlist=[], use_marked=False):
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
                    if sat['sat'] == plot_sat and (sat['mark'] == 'used' or use_marked == True):
                        if xaxis == 't': 
                            x.append(t)
                        elif xaxis == 'el':
                            x.append(np.degrees(sat[xaxis]))
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

def tropo(t, lat, lon, hgt, gpt3_grid):
    if gpt3_grid is None: return 0., 0.
    meteo = gpt3.gpt3(t, lon, lat, hgt, gpt3_grid)
    zhd = gpt3.saastamoinen_zhd(lat, hgt, meteo['p'])
    zwd = gpt3.askne_zwd(meteo['e'], meteo['Tm'], meteo['la'])
    return zhd, zwd

def dtrp(t, lat, lon, hgt, el, zhd, zwd):
    gmfh, gmfw = gmf.gmf(t, lat, lon, hgt, np.pi/2-el)
    return zhd * gmfh + zwd * gmfw

def weight(el):
    """ Higher elevation means higher weight !
    """
    assert el >= 0. and el <= np.pi / 2
    return np.sin(el) * np.sin(el)
    # return 1. / np.cos(el)

def main() -> int:

    #try:
        args = parser.parse_args()

        if not os.path.isfile(args.sp3):
            print('Error. Failed to locate sp3 file {:}'.format(args.sp3), file=sys.stderr)
            return 1

# construct an Interpolator
        intrp = interpolator.Sp3Interpolator(args.sp3, [args.satsys], 3600, 10, args.interp_type, True, ['M', 'E'])

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

# compute ZHD and ZWD (if needed)
        zhd, zwd = tropo(rnx.time_first_obs, lat, lon, hgt, args.gpt3)

# LS matrices and info
        dl = []; J = []; w = []; x0 = site_xyz + [0.]
        num_obs = 0; num_obs_used = 0;
        sigma0 = args.sigma0
        var0 = sigma0 * sigma0
        rawobs = {}
        sats_used = []
        num_epochs = 0

# partials w.r.t receiver clock polynomial model: 
# a0 + a1*dt + a2*dt^2 + ae*dt^3 + ...
# returned as [1, dt, dt^2 ...]
        clk_partials = lambda dt : [ np.power(dt,n) for n in range(args.rcvr_clk_order+1) ]
        clk_value_at = lambda x,dt : sum([a*np.power(dt,n) for n,a in enumerate(x[3:])])

        f1 = gs.GPS_L1_FREQ
        f2 = gs.GPS_L2_FREQ
        f31 = f1*f1/(f1*f1 - f2*f2)
        f32 = f2*f2/(f1*f1 - f2*f2)

# Loop through the RINEX file, block-by-block and collect infor for further
# processing. While loopoing, formulate the LS matrices to be solved for once 
# we are through with the file.
        for block in rnx:
                t = block.t() ## current epoch (of obs/block)
                epoch_used = False
# consider only GPS satellite; loop through each GPS satellite in the block
                for sat, obs in block.filter_satellite_system("gps", False):
# if the satellite is not in eclusion list,
                    if sat not in args.sat_exclude:
                        num_obs += 1 ## total number of observations
# get code observations for L1 and L2
                        p1 = fetch(obs, 'C1W', 'C1C', 'C1X')
                        p2 = fetch(obs, 'C2W', 'C2C', 'C2D', 'C2P', 'C2X', 'C2S', 'C2L')
# get satellite position in ITRF (interpolation)
                        try:
                            x,y,z,clk = intrp.sat_at(sat, t)
                        except:
                            x = None
# if we do have P1 and P2 AND we have satellite coordinates and clock 
# correction, carry on ...
                        if (p1 is not None and p2 is not None) and (x is not None and clk is not None):
                            p1 = p1['value']
                            p2 = p2['value']
# compute azimouth and elevation (receiver-to-satellite)
                            az, el = azele(np.array(site_xyz), np.array((x,y,z)), R)
# if we are above cut-off angle ...
                            if el > np.radians(args.cutoff):
# make P3 linear combination
                                #p3 = (gs.GPS_L1_FREQ * gs.GPS_L1_FREQ * p1 - gs.GPS_L2_FREQ * gs.GPS_L2_FREQ * p2) / (gs.GPS_L1_FREQ*gs.GPS_L1_FREQ - gs.GPS_L2_FREQ * gs.GPS_L2_FREQ)
                                p3 = f31 * p1 - f32 * p2
# compute receiver-to-satellite distance and partials w.r.t. r_rec
                                r, drdx, drdy, drdz = pseudorange(np.array(site_xyz), np.array((x,y,z)))
# tropospheric correction
                                dT = dtrp(t, lat, lon, hgt, el, zhd, zwd)
# time from t0 (i.e. start of RINEX) for receiver clock correction
                                dsec = float((t-rnx.time_first_obs).total_nanoseconds()) * 1e-9
# residual: observed - computed
                                p3res = p3 + gs.C * clk - (r + x0[3] + dT)
                                if abs(p3res) < 8e3:
                                    num_obs_used += 1 # count used observations
                                    epoch_used = True # epoch is used
                                    dl.append(p3res)  # append residual to vector
# Jacobian
                                    J.append([drdx, drdy, drdz]+clk_partials(dsec))
                                    w.append(weight(el)) # weight of obs
# store information of this observation for later use
                                    if t not in rawobs: rawobs[t] = []
                                    rawobs[t].append({'sat':sat, 'xyzsat': np.array((x,y,z)), 'p3': p3, 'el': el, 'clksat': clk, 'mark': 'used', 'res': dl[-1], 'dsec': dsec})
# add satellite to used satellites
                                    if not sat in sats_used: sats_used.append(sat)
                                else:
                                    print("Residual too high! Observation skipped, sat={:}@{:}, res={:}km".format(sat, t, p3res*1e-3))
                        else:
                            # print("Failed to find observable for sat {:} at {:}".format(sat, t))
                            pass
# augment num_epochs if any observation for this epoch was used
                num_epochs += int(epoch_used)

# Plots
        # infoplot('t', 'res', rawobs, 'Raw Observation Residuals', [])
        # infoplot('el', 'res', rawobs, 'Raw Observation Residuals', [])
        print("Used {:} out of {:} observations ~{:.1f}%".format(num_obs_used, num_obs, num_obs_used * 100. / num_obs))

# iterative Least Squares solution:
# ---------------------------------------------------------------------------
        P = lambda var, weights: np.diag([xj * xj / 10. for xj in weights])
        x0 = np.concatenate((np.array(x0[0:3]), np.zeros(args.rcvr_clk_order+1)))
        dx = np.ones(x0.shape[0])
# convergence criteria
        MAX_ITERATIONS = 10
        iteration = 0
        max_diffs = [1e-2, 1e-2, 1e-2] + [1e-9 for i in range(args.rcvr_clk_order+1)]
        # print(max_diffs)
        # print(np.absolute(dx))
        while np.any(np.greater_equal(np.absolute(dx), max_diffs)) and iteration < MAX_ITERATIONS:
            iteration += 1
            J  = np.array(J)  # to numpy matrix
            dl = np.array(dl) # to numpy vector
# var-covar matrix (un-scaled)
            N = np.linalg.inv(J.transpose() @ P(var0, w) @ J)
# estimate corrections
            dx = (N @ J.transpose() @ P(var0, w) @ dl)
# compute new estimates x_new <- x_apriori + dx
            x0 = x0 + dx
# get the cartesian to topocentric rotation matrix (again)
            lat, lon, hgt = transformations.car2ell(x0[0], x0[1], x0[2])
            Rt = transformations.geodetic2lvlh(lat, lon)
            R = Rt.transpose()
# compute residuals for all observations that are marked 'used'
            resi = []
            for t, v in rawobs.items():
                for satobs in v:
                    if satobs['mark'] == 'used':
                        r, _, _, _ = pseudorange(np.array((x0[0], x0[1], x0[2])), satobs['xyzsat'])
                        dT = dtrp(t, lat, lon, hgt, satobs['el'], zhd, zwd)
                        res = satobs['p3'] + gs.C * satobs['clksat'] - (r + dT + clk_value_at(x0, satobs['dsec']))
                        resi.append(res)
                        satobs['res'] = res
                        # print("r={:5.1f}km el={:5.1f}[deg] p={:.3f}".format(res*1e-3, np.degrees(satobs['el']), weight(satobs['el'])*weight(satobs['el']))) 
# variance of unit weight
            numobsi = dl.shape[0]
            numpars = 3 + args.rcvr_clk_order + 1
            u = np.array((resi))
            sigma0 = np.sqrt(u @ P(var0, w) @ u.transpose() / (numobsi - numpars))
            var0 = sigma0 * sigma0
            print("LS iteration {:}".format(iteration))
            print("-----------------------------------------------")
            print("Num of observations: {:}".format(numobsi))
            print("Sigma0 = {:.1f}[m]".format(sigma0))

# variance of residuals
            # Vres = var0 * P(1./var0, [1./wi for wi in w]) - var0 * (J @ N @ J.transpose())

# remove observations with large residuals and prepare matrices for new iteration
            J = []; dl = []; w = [];
            obs_flagged = 0;
            for t, v in rawobs.items():
                for satobs in v:
                    if satobs['mark'] == 'used':
                        xyzrec = np.array((x0[0], x0[1], x0[2]))
                        r, drdx, drdy, drdz = pseudorange(xyzrec, satobs['xyzsat'])
                        az, el = azele(xyzrec, satobs['xyzsat'], R)
                        satobs['el'] = el
                        if (abs(satobs['res']) <= 10e0 * sigma0):
                            dl.append(satobs['res'])
                            J.append([drdx, drdy, drdz]+clk_partials(satobs['dsec']))
                            w.append(weight(satobs['el'])) # weight of obs
                        else:
                            satobs['mark'] = 'outlier'
                            obs_flagged += 1
            print("Num. observations marked: {:}".format(obs_flagged))
# print estimates
            print("Estimates:") 
            print("{:15.3f} +- {:7.3f} ({:+6.3f}) [m]".format(x0[0], np.sqrt(var0*N[0,0]), dx[0]))
            print("{:15.3f} +- {:7.3f} ({:+6.3f}) [m]".format(x0[1], np.sqrt(var0*N[1,1]), dx[1]))
            print("{:15.3f} +- {:7.3f} ({:+6.3f}) [m]".format(x0[2], np.sqrt(var0*N[2,2]), dx[2]))
            for j in range(args.rcvr_clk_order+1): print("{:.3f} +- {:.3f}[nsec]".format((x0[3+j]/gs.C)*1e9, np.sqrt(var0*N[3+j,3+j])), end=' ')
            print("")

# plot residuals for this iteration
            if iteration >= 2:
                # infoplot('t', 'res', rawobs,  'Observation Residuals iteration {:}'.format(i+1), [])
                infoplot('el', 'res', rawobs, 'Observation Residuals iteration {:}'.format(iteration+1), [])

    #except Exception as err:
    #    print("Error. Exception caught:", err)
    #    return 100

if __name__ == "__main__":
    sys.exit(main())
