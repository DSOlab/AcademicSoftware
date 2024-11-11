import numpy as np
import datetime
import attotime

# mean gravity in m/s**2
gm = 9.80665e0
# molar mass of dry air in kg/mol
dMtr = 28.965e-3
# universal gas constant in J/K/mol
Rg = 8.3143e0

def saastamoinen_zhd(lat, h, p):
    """ Input parameters:
        lat: Ellipsoidal latitide in [rad]
        h  : Height above the ellipsoid [m]
        p  : Total surface pressure in [hPa]
        Returns:
        ZHD in [m]
    """
    h *= 1e-3 
    return 22768e-7 * p / (1. - 266e-5 * np.cos(2.*lat) - 28e-5*h)

def askne_zwd(e, Tm, wvlr):
    """ Input parameters:
        e:    water vapor pressure in hPa 
        Tm:   mean temperature in Kelvin
        wvlr: water vapor lapse rate
        Returns:
        ZWD in [m]
    """
    k1  = 77.604
    k2 = 64.79
    k2p = k2 - k1*18.0152/28.9644
    k3  = 377600.

    # specific gas constant for dry consituents
    Rd = Rg/dMtr
    return 1e-6*(k2p + k3/Tm)*Rd/(wvlr + 1.)/gm*e

def parse_gpt3_line(line):
    data_keys = ['lat', 'lon', 'p:', 'T:', 'Q:', 'dT:', 'undu', 'Hs', 'a_h:', 'a_w:', 'lambda:', 'Tm:', 'Gn_h:', 'Ge_h:', 'Gn_w:', 'Ge_w:']
    key_coeffs = ['a0', 'A1', 'B1', 'A2', 'B2']
    data = [float(x) for x in line.split()]
    idx = 0
    dct = {}
    for key in data_keys:
        if key.endswith(':'):
            scale = 1.
            if key[0:-1] in ['Q', 'dT', 'a_h', 'a_w']:
                scale = 1e-3
            elif key[0:-1] in ['Gn_h', 'Ge_h', 'Gn_w', 'Ge_w']:
                scale = 1e-6
            dct[key[0:-1]] = {}
            for kf in key_coeffs:
                dct[key[0:-1]][kf] = data[idx] * scale
                idx += 1
        else:
            dct[key] = data[idx]
            idx += 1
    return dct

def bilinear_interpolation(doy, lat, lon, hgt, tld, trd, bld, brd):
    cfy = np.cos(doy/365.25*2*np.pi)   # coefficient for A1
    chy = np.cos(doy/365.25*4*np.pi)   # coefficient for B1
    sfy = np.sin(doy/365.25*2*np.pi)   # coefficient for A2
    shy = np.sin(doy/365.25*4*np.pi)   # coefficient for B2
# transforming ellipsoidal height to orthometric height
    hgt = [ hgt - d['undu'] for d in [tld, trd, bld, brd] ]
    fun = lambda dct, arg: dct[arg]['a0'] + dct[arg]['A1'] * cfy + dct[arg]['B1'] * sfy + dct[arg]['A2'] * chy + dct[arg]['B2'] * shy
# temperature at the height of the grid
    T0 = [ fun(d,'T') for d in [tld, trd, bld, brd] ]
# pressure at the height of the grid
    p0 = [ fun(d,'p') for d in [tld, trd, bld, brd] ]
# humidity
    Ql = [ fun(d,'Q') for d in [tld, trd, bld, brd] ]
# reduction = stationheight - gridheight
    redh = [ h - d['Hs'] for (h,d) in zip(hgt,[tld, trd, bld, brd]) ]
# lapse rate of the temperature in degree / m
    dTl = [ fun(d,'dT') for d in [tld, trd, bld, brd] ]
# temperature reduction to station height
    Tl = [ t0 + dtl*redh_ - 273.15 for (t0,dtl,redh_) in zip(T0,dTl,redh) ]
# virtual temperature
    Tv = [ t0*(1+0.6077 * ql) for (t0,ql) in zip(T0,Ql) ]
    c =  [ gm*dMtr/(Rg*tv) for tv in Tv ]
# pressure in hPa
    pl = [ (p_*np.exp(-c_*redh_))/100. for (p_, c_, redh_) in zip(p0,c,redh) ]
# hydrostatic coefficient ah
    ahl = [ fun(d,'a_h') for d in [tld, trd, bld, brd] ]
# wet coefficient aw
    awl = [ fun(d,'a_w') for d in [tld, trd, bld, brd] ]
# water vapor decrease factor la
    lal = [ fun(d,'lambda') for d in [tld, trd, bld, brd] ]
# mean temperature of the water vapor Tm
    Tml = [ fun(d,'Tm') for d in [tld, trd, bld, brd] ]
# water vapor pressure in hPa 
    e0 = [ql*p_/(0.622+0.378*ql)/100. for (ql,p_) in zip(Ql,p0) ]
    el = [e_*(100.*pl_/p_)**(lal_+1.) for (e_,pl_,p_,lal_) in zip(e0,pl,p0,lal)]
# gradients
    Gn_hl = [ fun(d,'Gn_h') for d in [tld, trd, bld, brd] ]
    Ge_hl = [ fun(d,'Ge_h') for d in [tld, trd, bld, brd] ]
    Gn_wl = [ fun(d,'Gn_w') for d in [tld, trd, bld, brd] ]
    Ge_wl = [ fun(d,'Ge_w') for d in [tld, trd, bld, brd] ]
# interpolate
    y2 = tld['lat']
    y1 = bld['lat']
    y  = np.degrees(lat)
    x2 = trd['lon']
    x1 = tld['lon']
    x  = np.degrees(lon)
    assert x2 >= x and x >= x1
    assert y2 >= y and y >= y1
    bi = lambda Q12,Q22,Q11,Q21: (1./((x2-x1)*(y2-y1))) * (Q11*(x2-x)*(y2-y) + Q21*(x-x1)*(y2-y) + Q12*(x2-x)*(y-y1) + Q22*(x-x1)*(y-y1))
    return {'p': bi(*pl), 'T': bi(*Tl), 'dT': bi(*dTl)*1e3, 'e': bi(*el), 'ah': bi(*ahl), 'aw': bi(*awl), 'la': bi(*lal), 'Tm': bi(*Tml), 'Gn_h': bi(*Gn_hl), 'Ge_h': bi(*Ge_hl), 'Gn_w': bi(*Gn_wl), 'Ge_w': bi(*Ge_wl)}


def get_grid_nodes(lon, lat, grid):
    """ Grid file: https://vmf.geo.tuwien.ac.at/codes/gpt3_5.grd
    """

# we want positive longitude for the grid
    dlon = np.degrees(lon) + 360. if np.degrees(lon)<0. else np.degrees(lon)
    if dlon < 2.5 or dlon > 360.-2.5:
        print("Error. Longitude out of gpt3 grid file range.")
        raise RuntimeError
    dlat = np.degrees(lat)
    if dlat < -87.5 or dlat > 87.5:
        print("Error. Latitude out of gpt3 grid file range.")
        raise RuntimeError

    lons_per_lat = int(np.floor(((360.0-2.5) - 2.5) / 5.)) + 1
    lon_index    = int(np.floor((dlon - 2.5) / 5.))
    lat_index    = int(np.floor((dlat - 87.5) / (-5.))) 
    # print("lpl={:}, lat idx={:} lon idx={:}".format(lons_per_lat, lat_index, lon_index))
    # print("Lat: 87.5 + {:}*(-5) = {:}".format(lat_index, 87.5+lat_index*(-5)))
    # print("Lon: 02.5 + {:}*(-5) = {:}".format(lon_index, 2.5+lon_index*(-5)))

    tl = lat_index*lons_per_lat+lon_index + 1
    tld = {}; trd = {}; bld = {}; brd = {};
    
# read respective lines off from the grid file
    with open(grid, 'r') as fin:
# read and validate first line
        line = fin.readline()
        assert line.strip().replace(" ", "") == "%latlonp:a0A1B1A2B2T:a0A1B1A2B2Q:a0A1B1A2B2dT:a0A1B1A2B2unduHsa_h:a0A1B1A2B2a_w:a0A1B1A2B2lambda:a0A1B1A2B2Tm:a0A1B1A2B2Gn_h:a0A1B1A2B2Ge_h:a0A1B1A2B2Gn_w:a0A1B1A2B2Ge_w:a0A1B1A2B2"
        #keys = line[1:].strip().split()
        for i in range(tl): 
            line = fin.readline()
        tld = parse_gpt3_line(line)
        line = fin.readline()
        trd = parse_gpt3_line(line)
        for i in range(lons_per_lat - 1): line = fin.readline()
        bld = parse_gpt3_line(line)
        line = fin.readline()
        brd = parse_gpt3_line(line)

        # print("({:}, {:})   {:}, {:})".format(tld['lat'], tld['lon'], trd['lat'], trd['lon']))
        # print("({:}, {:})   {:}, {:})".format(bld['lat'], bld['lon'], brd['lat'], brd['lon']))

# checks
    assert tld['lat'] == trd['lat'] and tld['lon']  < trd['lon']
    assert bld['lat'] == brd['lat'] and bld['lon']  < brd['lon']
    assert tld['lat']  > bld['lat'] and tld['lon'] == bld['lon']
    assert trd['lat']  > brd['lat'] and trd['lon'] == brd['lon']
    # print("{:} >= {:} and {:} < {:}".format(tld['lat'], np.degrees(lat), bld['lat'], np.degrees(lat)))
    assert (tld['lat'] >= np.degrees(lat)) and (bld['lat'] < np.degrees(lat))
    assert (tld['lon'] <= np.degrees(lon)) and (brd['lon'] > np.degrees(lon))

    return tld, trd, bld, brd

def gpt3(t, lon, lat, hgt, grid):
# first get the four surrounding ndes
    tld, trd, bld, brd = get_grid_nodes(lon, lat, grid)
# bilinear interpolation (we must first compute fractional doy)
    doy = int(t.strftime('%j'))
    try:
        t0 = datetime.datetime(t.year, t.month, t.day)
        doy += (t-t0).total_seconds() / 86400.
    except: # is this maybe an attotime instance ?
        t0 = attotime.attodatetime(t.year, t.month, t.day)
        doy += float((t-t0).total_seconds()) / 86400.

    meteo = bilinear_interpolation(doy, lat, lon, hgt, tld, trd, bld, brd)
    return meteo
