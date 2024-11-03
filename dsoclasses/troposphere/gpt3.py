import numpy as np

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

def saastamoinen_zwd():
    pass

def parse_gpt3_line(line):
    data_keys = ['lat', 'lon', 'p:', 'T:', 'Q:', 'dT:', 'undu', 'Hs', 'a_h:', 'a_w:', 'lambda:', 'Tm:', 'Gn_h:', 'Ge_h:', 'Gn_w:', 'Ge_w:']
    key_coeffs = ['a0', 'A1', 'B1', 'A2', 'B2']
    data = [float(x) for x in line.split()]
    idx = 0
    dct = {}
    for key in data_keys:
        if key.endswith(':'):
            dct[key[0:-1]] = {}
            for kf in key_coeffs:
                dct[key[0:-1]][kf] = data[idx]
                idx += 1
        else:
            dct[key] = data[idx]
            idx += 1
    return dct

def interpolate(doy, lon, lat, hgt, tld, trd, bld, brd):
    cfy = np.cos(doy/365.25*2*np.pi)   # coefficient for A1
    chy = np.cos(doy/365.25*4*np.pi)   # coefficient for B1
    sfy = np.sin(doy/365.25*2*np.pi)   # coefficient for A2
    shy = np.sin(doy/365.25*4*np.pi)   # coefficient for B2
# transforming ellipsoidal height to orthometric height
    hgt = [ hgt - d['undu'] for d in [tld, trd, bld, brd] ]
# pressure, temperature at the height of the grid
    fun = lambda dct, arg: dct[arg]['a0'] + dct[arg]['A1'] * cfy + dct[arg]['B1'] * sfy + dct[arg]['A2'] * chy + dct[arg]['B2'] * cfy
    T0 = [ fun(d,'T') for d in [tld, trd, bld, brd] ]

def gpt3grid(lon, lat, grid):
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
