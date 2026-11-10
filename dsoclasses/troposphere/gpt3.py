import numpy as np
import datetime

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
    x  = np.degrees(lon) % 360.0
    if x < x1:
        x += 360.0
    if not (x1 - 1e-9 <= x <= x2 + 1e-9):
        raise ValueError("Longitude is outside selected GPT3 interpolation cell")
    if not (y1 - 1e-9 <= y <= y2 + 1e-9):
        raise ValueError("Latitude is outside selected GPT3 interpolation cell")
    bi = lambda Q12,Q22,Q11,Q21: (1./((x2-x1)*(y2-y1))) * (Q11*(x2-x)*(y2-y) + Q21*(x-x1)*(y2-y) + Q12*(x2-x)*(y-y1) + Q22*(x-x1)*(y-y1))
    return {'p': bi(*pl), 'T': bi(*Tl), 'dT': bi(*dTl)*1e3, 'e': bi(*el), 'ah': bi(*ahl), 'aw': bi(*awl), 'la': bi(*lal), 'Tm': bi(*Tml), 'Gn_h': bi(*Gn_hl), 'Ge_h': bi(*Ge_hl), 'Gn_w': bi(*Gn_wl), 'Ge_w': bi(*Ge_wl)}


def get_grid_nodes(lon, lat, grid):
    """Read the four grid points bracketing a geodetic position.

    Longitude is periodic: the 357.5° / 2.5° seam must work just like any
    interior cell. The published 5-degree grid contains nodes from -87.5°
    to 87.5° latitude; higher latitudes are rejected rather than silently
    extrapolating across the poles.

    This intentionally scans the text file on each call; clarity and strict
    grid validation are prioritized over speed.
    """
    import math

    if not math.isfinite(float(lon)) or not math.isfinite(float(lat)):
        raise ValueError('GPT3 latitude and longitude must be finite')
    longitude = math.degrees(lon) % 360.0
    latitude = math.degrees(lat)
    if abs(latitude) > 87.5 + 1e-9:
        raise ValueError('GPT3 5-degree grid does not cover latitudes beyond +/-87.5 degrees')

    nodes = {}
    with open(grid, encoding='ascii', errors='replace') as stream:
        header = stream.readline()
        if not header.lstrip().startswith('%lat lon') and not header.lstrip().startswith('% lat lon'):
            # Most published grids start with the canonical "% lat lon ..." header.
            # There is no need to compare every column token here; validate
            # each required numeric record below instead.
            if 'lat' not in header.lower() or 'lon' not in header.lower():
                raise ValueError('GPT3 grid header is missing latitude/longitude columns')
        for line_no, line in enumerate(stream, start=2):
            if not line.strip() or line.startswith('%'):
                continue
            try:
                node = parse_gpt3_line(line)
            except (ValueError, IndexError) as exc:
                raise ValueError(f'Malformed GPT3 grid row at line {line_no}') from exc
            key = (round(node['lat'], 6), round(node['lon'] % 360.0, 6))
            if key in nodes:
                raise ValueError(f'Duplicate GPT3 grid point {key}')
            nodes[key] = node

    if not nodes:
        raise ValueError('GPT3 grid has no data rows')

    # Compute longitude cell at 2.5, 7.5, ..., 357.5 degrees. The upper
    # boundary is allowed beyond 360 so bilinear interpolation remains simple.
    lon_lower = 2.5 + 5.0 * math.floor((longitude - 2.5) / 5.0)
    lon_upper = lon_lower + 5.0
    lon_query = longitude
    if longitude < 2.5:
        lon_query += 360.0
    # Compute latitude cell with a well-defined choice at exact boundaries.
    lat_lower = -87.5 + 5.0 * math.floor((latitude + 87.5) / 5.0)
    if lat_lower >= 87.5:
        lat_lower = 82.5
    lat_upper = lat_lower + 5.0

    def lookup(la, lo):
        key = (round(la, 6), round(lo % 360.0, 6))
        try:
            item = nodes[key].copy()
        except KeyError as exc:
            raise ValueError(f'GPT3 grid lacks node {key}; check grid resolution and completeness') from exc
        item['lon'] = lo  # Unwrap seam for the subsequent bilinear formula.
        return item

    # Returned order used by bilinear_interpolation: top-left, top-right,
    # bottom-left, bottom-right.
    tld = lookup(lat_upper, lon_lower)
    trd = lookup(lat_upper, lon_upper)
    bld = lookup(lat_lower, lon_lower)
    brd = lookup(lat_lower, lon_upper)
    return tld, trd, bld, brd


def gpt3(t, lon, lat, hgt, grid):
    """Return GPT3 meteorology at an epoch, geodetic position and grid path.

    ``lon`` and ``lat`` are radians, ``hgt`` is ellipsoidal metres. Epochs
    may be datetime or attotime; no datetime-type-specific subtraction is
    necessary to compute the fractional day-of-year.
    """
    tld, trd, bld, brd = get_grid_nodes(lon, lat, grid)
    # Use the same epoch normalization as the orbit and RINEX modules.
    # MJD and day-of-year are calendar quantities; no implicit GPS/UTC
    # conversion takes place, and sub-microsecond RINEX labels are accepted.
    from dsoclasses.time.calmjd import cal2mjd
    from dsoclasses.time.pyattotime import at2pt
    _, seconds_of_day = cal2mjd(t)
    doy = at2pt(t).timetuple().tm_yday + seconds_of_day / 86400.0
    return bilinear_interpolation(doy, lat, lon, hgt, tld, trd, bld, brd)
