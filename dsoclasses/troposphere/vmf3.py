import attotime
from datetime import datetime, timezone

class SiteVmf3Block:

    def __init__(self):
        self._t = 
        self._records = {}

class SiteVmf3:

    def __init__(self, data_dir, site_list=[]):
        self._data_dir = data_dir
        srlf._site_list = site_list
        self._block0 = SiteVmf3Block()
        self._block1 = SiteVmf3Block()

    def getBlockForEpoch(t):
        

    def mf(site_list, t, els):
        for site in site_list:
            if site not in self._site_list:
                raise RuntimeError('SiteVmf3: Cannot compute vmf3 for site '+site'; not included in list!')

def _to_utc(dt: datetime) -> datetime:
    """Return dt in UTC. Naive datetimes are assumed to already be UTC."""
    if dt.tzinfo is None:
        return dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc)

def datetime_to_mjd(dt: datetime) -> float:
    """Convert a (UTC) datetime to Modified Julian Date."""
    dt = _to_utc(dt)
    y, m, d = dt.year, dt.month, dt.day
    hh = dt.hour
    mm = dt.minute
    ss = dt.second + dt.microsecond / 1e6

    # Gregorian calendar to JD (Meeus)
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + A // 4
    jd0 = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5
    frac_day = (hh + mm / 60.0 + ss / 3600.0) / 24.0
    jd = jd0 + frac_day
    return jd - 2400000.5  # MJD
"""
Python translation of vmf3.m (Vienna Mapping Function 3, 'gridless' seasonal model).

Function:
    vmf3(ah, aw, mjd, lat, lon, zd) -> (mfh, mfw)

Inputs:
    ah, aw : hydrostatic/wet continued-fraction 'a' coefficients (dimensionless), e.g. from a grid model.
    mjd    : Modified Julian Date (UTC)
    lat    : geodetic latitude [rad]
    lon    : longitude [rad]
    zd     : zenith distance [rad] (so elevation el = pi/2 - zd)

Outputs:
    mfh, mfw : hydrostatic/wet mapping functions (dimensionless)

Notes:
  - This module expects coefficient tables (anm_*, bnm_*) in 'vmf3_coeffs.npz' next to the script.
  - Coefficients are from the uploaded MATLAB source; degree/order up to n=12 (91 terms).
"""
from __future__ import annotations
import math
from pathlib import Path
from typing import Tuple, Dict
import numpy as np

DEG2RAD = math.pi / 180.0
RAD2DEG = 180.0 / math.pi

def _load_coeffs(path: str | Path | None = None) -> Dict[str, np.ndarray]:
    if path is None:
        path = Path(__file__).with_name("vmf3_coeffs.npz")
    data = np.load(path)
    return {k: data[k] for k in data.files}

def _ymd_from_mjd(mjd: float) -> Tuple[int,int,int,float]:
    """Return (year, month, day, frac_of_day)."""
    jd = mjd + 2400000.5
    # Split into integer and fractional days around noon
    Z = math.floor(jd + 0.5)
    F = (jd + 0.5) - Z
    # Gregorian calendar conversion (Fliegel & Van Flandern)
    A = Z
    alpha = math.floor((A - 1867216.25) / 36524.25)
    A = A + 1 + alpha - math.floor(alpha / 4)
    B = A + 1524
    C = math.floor((B - 122.1) / 365.25)
    D = math.floor(365.25 * C)
    E = math.floor((B - D) / 30.6001)
    day = B - D - math.floor(30.6001 * E) + F
    month = E - 1 if E < 14 else E - 13
    year = C - 4716 if month > 2 else C - 4715
    d_int = int(math.floor(day))
    frac = float(day - d_int)
    return int(year), int(month), d_int, frac

def _doy_fraction(mjd: float) -> float:
    y, m, d, frac = _ymd_from_mjd(mjd)
    leap = (y % 4 == 0 and y % 100 != 0) or (y % 400 == 0)
    mdays = [31, 29 if leap else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    doy0 = sum(mdays[:m-1]) + d
    return doy0 + frac

def _legendre_VW(x: float, y: float, z: float, nmax: int = 12):
    """Recurrence for (fully normalized) V/W as in MATLAB vmf3.m."""
    V = np.zeros((nmax+2, nmax+2), dtype=float)
    W = np.zeros_like(V)
    V[1,1] = 1.0; W[1,1] = 0.0
    V[2,1] = z*V[1,1]; W[2,1] = 0.0
    for n in range(2, nmax+1):
        V[n+1,1] = ((2*n-1)*z*V[n,1] - (n-1)*V[n-1,1]) / n
        W[n+1,1] = 0.0
    for m in range(1, nmax+1):
        V[m+1,m+1] = (2*m-1) * (x*V[m,m] - y*W[m,m])
        W[m+1,m+1] = (2*m-1) * (x*W[m,m] + y*V[m,m])
        if m < nmax:
            V[m+2,m+1] = (2*m+1) * z * V[m+1,m+1]
            W[m+2,m+1] = (2*m+1) * z * W[m+1,m+1]
        for n in range(m+2, nmax+1):
            V[n+1,m+1] = ((2*n-1)*z*V[n,m+1] - (n+m-1)*V[n-1,m+1]) / (n-m)
            W[n+1,m+1] = ((2*n-1)*z*W[n,m+1] - (n+m-1)*W[n-1,m+1]) / (n-m)
    return V, W

def _seasonal_series(V: np.ndarray, W: np.ndarray, anm: np.ndarray, bnm: np.ndarray):
    """
    Combine spherical-harmonic coefficients into seasonal series:
    returns (A0, A1, B1, A2, B2) for one field.
    anm/bnm are arrays of shape (91,5), ordered by n=0..12, m=0..n.
    """
    nmax = 12
    A0=A1=B1=A2=B2 = 0.0
    i = 0
    for n in range(0, nmax+1):
        for m in range(0, n+1):
            vn = V[n+1, m+1]
            wn = W[n+1, m+1]
            a = anm[i]; b = bnm[i]
            A0 += a[0]*vn + b[0]*wn
            A1 += a[1]*vn + b[1]*wn
            B1 += a[2]*vn + b[2]*wn
            A2 += a[3]*vn + b[3]*wn
            B2 += a[4]*vn + b[4]*wn
            i += 1
    return A0, A1, B1, A2, B2

def vmf3(ah: float, aw: float, mjd: float, lat: float, lon: float, zd: float,
         coeffs_path: str | Path | None = None) -> Tuple[float, float]:
    """
    Compute VMF3 hydrostatic/wet mapping functions at a site and epoch.

    Parameters
    ----------
    ah, aw : float
        Continued-fraction 'a' parameters (hydrostatic, wet).
    mjd : float
        Modified Julian Date (UTC).
    lat, lon : float
        Latitude and longitude in radians.
    zd : float
        Zenith distance in radians.

    Returns
    -------
    mfh, mfw : floats
        Hydrostatic and wet mapping functions.
    """
    coeffs = _load_coeffs(coeffs_path)
    # Geometry
    el = 0.5*math.pi - zd
    polDist = 0.5*math.pi - lat  # co-latitude
    x = math.sin(polDist) * math.cos(lon)
    y = math.sin(polDist) * math.sin(lon)
    z = math.cos(polDist)
    V, W = _legendre_VW(x, y, z, nmax=12)

    # Seasonal series for b and c
    bh_A0,bh_A1,bh_B1,bh_A2,bh_B2 = _seasonal_series(V, W, coeffs['anm_bh'], coeffs['bnm_bh'])
    bw_A0,bw_A1,bw_B1,bw_A2,bw_B2 = _seasonal_series(V, W, coeffs['anm_bw'], coeffs['bnm_bw'])
    ch_A0,ch_A1,ch_B1,ch_A2,ch_B2 = _seasonal_series(V, W, coeffs['anm_ch'], coeffs['bnm_ch'])
    cw_A0,cw_A1,cw_B1,cw_A2,cw_B2 = _seasonal_series(V, W, coeffs['anm_cw'], coeffs['bnm_cw'])

    # Time dependence via day-of-year (fractional)
    doy = _doy_fraction(mjd)
    w = 2.0*math.pi * (doy/365.25)
    bh = bh_A0 + bh_A1*math.cos(w) + bh_B1*math.sin(w) + bh_A2*math.cos(2*w) + bh_B2*math.sin(2*w)
    bw = bw_A0 + bw_A1*math.cos(w) + bw_B1*math.sin(w) + bw_A2*math.cos(2*w) + bw_B2*math.sin(2*w)
    ch = ch_A0 + ch_A1*math.cos(w) + ch_B1*math.sin(w) + ch_A2*math.cos(2*w) + ch_B2*math.sin(2*w)
    cw = cw_A0 + cw_A1*math.cos(w) + cw_B1*math.sin(w) + cw_A2*math.cos(2*w) + cw_B2*math.sin(2*w)

    # Continued-fraction mapping (Kouba style)
    s = math.sin(el)
    # hydrostatic
    mfh = (1.0 + ah/(1.0 + bh/(1.0 + ch))) / (s + ah/(s + bh/(s + ch)))
    # wet
    mfw = (1.0 + aw/(1.0 + bw/(1.0 + cw))) / (s + aw/(s + bw/(s + cw)))
    return mfh, mfw
