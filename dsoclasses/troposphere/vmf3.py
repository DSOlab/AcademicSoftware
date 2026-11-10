"""VMF3 site interpolation and global Vienna mapping function implementation.

The file reader holds a small table of atmospheric records at known epochs
rather than maintaining fragile mutable stream positions. This deliberately
trades memory for straightforward validation and reproducibility.
"""
from __future__ import annotations

from bisect import bisect_left
from dsoclasses.time.calmjd import cal2fmjd, mjd2cal
from dsoclasses.time.pyattotime import to_attoseconds
import datetime
import math
import warnings
from pathlib import Path
from typing import Tuple, Dict
import numpy as np

_SITE_FIELDS = ('ah', 'aw', 'zhd', 'zwd', 'p', 't', 'wvp')


class SiteVmf3:
    """Interpolated VMF3 site coefficients from a standard 6-hourly file.

    Values between consecutive epochs are linearly interpolated; out-of-range
    queries are rejected. Missing requested sites and duplicate site/epoch
    records are errors, not silently substituted by neighboring stations.
    """

    def __init__(self, filename, site_list=None):
        self._fn = str(filename)
        self._site_list = set(site_list or [])
        records = {}
        # A small number of records in genuine site-wise VMF3 products have
        # nonpositive wet mapping coefficients (a_w).  Such a coefficient is
        # not usable by the VMF3 mapping function, but one bad six-hour record
        # must not make an otherwise valid annual file impossible to open.
        # Remember these records so we can tell the caller what was skipped.
        invalid_aw = []
        with open(filename, encoding='ascii', errors='replace') as stream:
            for line_no, line in enumerate(stream, start=1):
                if not line.strip() or line.lstrip().startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 9:
                    raise ValueError(f'Malformed VMF3 record at {filename}:{line_no}')
                name = parts[0]
                if self._site_list and name not in self._site_list:
                    continue
                try:
                    epoch = mjd2cal(float(parts[1]))
                    values = tuple(float(x) for x in parts[2:9])
                except (ValueError, OverflowError) as exc:
                    raise ValueError(f'Invalid VMF3 values at line {line_no}') from exc
                if not np.all(np.isfinite(values)):
                    raise ValueError(f'Non-finite VMF3 values at line {line_no}')
                ah, aw, zhd, zwd, pressure, temperature, vapour = values
                # Negative temperatures (degrees Celsius) are entirely valid.
                # The other quantities must satisfy their physical sign
                # constraints.  A negative ZHD or pressure suggests damaged
                # input and should not be quietly hidden by the reader.
                if min(ah, zhd, pressure) <= 0 or zwd < 0 or vapour < 0:
                    raise ValueError(f'Unphysical VMF3 fields at line {line_no}')
                if aw <= 0:
                    # Real VMF3 site products occasionally include anomalous
                    # negative a_w values.  They cannot safely be used in the
                    # mapping function.  Skip the *whole* station/epoch record
                    # instead of changing a_w, substituting another site's
                    # data, or losing the entire annual dataset.
                    invalid_aw.append((line_no, name, epoch, aw))
                    continue
                slot = records.setdefault(epoch, {})
                if name in slot:
                    raise ValueError(f'Duplicate VMF3 data for {name} at {epoch}')
                slot[name] = dict(zip(_SITE_FIELDS, values))
        if invalid_aw:
            first_line, first_site, first_epoch, first_aw = invalid_aw[0]
            warnings.warn(
                f'Skipped {len(invalid_aw)} VMF3 record(s) with nonpositive '
                f'wet mapping coefficient a_w in {filename}; first: '
                f'{first_site} at {first_epoch} (line {first_line}, '
                f'a_w={first_aw:g}). Interpolation across missing six-hour '
                f'records is not permitted.',
                RuntimeWarning,
                stacklevel=2,
            )
        if not records:
            if self._site_list:
                raise LookupError(f'VMF3 stations not present in file: {sorted(self._site_list)}')
            raise ValueError(f'No usable VMF3 data in {filename}')
        self._records = records
        self._epochs = sorted(records)
        # Normalized times are used only for searching, never as floating
        # Julian dates; queries may be native datetime or attodatetime.
        self._epoch_as = tuple(to_attoseconds(t) for t in self._epochs)
        if self._site_list:
            missing = self._site_list - {site for d in records.values() for site in d}
            if missing:
                raise LookupError(f'VMF3 stations not present in file: {sorted(missing)}')

    def _bracket(self, t):
        """Return a bracketing epoch pair (or same epoch for an exact hit)."""
        query = to_attoseconds(t)
        if query < self._epoch_as[0] or query > self._epoch_as[-1]:
            raise ValueError(f'VMF3 epoch {t} is outside file coverage '
                             f'[{self._epochs[0]}, {self._epochs[-1]}]')
        idx = bisect_left(self._epoch_as, query)
        if idx < len(self._epoch_as) and self._epoch_as[idx] == query:
            return self._epochs[idx], self._epochs[idx]
        return self._epochs[idx - 1], self._epochs[idx]

    def interpolate_site(self, site, t):
        """Return time-interpolated atmospheric fields for one station."""
        if self._site_list and site not in self._site_list:
            raise ValueError(f'Site {site!r} was not selected when loading VMF3')
        before, after = self._bracket(t)
        try:
            d0 = self._records[before][site]
            d1 = self._records[after][site]
        except KeyError as exc:
            raise LookupError(f'Missing VMF3 data for site {site!r} at a bracketing epoch') from exc
        if before == after:
            return d0.copy()
        # Site VMF3 products are normally six-hourly. Do not interpolate
        # silently across missing days, which would fabricate atmospheric
        # conditions unrelated to the observation epoch.
        if to_attoseconds(after) - to_attoseconds(before) > 6 * 3600 * 10**18 + 10**12:
            raise ValueError(f'VMF3 gap for {site!r} exceeds six hours')
        # Compare attodatetime to datetime without mixing concrete types.
        alpha = ((to_attoseconds(t) - to_attoseconds(before))
                 / (to_attoseconds(after) - to_attoseconds(before)))
        return {key: (1 - alpha) * d0[key] + alpha * d1[key] for key in _SITE_FIELDS}

    def getBlocksForEpoch(self, t):
        """Compatibility method returning the two bracketing datetimes."""
        return self._bracket(t)

    def vmf3(self, site_name, site_lat, site_lon, el, t):
        """Return interpolated site data plus hydrostatic/wet mapping factors."""
        if not all(math.isfinite(float(v)) for v in (site_lat, site_lon, el)):
            raise ValueError('VMF3 latitude, longitude and elevation must be finite')
        if not -math.pi/2 <= site_lat <= math.pi/2:
            raise ValueError('VMF3 latitude must be in -90..90 degrees')
        if not (0 < el <= math.pi / 2):
            raise ValueError('VMF3 requires elevation 0 < el <= pi/2 radians')
        fields = self.interpolate_site(site_name, t)
        fields['mfh'], fields['mfw'] = vmf3(
            fields['ah'], fields['aw'], cal2fmjd(t),
            site_lat, site_lon, math.pi / 2 - el,
        )
        return fields

    def tropo_delay(self, site_name, site_lat, site_lon, el, t):
        d = self.vmf3(site_name, site_lat, site_lon, el, t)
        return d['zhd'] * d['mfh'] + d['zwd'] * d['mfw']

    def tropo_dealy_info(self, site_name, site_lat, site_lon, el, t):
        """Backwards-compatible historical misspelling of ``vmf3``."""
        return self.vmf3(site_name, site_lat, site_lon, el, t)

def _legendre_VW(x: float, y: float, z: float, nmax: int = 12):
    """Recurrence for (fully normalized) V/W as in MATLAB vmf3.m."""
    V = np.zeros((nmax + 2, nmax + 2), dtype=float)
    W = np.zeros_like(V)
    V[1, 1] = 1.0
    W[1, 1] = 0.0
    V[2, 1] = z * V[1, 1]
    W[2, 1] = 0.0
    for n in range(2, nmax + 1):
        V[n + 1, 1] = ((2 * n - 1) * z * V[n, 1] - (n - 1) * V[n - 1, 1]) / n
        W[n + 1, 1] = 0.0
    for m in range(1, nmax + 1):
        V[m + 1, m + 1] = (2 * m - 1) * (x * V[m, m] - y * W[m, m])
        W[m + 1, m + 1] = (2 * m - 1) * (x * W[m, m] + y * V[m, m])
        if m < nmax:
            V[m + 2, m + 1] = (2 * m + 1) * z * V[m + 1, m + 1]
            W[m + 2, m + 1] = (2 * m + 1) * z * W[m + 1, m + 1]
        for n in range(m + 2, nmax + 1):
            V[n + 1, m + 1] = (
                (2 * n - 1) * z * V[n, m + 1] - (n + m - 1) * V[n - 1, m + 1]
            ) / (n - m)
            W[n + 1, m + 1] = (
                (2 * n - 1) * z * W[n, m + 1] - (n + m - 1) * W[n - 1, m + 1]
            ) / (n - m)
    return V, W


def _seasonal_series(V: np.ndarray, W: np.ndarray, anm: np.ndarray, bnm: np.ndarray):
    """
    Combine spherical-harmonic coefficients into seasonal series:
    returns (A0, A1, B1, A2, B2) for one field.
    anm/bnm are arrays of shape (91,5), ordered by n=0..12, m=0..n.
    """
    nmax = 12
    A0 = A1 = B1 = A2 = B2 = 0.0
    i = 0
    for n in range(0, nmax + 1):
        for m in range(0, n + 1):
            vn = V[n + 1, m + 1]
            wn = W[n + 1, m + 1]
            a = anm[i]
            b = bnm[i]
            A0 += a[0] * vn + b[0] * wn
            A1 += a[1] * vn + b[1] * wn
            B1 += a[2] * vn + b[2] * wn
            A2 += a[3] * vn + b[3] * wn
            B2 += a[4] * vn + b[4] * wn
            i += 1
    return A0, A1, B1, A2, B2


def _load_coeffs(path: str | Path | None = None) -> Dict[str, np.ndarray]:
    if path is None:
        path = Path(__file__).with_name("vmf3_coeffs.npz")
    data = np.load(path)
    return {k: data[k] for k in data.files}


def _ymd_from_mjd(mjd: float) -> Tuple[int, int, int, float]:
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
    doy0 = sum(mdays[: m - 1]) + d
    return doy0 + frac


def vmf3(
    ah: float,
    aw: float,
    mjd: float,
    lat: float,
    lon: float,
    zd: float,
    coeffs_path: str | Path | None = None,
) -> Tuple[float, float]:
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
    if not all(math.isfinite(float(v)) for v in (ah, aw, mjd, lat, lon, zd)):
        raise ValueError('VMF3 inputs must be finite')
    if not -math.pi/2 <= lat <= math.pi/2:
        raise ValueError('VMF3 latitude must be in -90..90 degrees')
    if ah <= 0 or aw <= 0 or not 0 <= zd < math.pi/2:
        raise ValueError('VMF3 requires positive a coefficients and above-horizon geometry')
    coeffs = _load_coeffs(coeffs_path)
    # Geometry
    el = 0.5 * math.pi - zd
    polDist = 0.5 * math.pi - lat  # co-latitude
    x = math.sin(polDist) * math.cos(lon)
    y = math.sin(polDist) * math.sin(lon)
    z = math.cos(polDist)
    V, W = _legendre_VW(x, y, z, nmax=12)

    # Seasonal series for b and c
    bh_A0, bh_A1, bh_B1, bh_A2, bh_B2 = _seasonal_series(
        V, W, coeffs["anm_bh"], coeffs["bnm_bh"]
    )
    bw_A0, bw_A1, bw_B1, bw_A2, bw_B2 = _seasonal_series(
        V, W, coeffs["anm_bw"], coeffs["bnm_bw"]
    )
    ch_A0, ch_A1, ch_B1, ch_A2, ch_B2 = _seasonal_series(
        V, W, coeffs["anm_ch"], coeffs["bnm_ch"]
    )
    cw_A0, cw_A1, cw_B1, cw_A2, cw_B2 = _seasonal_series(
        V, W, coeffs["anm_cw"], coeffs["bnm_cw"]
    )

    # Time dependence via day-of-year (fractional)
    doy = _doy_fraction(mjd)
    w = 2.0 * math.pi * (doy / 365.25)
    bh = (
        bh_A0
        + bh_A1 * math.cos(w)
        + bh_B1 * math.sin(w)
        + bh_A2 * math.cos(2 * w)
        + bh_B2 * math.sin(2 * w)
    )
    bw = (
        bw_A0
        + bw_A1 * math.cos(w)
        + bw_B1 * math.sin(w)
        + bw_A2 * math.cos(2 * w)
        + bw_B2 * math.sin(2 * w)
    )
    ch = (
        ch_A0
        + ch_A1 * math.cos(w)
        + ch_B1 * math.sin(w)
        + ch_A2 * math.cos(2 * w)
        + ch_B2 * math.sin(2 * w)
    )
    cw = (
        cw_A0
        + cw_A1 * math.cos(w)
        + cw_B1 * math.sin(w)
        + cw_A2 * math.cos(2 * w)
        + cw_B2 * math.sin(2 * w)
    )

    # Continued-fraction mapping (Kouba style)
    s = math.sin(el)
    # hydrostatic
    mfh = (1.0 + ah / (1.0 + bh / (1.0 + ch))) / (s + ah / (s + bh / (s + ch)))
    # wet
    mfw = (1.0 + aw / (1.0 + bw / (1.0 + cw))) / (s + aw / (s + bw / (s + cw)))
    if not math.isfinite(float(mfh)) or not math.isfinite(float(mfw)) or min(mfh, mfw) < 1 - 1e-12:
        raise ValueError('VMF3 returned nonphysical mapping factors')
    return mfh, mfw
