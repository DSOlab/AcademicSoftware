from dsoclasses.time.calmjd import cal2fmjd, mjd2cal
import datetime
import math
from pathlib import Path
from typing import Tuple, Dict
import numpy as np


class DummyStream:
    def __init__(self, fn):
        # some comments contain on-ascii chars ...
        self._stream = open(fn, "r", encoding="ascii", errors="ignore")

    def stream(self):
        return self._stream

    def stream_pos(self):
        return self._stream.tell()

    def set_stream_pos(self, pos):
        self._stream.seek(pos)

    def rewind(self, pos=0):
        self._stream.seek(pos)


class SiteVmf3Block:

    def __init__(self):
        self._t = datetime.datetime.min
        self._records = {}

    def site_data(self, site):
        return self._records[site]

    @staticmethod
    def queryNextBlockEpoch(fin: DummyStream) -> datetime.datetime:
        pos = fin.stream_pos()
        while True:
            line = fin.stream().readline()
            if not line:  # EOF
                raise EOFError("Failed quering next block epoch; EOF encountered")
            if line[0] != "#":
                l = line.split()
                return mjd2cal(float(l[1]))
        fin.set_stream_pos(pos)

    @staticmethod
    def skipNextBlock(fin: DummyStream):
        current_epoch = -999.0
        # self._records = {}
        pos = fin.stream_pos()
        while True:
            line = fin.stream().readline()
            if not line:  # EOF
                raise EOFError("Failed skipping next block epoch; EOF encountered")
            if line[0] != "#":
                mjd = float(line.split()[1])
                current_epoch = mjd if current_epoch == -999.0 else current_epoch
                if mjd == current_epoch:
                    pos = fin.stream_pos()
                elif mjd > current_epoch:
                    break
        # self._t = mjd2cal(current_epoch)
        fin.set_stream_pos(pos)

    def parseNextBlock(self, fin: DummyStream, site_list=[]):
        """
        # columns:
        # -------
        # (1) station name
        # (2) modified Julian date
        # (3) hydrostatic mf coefficient a_h
        # (4) wet mf coefficient a_w
        # (5) zenith hydrostatic delay (m)
        # (6) zenith wet delay (m)
        # (7) pressure at the site (hPa)
        # (8) temperature at the site (�C)
        # (9) water vapour pressure at the site (hPa)
        #
        #    ADEA      58119.00  0.00121504  0.00074311  2.2149  0.0364   973.84  -0.21   2.28
        """
        current_epoch = -999.0
        self._records = {}
        pos = fin.stream_pos()
        while True:
            line = fin.stream().readline()
            if not line:  # EOF
                raise EOFError("Failed parsing next block epoch; EOF encountered")
            l = line.split()
            if line[0] != "#":
                if current_epoch == -999.0 or float(l[1]) == current_epoch:
                    current_epoch = float(l[1])
                    if l[0] in site_list or site_list == []:
                        self._records[l[0]] = {
                            "ah": float(l[2]),
                            "aw": float(l[3]),
                            "zhd": float(l[4]),
                            "zwd": float(l[5]),
                            "p": float(l[6]),
                            "t": float(l[7]),
                            "wvp": float(l[8]),
                        }
                    pos = fin.stream_pos()
                elif float(l[1]) > current_epoch:
                    break
        self._t = mjd2cal(current_epoch)
        fin.set_stream_pos(pos)


class SiteVmf3:

    def __init__(self, fn, site_list=[]):
        self._fn = fn
        self._fin = DummyStream(self._fn)
        self._site_list = site_list
        self._block0 = SiteVmf3Block()
        self._block1 = SiteVmf3Block()

    def interpolate_site(self, site, t):
        data0 = self._block0.site_data(site)
        data1 = self._block1.site_data(site)
        x1mx0 = (self._block1._t - self._block0._t).total_seconds()
        xmx0 = (t - self._block0._t).total_seconds()
        x1mx = (self._block1._t - t).total_seconds()
        w0 = x1mx / x1mx0
        w1 = xmx0 / x1mx0
        return {
            k: data0[k] * w0 + data1[k] * w1
            for k in ["ah", "aw", "zhd", "zwd", "p", "t", "wvp"]
        }

    def getBlocksForEpoch(self, t):
        # quick return
        if t >= self._block0._t and t < self._block1._t:
            return

        # for each day we have four data blocks, starting at [mjd.0, mjd.25, mjd.50, mjd.75]
        mjd = cal2fmjd(t)
        ip = int((mjd - int(mjd)) * 100) // 25  # previous range in [0,3]
        np = ip + 1 if (ip != 3) else 0
        t0 = mjd2cal(int(mjd) + (ip * 25) / 100.0)
        t1 = mjd2cal(int(mjd + int(np == 0)) + (np * 25) / 100.0)

        # easy case, we just move forward one epoch ...
        if self._block1._t == t0:
            self._block0._t = self._block1._t
            self._block1.parseNextBlock(self._fin, self._site_list)
            assert self._block1._t == t1 and self._block0._t == t0
            return

        # else rewind file and start looking for required epochs ...
        self._fin.rewind()
        next_epoch = SiteVmf3Block.queryNextBlockEpoch(self._fin)
        while next_epoch < t0:
            SiteVmf3Block.skipNextBlock(self._fin)
            next_epoch = SiteVmf3Block.queryNextBlockEpoch(self._fin)

        self._block0.parseNextBlock(self._fin, self._site_list)
        self._block1.parseNextBlock(self._fin, self._site_list)
        assert self._block1._t == t1 and self._block0._t == t0
        return

    def vmf3(self, site_name, site_lat, site_lon, el, t):
        if (site_name not in self._site_list) and (self._site_list != []):
            raise RuntimeError(
                f"Site {site_name} not included in SiteVmf3 instance! cannot compute vmf3."
            )
        self.getBlocksForEpoch(t)
        site_dct = self.interpolate_site(site_name, t)
        mfh, mfw = vmf3(
            site_dct["ah"],
            site_dct["aw"],
            cal2fmjd(t),
            site_lat,
            site_lon,
            math.pi / 2 - el,
        )
        site_dct["mfh"] = mfh
        site_dct["mfw"] = mfw
        return site_dct

    def tropo_delay(self, site_name, site_lat, site_lon, el, t):
        d = self.vmf3(site_name, site_lat, site_lon, el, t)
        return d["zhd"] * d["mfh"] + d["zwd"] * d["mfw"]

    def tropo_dealy_info(self, site_name, site_lat, site_lon, el, t):
        d = self.vmf3(site_name, site_lat, site_lon, el, t)
        return d


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
    return mfh, mfw
