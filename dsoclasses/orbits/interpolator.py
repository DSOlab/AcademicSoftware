import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
from dsoclasses.time.calmjd import cal2fmjd
from dsoclasses.orbits.sp3c import Sp3


def flag_is_on(flag_str, flag_list):
    flag_list = [f.lower() for f in flag_list]
    for f in flag_str.strip():
        if f.lower() in flag_list:
            return True
    return False


class OrbitInterpolator:

    def __init__(
        self,
        satid,
        dct,
        interval_in_sec=1800,
        min_data_pts=4,
        itype="Polynomial",
        include_clock=False,
        owns_t=True,
        exclude_missing_clock_values=False,
        exclude_flag_events=[],
    ):
        self._satellite = satid
        self._type = itype
        self._dsec = interval_in_sec
        self._minpts = min_data_pts
        self._t = None
        self._x = []
        self._y = []
        self._z = []
        self._c = []

        # sort dictionary in chronological order
        din = dict(sorted(dct.items()))

        # get time and position in individual arrays
        if owns_t:
            for k, v in din.items():
                self._t.append(cal2fmjd(k))
                self._x.append(v["x"])
                self._y.append(v["y"])
                self._z.append(v["z"])
                self._c.append(v["c"])
                if len(self.t) >= 2:
                    assert self.t[-1] > self.t[-2]

        else:  # this is probably a call from Sp3Interpolator
            for k, v in din.items():
                try:
                    x, y, z, c, f = v[satid]
                    if (not exclude_missing_clock_values) or (
                        exclude_missing_clock_values == True and not np.isnan(c)
                    ):
                        if (
                            f == ""
                            or exclude_flag_events == []
                            or not flag_is_on(f, exclude_flag_events)
                        ):
                            self._x.append(x)
                            self._y.append(y)
                            self._z.append(z)
                            self._c.append(c)
                        else:
                            print(
                                "Skipping sp3 record for satellite {:}: event flag is {:}".format(
                                    satid, f
                                )
                            )
                    else:
                        print(
                            "Skipping sp3 record for satellite {:}: missing clock value ([{:}])".format(
                                satid, c
                            )
                        )
                except:
                    print(
                        "Error. Failed finding satellite{:} for epoch {:}".format(
                            satid, k
                        )
                    )
        # prepare interpolators if needed
        if itype != "Polynomial":
            self._last_start = -1
            self._last_stop = -1

    def create_nonpoly_interpolators(self, start, stop, tarray=None):
        tarray = self._t if tarray is None else tarray
        if self._type == "CubicSpline":
            xspl = CubicSpline(tarray[start:stop], self._x[start:stop])
            yspl = CubicSpline(tarray[start:stop], self._y[start:stop])
            zspl = CubicSpline(tarray[start:stop], self._z[start:stop])
            cspl = CubicSpline(tarray[start:stop], self._c[start:stop])
        elif self._type == "PchipInterpolator":
            xspl = PchipInterpolator(tarray[start:stop], self._x[start:stop])
            yspl = PchipInterpolator(tarray[start:stop], self._y[start:stop])
            zspl = PchipInterpolator(tarray[start:stop], self._z[start:stop])
            cspl = PchipInterpolator(tarray[start:stop], self._c[start:stop])
        return xspl, yspl, zspl, cspl

    def find_interval(self, t, tarray=None):
        """Given a (fractional MJD) array (tarray) and a datetime instance (t),
        this function will return the two indexes (of tarray) for which:
        i1: (t-tarray[i1]) [sec] <= self._dsec
        i2: (tarray[i2]-t) [sec] > self._dsec
        so that any entry in tarray[i1:i2] lies within self._dsec from t

        If the parameter tarray is set to None, the instance's self._t array
        will be used.

        If no suitable interval is found, the corresponding index is set to
        None (i.e. i1, i2, or both).
        """
        tarray = self._t if tarray is None else tarray
        mjd = cal2fmjd(t)
        start = None
        stop = None
        for idx, mjdi in enumerate(tarray):
            if (mjd - mjdi) * 86400.0 <= self._dsec:
                start = idx
                break
        for idx, mjdi in enumerate(tarray[start:]):
            if (mjdi - mjd) * 86400.0 > self._dsec:
                stop = start + idx
                break
        return start, stop

    def sat_at(self, t, tarray=None):
        """Get satellite position and clock at a given epoch (t)."""
        start, stop = self.find_interval(t, tarray)
        if start is None or stop is None or (stop - start < self._minpts):
            # print("Warning cannot interpolate at date {:}".format(t))
            raise RuntimeError(
                f"ERROR Cannot find suitable interval for interpolating satellite orbit at {t}"
            )

        mjd = cal2fmjd(t)
        if self._type == "Polynomial":
            x = np.interp(mjd, tarray[start:stop], self._x[start:stop])
            y = np.interp(mjd, tarray[start:stop], self._y[start:stop])
            z = np.interp(mjd, tarray[start:stop], self._z[start:stop])
            c = np.interp(mjd, tarray[start:stop], self._c[start:stop])
            return x, y, z, c

        # non-polynomial interpolation
        if start != self._last_start or stop != self._last_stop:
            self._last_start, self._last_stop = start, stop
            self._xspl, self._yspl, self._zspl, self._cspl = (
                self.create_nonpoly_interpolators(start, stop, tarray)
            )
        return self._xspl(mjd), self._yspl(mjd), self._zspl(mjd), self._cspl(mjd)


class Sp3Interpolator:

    def __init__(
        self,
        sp3fn,
        sat_systems,
        interval_in_sec=1800,
        min_data_pts=4,
        itype="Polynomial",
        exclude_missing_clock_values=False,
        exclude_flag_events=[],
    ):
        sp3 = Sp3(sp3fn)
        self._time_sys = sp3.time_sys
        data = sp3.get_system_pos(sat_systems, True)
        # get the date array, both in mjd an python datetime
        self._tpyd = []
        self._tmjd = []
        for k in data:
            self._tmjd.append(cal2fmjd(k))
            self._tpyd.append(k)
            if len(self._tmjd) >= 2:
                assert self._tmjd[-1] > self._tmjd[-2]
        # for each satellite, create an OrbitInterpolator instance
        self._interpolators = {}
        for sat in sp3.sat_ids:
            if sat[0].lower() in [s.lower() for s in sat_systems]:
                self._interpolators[sat] = OrbitInterpolator(
                    sat,
                    data,
                    interval_in_sec,
                    min_data_pts,
                    itype,
                    True,
                    False,
                    exclude_missing_clock_values,
                    exclude_flag_events,
                )

    def sat_at(self, satid, t, tarray=None):
        return self._interpolators[satid].sat_at(t, self._tmjd)
