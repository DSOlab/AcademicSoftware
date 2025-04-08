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
        self.satellite = satid
        self.type = itype
        self.dsec = interval_in_sec
        self.minpts = min_data_pts
        self.has_clock = include_clock
        # sort dictionary in chronological order
        din = dict(sorted(dct.items()))
        # get time and position in individual arrays
        if owns_t:
            self.x = []
            self.y = []
            self.z = []
            self.t = []
            self.c = []
            for k, v in din.items():
                self.t.append(cal2fmjd(k))
                self.x.append(v["x"])
                self.y.append(v["y"])
                self.z.append(v["z"])
                self.c.append(v["c"])
                if len(self.t) >= 2:
                    assert self.t[-1] > self.t[-2]
        else:  # this is probably a call from Sp3Interpolator
            self.x = []
            self.y = []
            self.z = []
            self.c = []
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
                            self.x.append(x)
                            self.y.append(y)
                            self.z.append(z)
                            self.c.append(c)
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
            self.last_start = -1
            self.last_stop = -1

    def create_interpolators(self, start, stop, tarray=None):
        if not tarray:
            tarray = self.t
        if self.type == "CubicSpline":
            xspl = CubicSpline(tarray[start:stop], self.x[start:stop])
            yspl = CubicSpline(tarray[start:stop], self.y[start:stop])
            zspl = CubicSpline(tarray[start:stop], self.z[start:stop])
            if self.has_clock:
                cspl = CubicSpline(tarray[start:stop], self.c[start:stop])
        elif self.type == "PchipInterpolator":
            xspl = PchipInterpolator(tarray[start:stop], self.x[start:stop])
            yspl = PchipInterpolator(tarray[start:stop], self.y[start:stop])
            zspl = PchipInterpolator(tarray[start:stop], self.z[start:stop])
            if self.has_clock:
                cspl = PchipInterpolator(tarray[start:stop], self.c[start:stop])
        if self.has_clock:
            return xspl, yspl, zspl, cspl
        return xspl, yspl, zspl

    def find_interval(self, t, tarray=None):
        if not tarray:
            tarray = self.t
        mjd = cal2fmjd(t)
        start = None
        stop = None
        for idx, mjdi in enumerate(tarray):
            if (mjd - mjdi) * 86400.0 <= self.dsec:
                start = idx
                break
        for idx, mjdi in enumerate(tarray[start:]):
            if (mjdi - mjd) * 86400.0 > self.dsec:
                stop = start + idx
                break
        return start, stop

    def sat_at(self, t, tarray=None):
        start, stop = self.find_interval(t, tarray)
        if start is None or stop is None or stop - start < self.minpts:
            # print("Warning cannot interpolate at date {:}".format(t))
            raise RuntimeError
        mjd = cal2fmjd(t)
        if self.type == "Polynomial":
            x = np.interp(mjd, tarray[start:stop], self.x[start:stop])
            y = np.interp(mjd, tarray[start:stop], self.y[start:stop])
            z = np.interp(mjd, tarray[start:stop], self.z[start:stop])
            if not self.has_clock:
                return x, y, z
            c = np.interp(mjd, tarray[start:stop], self.c[start:stop])
            return x, y, z, c

        if start != self.last_start or stop != self.last_stop:
            self.last_start, self.last_stop = start, stop
            if not self.has_clock:
                self.xspl, self.yspl, self.zspl = self.create_interpolators(
                    start, stop, tarray
                )
                # print("Debug. no clocks used ...")
            else:
                self.xspl, self.yspl, self.zspl, self.cspl = self.create_interpolators(
                    start, stop, tarray
                )
        if self.has_clock:
            return self.xspl(mjd), self.yspl(mjd), self.zspl(mjd), self.cspl(mjd)
        return self.xspl(mjd), self.yspl(mjd), self.zspl(mjd), None


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
        self.time_sys = sp3.time_sys
        data = sp3.get_system_pos(sat_systems, True)
        # get the date arrayi, both in mjd an python datetime
        self.tpyd = []
        self.tmjd = []
        for k in data:
            self.tmjd.append(cal2fmjd(k))
            self.tpyd.append(k)
            if len(self.tmjd) >= 2:
                assert self.tmjd[-1] > self.tmjd[-2]
        # for each satellite, create an OrbitInterpolator instance
        self.interpolators = {}
        for sat in sp3.sat_ids:
            if sat[0].lower() in [s.lower() for s in sat_systems]:
                self.interpolators[sat] = OrbitInterpolator(
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
        return self.interpolators[satid].sat_at(t, self.tmjd, tarray)
