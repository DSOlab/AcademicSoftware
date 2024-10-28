import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
from dsoclasses.time.calmjd import cal2fmjd

class OrbitInterpolator:

    def __init__(self, satid, dct, interval_in_sec=1800, min_data_pts=4, itype='Polynomial'):
        self.satellite = satid
        self.type = itype
        self.dsec = interval_in_sec
        self.minpts = min_data_pts
# sort dictionary in chronological order
        din = dict(sorted(dct.items()))
# get time and position in individual arrays
        self.x=[];self.y=[];self.z=[];self.t=[];
        for k, v in din.items():
            self.t.append(cal2fmjd(k))
            self.x.append(v['x']); self.y.append(v['y']); self.z.append(v['z']);
            if len(self.t) >= 2: assert self.t[-1] > self.t[-2]
# prepare interpolators if needed
        if itype != 'Polynomial':
            self.last_start = 0
            self.last_stop = min_data_pts
            self.xspl, self.yspl, self.zspl = self.create_interpolators(0, min_data_pts)

    def create_interpolators(self, start, stop):
        if self.type == 'CubicSpline':
            xspl = CubicSpline(self.t[start:stop], self.x[start:stop])
            yspl = CubicSpline(self.t[start:stop], self.y[start:stop])
            zspl = CubicSpline(self.t[start:stop], self.z[start:stop])
        elif self.type == 'PchipInterpolator':
            xspl = PchipInterpolator(self.t[start:stop], self.x[start:stop])
            yspl = PchipInterpolator(self.t[start:stop], self.y[start:stop])
            zspl = PchipInterpolator(self.t[start:stop], self.z[start:stop])
        return xspl, yspl, zspl

    def find_interval(self, t):
        mjd = cal2fmjd(t)
        start = None
        stop = None
        for idx, mjdi in enumerate(self.t):
            if (mjd-mjdi)*86400. <= self.dsec:
                start = idx
                break
        for idx, mjdi in enumerate(self.t[start:]):
            if (mjdi-mjd)*86400. > self.dsec:
                stop = start + idx
                break
        return start, stop

    def at(self, t):
        start, stop = self.find_interval(t)
        if start is None or stop is None or stop-start<self.minpts:
            # print("Warning cannot interpolate at date {:}".format(t))
            raise RuntimeError
        mjd = cal2fmjd(t)
        if self.type == 'Polynomial':
            x = np.interp(mjd, self.t[start:stop], self.x[start:stop])
            y = np.interp(mjd, self.t[start:stop], self.y[start:stop])
            z = np.interp(mjd, self.t[start:stop], self.z[start:stop])
            return x,y,z
        
        if start != self.last_start or stop != last_stop:
            self.xspl, self.yspl, self.zspl = self.create_interpolators(start, stop)
            self.last_start, self.last_stop = start, stop
        return self.xspl(mjd), self.yspl(mjd), self.zspl(mjd)
        
        #elif self.type == 'CubicSpline': 
        #    if start != self.last_start or stop != last_stop:
        #        self.xspl = CubicSpline(self.t[start:stop], self.x[start:stop])
        #        self.yspl = CubicSpline(self.t[start:stop], self.y[start:stop])
        #        self.zspl = CubicSpline(self.t[start:stop], self.z[start:stop])
        #        self.last_start, self.last_stop = start, stop
        #elif self.type == 'CubicSpline': 
        #    if start != self.last_start or stop != last_stop:
        #        self.xspl = CubicSpline(self.t[start:stop], self.x[start:stop])
        #        self.yspl = CubicSpline(self.t[start:stop], self.y[start:stop])
        #        self.zspl = CubicSpline(self.t[start:stop], self.z[start:stop])
        #        self.last_start, self.last_stop = start, stop
        #    x, y, z = self.xspl(mjd), self.yspl(mjd), self.zspl(mjd)
        #return x,y,z
