#! /usr/bin/python

import sys
import datetime
from dsoclasses.itrf_realization import ItrfSiteMotion

s = ItrfSiteMotion('ARFB', sys.argv[1])
x, y, z = s.at(datetime.datetime.now())
print("ARFB (noPsd) {:+.4f} {:+.4f} {:+.4f}".format(x, y, z))

s = ItrfSiteMotion('ARFB', sys.argv[1], sys.argv[2])
x, y, z = s.at(datetime.datetime.now())
print("ARFB (w Psd) {:+.4f} {:+.4f} {:+.4f}".format(x, y, z))

s = ItrfSiteMotion('ARFB', sys.argv[1], sys.argv[2], sys.argv[3])
x, y, z = s.at(datetime.datetime.now())
print("ARFB (w Psd) {:+.4f} {:+.4f} {:+.4f}".format(x, y, z))

tl = []; xl = []; yl = []; zl = [];
t = datetime.datetime(2006, 1, 1)
while t<datetime.datetime.now():
    tx, ty, tz = s.at(datetime.datetime.now())
    tl.append(t); xl.append(x);yl.append(y);zl.append(z)
    t += datetime.timedelta(days=1)
