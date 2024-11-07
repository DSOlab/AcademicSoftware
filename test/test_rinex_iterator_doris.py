#! /usr/bin/python

import sys
from dsoclasses.rinex.doris.rinex import DorisRinex

rnx = DorisRinex(sys.argv[1])
print(rnx.beacons)

for block in rnx:
    for k,v in block:
        print(k,v)
