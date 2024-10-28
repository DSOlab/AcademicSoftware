#! /usr/bin/python

import sys
from dsoclasses.rinex.gnss.rinex import GnssRinex

rnx = GnssRinex(sys.argv[1])
# rit = iter(rnx)
# print( next(rit).satellite('C02') )
# next(rit)

for block in rnx:
    try:
        print(block.satellite('C02'))
    except:
        pass 
