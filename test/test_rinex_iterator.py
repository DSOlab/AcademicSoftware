#! /usr/bin/python

import sys
from dsoclasses.rinex.gnss.rinex import GnssRinex

rnx = GnssRinex(sys.argv[1])
print(rnx.glo_slots)
sys.exit(9)

# rit = iter(rnx)
# print( next(rit).satellite('C02') )
# next(rit)

for block in rnx:
    #try:
    # print(block.satellite('C02'))
    # b = block.filter_satellite_system('gps')
    # print(b.satellite('G02'))
    for k,v in block.filter_satellite_system('glonass'):
        print(k,v)
    #except:
    #    pass 
