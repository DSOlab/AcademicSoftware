#! /usr/bin/python

import sys
import datetime

rinex_file = sys.argv[1]
observation_type = sys.argv[2]

def resolve_date(dstr):
    l = dstr.split()
    y, m, d, h, mn = [int(x) for x in l[0:5]]
    sec = float(l[5])
    s = int(sec)
    fmicrosec = (sec - s) * 1e6
    microsec  = int((sec - s) * 1e6)
    return datetime.datetime(y,m,d,h,mn,s,microsec)

# open the RINEX file for reading
with open(rinex_file, 'r') as fin:
    line = fin.readline()

    # find the line that ends with "SYS / # / OBS TYPES"
    while line and (not line.rstrip().endswith("SYS / # / OBS TYPES")):
        line = fin.readline()

    # resolve line with observable types
    if line.strip()[60:] == "SYS / # / OBS TYPES":
        # D   10  L1  L2  C1  C2  W1  W2   F   P   T   H              SYS / # / OBS TYPES 
        l = line.split()
        assert l[0] == "D"
        num_types = int(l[1])
        types = line[8:60].split()
    else:
        print("Error. Failed finding line \"SYS / # / OBS TYPES\"")
        raise RuntimeError

    # get the index of the observable we want
    index = types.index(observation_type)

    # got to the line marked as "END OF HEADER"
    while line and (not line[60:].rstrip().endswith("END OF HEADER")):
        line = fin.readline()

    line = fin.readline()
    while line:
        # read data blocks; first line is something like:
        # > 2024 01 01 00 00 40.129948870  0  3        1.740551600 0 
        assert line[0] == '>'
        epoch = resolve_date(line[2:31])
        num_beacons = int(line[34:37])
        # for every beacon, read its observations and print the one we want
        for i in range(num_beacons):
            obs_to_follow = len(types)
            lines_per_beacon = obs_to_follow // 5
            obsline = fin.readline().rstrip('\n')
            beaconid = obsline[0:3]
            for i in range(lines_per_beacon-1): obsline += fin.readline()[3:].rstrip('\n')
            # we want the observable at index
            idx = 3 + index * 16
            if obsline[idx:idx+16].strip() != "":
                value = float(obsline[idx:idx+14])
                print("{:} {:} {:.3f}".format(beaconid, epoch, value))
        # done with this block, read next block header
        line = fin.readline()
