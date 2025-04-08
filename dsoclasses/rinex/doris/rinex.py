import datetime
import math
import sys
import attotime

class DorisRinex:
    
    def resolve_date(self, dstr):
        l = dstr.split()
        y, m, d, h, mn = [int(x) for x in l[0:5]]
        sec = float(l[5])
        s = int(sec)
        fmicrosec = (sec - s) * 1e6
        microsec  = int((sec - s) * 1e6)
        fnanosec  = (fmicrosec - microsec) * 1e3
        nanosec = int(fnanosec)
        tstr = "{:4d}:{:02d}:{:02d} {:02d}:{:02d}:{:02d}".format(y, m, d, h, mn, s)
        pyt = datetime.datetime.strptime(tstr, "%Y:%m:%d %H:%M:%S")
        return attotime.attodatetime(pyt.year, pyt.month, pyt.day, pyt.hour, pyt.minute, pyt.second, microsec, nanosec)

    def parse_header(self, fn):
        self.beacons = []
        self.time_ref_beacons = []
        with open(fn, 'r') as fin:
            line  = fin.readline()
            self.version = float(line[0:9])
            self.type = line[20]
            self.system = line[40]
            assert self.system == 'D'
            line = fin.readline()
            self.obscodes = {}
            while line and line.strip() != "END OF HEADER":
                if line[60:].strip() == "MARKER NAME":
                    self.marker_name = line[0:60].strip()
                elif line[60:].strip() == "MARKER NUMBER":
                    self.marker_number = line[0:60].strip()
                elif line[60:].strip() == "MARKER TYPE":
                    self.marker_type = line[0:60].strip()
                elif line[60:].strip() == "OBSERVER / AGENCY":
                    self.observer, self.agency = line[0:20].strip(), line[20:60].strip()
                elif line[60:].strip() == "REC # / TYPE / VERS":
                    self.receiver_number, self.receiver_type, self.receiver_version = line[0:20].strip(), line[20:40].strip(), line[40:60].strip()
                elif line[60:].strip() == "ANT # / TYPE":
                    self.antenna_number, self.antenna_type = line[0:20].strip(), line[20:40].strip()
                elif line[60:].strip() == "APPROX POSITION XYZ":
                    self.xapprox, self.yapprox, self.zapprox = [ float(x) for x in line[0:60].split() ] 
                elif line[60:].strip() == "ANTENNA: DELTA H/E/N":
                    self.dh, self.de, self.dn = [ float(x) for x in line[0:60].split() ]
                elif line[60:].strip() == "SYS / # / OBS TYPES":
                    system = line[0]
                    numobsc = int(line[3:7])
                    l = line[7:60].strip().split()
                    while len(l) < numobsc:
                        line = fin.readline()
                        l += line[7:60].strip().split()
                    self.obscodes[system] = l
                elif line[60:].strip() == "TIME OF FIRST OBS":
                    self.time_first_obs = self.resolve_date(line[0:44])
                    self.time_sys = line[44:60].strip()
                elif line[60:].strip() == "TIME OF LAST OBS":
                    self.time_last_obs = self.resolve_date(line[0:44])
                    self.time_last_obs_sys = line[44:60].strip()
                elif line[60:].strip() == "SATELLITE NAME":
                    self.sat_name = line[0:60].strip()
                elif line[60:].strip() == "COSPAR NUMBER":
                    self.sat_cospar = line[0:60].strip()
                elif line[60:].strip() == "# OF STATIONS":
                    self.num_stations = int(line[0:6].strip())
                elif line[60:].strip() == "# TIME REF STATIONS":
                    self.num_ref_stations = int(line[0:6].strip())
                elif line[60:].strip() == "TIME REF STATION":
                    num   = line[0:3]
                    bias  = float(line[5:21].strip())
                    shift = float(line[20:34].strip())
                    self.time_ref_beacons += [{'num': num, 'bias': bias, 'ref_shift': shift}]
                elif line[60:].strip() == "STATION REFERENCE":
                    num   = line[0:3]
                    id    = line[5:9]
                    name  = line[10:40].strip()
                    domes = line[40:50].strip()
                    type  = int(line[51])
                    freq  = line[52:60].strip()
                    if freq == "":
                        freq = 0
                    else:
                        freq = int(freq)
                    self.beacons += [{'num': num, 'id': id, 'name': name, 'domes': domes, 'type': type, 'freqshift': freq}]
                else:
                    pass
                line = fin.readline()

    def approx_cartesian(self): return [self.xapprox, self.yapprox, self.zapprox]

    class DataBlock: 

        def __init__(self, dct):
            self.dct = dct

        def t(self): return self.dct['epoch']
        def clock_offset(self): return self.dct['clock_offset']
        def flag(self): return self.dct['flag']
        def nbeacons(self): return self.dct['num_beacons']
        def beacon(self, beaconid): return self.dct[beaconid]
        def __iter__(self): 
            return ((key, value) for key, value in self.dct.items() if key not in ['epoch', 'flag', 'num_beacons', 'clock_offset'])
    
    def __init__(self, fn):
        self.filename = fn
        self.parse_header(fn)

    def __iter__(self):
        self.stream = open(self.filename, 'r')
        line = self.stream.readline()
        while line and line.strip() != "END OF HEADER":
            line = self.stream.readline()
        return self

    def __next__(self):
        # print("called next ...")
        line = self.stream.readline()
        if not line: raise StopIteration
        assert line[0] == '>'
        data_block = {}
        data_block['epoch'] = self.resolve_date(line[2:31])
        data_block['flag'] = int(line[31:34])
        data_block['num_beacons'] = int(line[34:37])
        # print(f"basic info {data_block['epoch'].strftime("%Y-%m-%d %H:%M:%S.%f")}, {data_block['flag']}, {data_block['num_beacons']}")
        try:
            data_block['clock_offset'] = float(line[42:56])
        except:
            pass
        for i in range(data_block['num_beacons']):
            obs_to_follow = len(self.obscodes['D'])
            lines_per_beacon = obs_to_follow // 5
            obsline = self.stream.readline().rstrip('\n')
            beaconid = obsline[0:3]
            for i in range(lines_per_beacon-1): obsline += self.stream.readline()[3:].rstrip('\n')
            data_block[beaconid] = {}
            for idx, obscode in enumerate(self.obscodes['D']):
                idx = 3 + idx * 16
                if obsline[idx:idx+16].strip() != "":
                    value = float(obsline[idx:idx+14])
                    m1 = int(obsline[idx+14]) if obsline[idx+14] != " " else None
                    m2 = int(obsline[idx+15]) if obsline[idx+15] != " " else None
                    data_block[beaconid][obscode] = {'value': value, 'm1': m1, 'm2': m2}
        return self.DataBlock(data_block)
