import datetime
import attotime
import math
import sys

def sysstr2sysid(sat_sys):
    if   sat_sys.strip().lower() in ['gps', 'g']: return 'G'
    elif sat_sys.strip().lower() in ['glonass', 'glo', 'r']: return 'R'
    elif sat_sys.strip().lower() in ['galileo', 'gal', 'e']: return 'E'
    elif sat_sys.strip().lower() in ['beidou', 'bds', 'c']: return 'C'
    elif sat_sys.strip().lower() in ['qzss', 'j']: return 'J'
    elif sat_sys.strip().lower() in ['irnss', 'navic/irnss', 'i']: return 'I'
    elif sat_sys.strip().lower() in ['sbas', 's']: return 'S'
    else:
        print("Error. Unknown satellite system string: {:}".format(sat_sys), file=sys.stderr)
        raise RuntimeError

def fetch(dct, *args):
    """ Given a dictionary containing e.g.
        R05 : {'C1C': {'value': 23539032.631, 'lli': None, 'ssi': 6}, 'L1C': {'value': 125829717.51, 'lli': 0, 'ssi': 6}, 'D1C': {'value': -4149.772, 'lli': None, 'ssi': 6}, 'S1C': {'value': 41.719, 'lli': None, 'ssi': None}, 'C1P': {'value': 23539032.74, 'lli': None, 'ssi': 6}, 'L1P': {'value': 125829714.502, 'lli': 0, 'ssi': 6}, 'D1P': {'value': -4149.698, 'lli': None, 'ssi': 6}, 'S1P': {'value': 41.062, 'lli': None, 'ssi': None}, 'C2P': {'value': 23539038.067, 'lli': None, 'ssi': 6}, 'L2P': {'value': 97867622.91, 'lli': 0, 'ssi': 6}, 'D2P': {'value': -3227.451, 'lli': None, 'ssi': 6}, 'S2P': {'value': 38.531, 'lli': None, 'ssi': None}, 'C2C': {'value': 23539037.837, 'lli': None, 'ssi': 6}, 'L2C': {'value': 97867623.908, 'lli': 0, 'ssi': 6}, 'D2C': {'value': -3227.359, 'lli': None, 'ssi': 6}, 'S2C': {'value': 38.531, 'lli': None, 'ssi': 6}}
        return the observation (dictionary) first encountered, matched by *args.
        E.g. if the above dictionary is stored in dct,
        fetch(dct, 'C1P', 'C1C', 'C2P')
        will return {'value': 23539032.74, 'lli': None, 'ssi': 6}, 'C1P'
    """
    for arg in args:
        if arg in dct:
            return dct[arg], arg
    raise RuntimeError("Error. Satellite observables do not include any of {:}".format(args))

def fetchv(dct, *args):
    """ Given a dictionary containing e.g.
        R05 : {'C1C': {'value': 23539032.631, 'lli': None, 'ssi': 6}, 'L1C': {'value': 125829717.51, 'lli': 0, 'ssi': 6}, 'D1C': {'value': -4149.772, 'lli': None, 'ssi': 6}, 'S1C': {'value': 41.719, 'lli': None, 'ssi': None}, 'C1P': {'value': 23539032.74, 'lli': None, 'ssi': 6}, 'L1P': {'value': 125829714.502, 'lli': 0, 'ssi': 6}, 'D1P': {'value': -4149.698, 'lli': None, 'ssi': 6}, 'S1P': {'value': 41.062, 'lli': None, 'ssi': None}, 'C2P': {'value': 23539038.067, 'lli': None, 'ssi': 6}, 'L2P': {'value': 97867622.91, 'lli': 0, 'ssi': 6}, 'D2P': {'value': -3227.451, 'lli': None, 'ssi': 6}, 'S2P': {'value': 38.531, 'lli': None, 'ssi': None}, 'C2C': {'value': 23539037.837, 'lli': None, 'ssi': 6}, 'L2C': {'value': 97867623.908, 'lli': 0, 'ssi': 6}, 'D2C': {'value': -3227.359, 'lli': None, 'ssi': 6}, 'S2C': {'value': 38.531, 'lli': None, 'ssi': 6}}
        return the observation (dictionary) first encountered, matched by *args.
        E.g. if the above dictionary is stored in dct,
        fetch(dct, 'C1P', 'C1C', 'C2P')
        will return 23539032.74
    """
    for arg in args:
        if arg in dct:
            return dct[arg]['value']
    raise RuntimeError("Error. Satellite observables do not include any of {:}".format(args))

class GnssRinex:
    
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
        with open(fn, 'r') as fin:
            line  = fin.readline()
            self.version = float(line[0:9])
            self.type = line[20]
            self.system = line[40]
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
                else:
                    pass
                line = fin.readline()

    def approx_cartesian(self): return [self.xapprox, self.yapprox, self.zapprox]

    class DataBlock: 

        def __init__(self, dct):
            self.dct = dct

        def t(self): return self.dct['epoch']
        def flag(self): return self.dct['flag']
        def nsats(self): return self.dct['num_sats']
        def satellite(self, satid): return self.dct[satid]
        def filter_satellite_system(self, sat_sys, include_block_info=True):
            if include_block_info:
                new_dct = {'epoch': self.dct['epoch'], 'flag': self.dct['flag'], 'num_sats': 0}
            else:
                new_dct = {}
            for k,v in self.dct.items():
                if k not in ['epoch', 'flag', 'num_sats']:
                    if sysstr2sysid(sat_sys) == sysstr2sysid(k[0]):
                        new_dct[k] = v
                        if include_block_info: new_dct['num_sats'] += 1
            return GnssRinex.DataBlock(new_dct)
        def __iter__(self): 
            return ((key, value) for key, value in self.dct.items() if key not in ['epoch', 'flag', 'num_sats'])
    
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
        line = self.stream.readline()
        if not line: raise StopIteration
        assert line[0] == '>'
        data_block = {}
        data_block['epoch'] = self.resolve_date(line[1:29])
        data_block['flag'] = int(line[29:32])
        data_block['num_sats'] = int(line[32:35])
        for i in range(data_block['num_sats']):
            line = self.stream.readline()
            satid = line[0:3]
            data_block[satid] = {}
            obs_to_follow = len(self.obscodes[satid[0]])
            for idx, obscode in enumerate(self.obscodes[satid[0]]):
                idx = 3 + idx * 16
                if line[idx:idx+16].strip() != "":
                    value = float(line[idx:idx+14])
                    try:
                        lli = int(line[idx+14]) if line[idx+14] != " " else None
                        ssi = int(line[idx+15]) if line[idx+15] != " " else None
                    except:
                        pass
                    data_block[satid][obscode] = {'value': value, 'lli': lli, 'ssi': ssi}
        return self.DataBlock(data_block)
