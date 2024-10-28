import datetime
import math

class GnssRinex:
    
    def resolve_date(self, dstr):
        l = dstr.split()
# maximum accuracy in datetime instances is microseconds
        assert l[-1][-1] == "0"
        y, m, d, h, mn = [int(x) for x in l[0:5]]
        sec = float(l[5])
        s = int(sec)
        fmicrosec = (sec - s) * 1e6
        tstr = "{:4d}:{:02d}:{:02d} {:02d}:{:02d}:{:02d}".format(y, m, d, h, mn, s)
        nofrag_t = datetime.datetime.strptime(tstr, "%Y:%m:%d %H:%M:%S")
        return nofrag_t + datetime.timedelta(microseconds=int(fmicrosec))

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
                    l = line[7:60].split()
                    while len(l) < numobsc:
                        line = fin.readline()
                        l += line[7:60].split()
                    self.obscodes[system] = l
                elif line[60:].strip() == "TIME OF FIRST OBS":
                    self.time_first_obs = self.resolve_date(line[0:44])
                    self.time_first_obs_sys = line[44:60].strip()
                elif line[60:].strip() == "TIME OF LAST OBS":
                    self.time_last_obs = self.resolve_date(line[0:44])
                    self.time_last_obs_sys = line[44:60].strip()
                else:
                    pass
                line = fin.readline()

    class DataBlock: 

        def __init__(self, dct):
            self.dct = dct

        def t(self): return self.dct['epoch']
        def flag(self): return self.dct['flag']
        def nsats(self): return self.dct['num_sats']
        def satellite(self, satid): return self.dct[satid]
    
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
        data_block['num_sats'] = int(line[32:])
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
