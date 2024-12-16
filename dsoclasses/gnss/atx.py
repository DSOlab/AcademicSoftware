import datetime
import attotime
import math
import sys
import numpy as np

class ReceiverAntennaPcv:
    def __init__(self, antenna):
        self.antenna = antenna
        self.pcv = {}
    def add_freq(self, freq, neu, z1z2dz, vals):
        d = {'neu': neu, 'z1z2dz': z1z2dz, 'pcv': vals}
        self.pcv[freq] = d
    def pco(self, freq): return self.pcv[freq]['neu']
    def pcv_hgt(self, freq, el): 
        za = np.pi / 2. - el
        assert za >= 0e0 and za <= np.pi / 2.
        cp = int((np.degrees(za) - self.pcv[freq]['z1z2dz'][0]) / self.pcv[freq]['z1z2dz'][2])
        np1 = cp + 1
        x0 = self.pcv[freq]['z1z2dz'][0] + cp * self.pcv[freq]['z1z2dz'][2]
        x1 = x0 + self.pcv[freq]['z1z2dz'][2]
        x  = np.degrees(za)
        assert np.degrees(za) >= x0 and np.degrees(za) < x1
        y0 = self.pcv[freq]['pcv'][cp]
        y1 = self.pcv[freq]['pcv'][np1]
        return (y0*(x1-x) + y1*(x-x0)) / (x1-x0)

class Atx:

    def parse_header(self, fn):
        with open(fn, 'r') as fin:
            line  = fin.readline()
            self.version = float(line[0:8])
            self.sat_sys = line[20]
            assert self.sat_sys in ['G', 'R', 'E', 'C', 'J', 'S', 'M']
            assert line[60:].rstrip() == 'ANTEX VERSION / SYST'
            while line and line.strip() != "END OF HEADER":
                if line[60:].strip() == "PCV TYPE / REFANT":
                    self.variation_type = line[0]
                    assert self.variation_type in ['A', 'R']
                    if self.variation_type == 'R':
                        print('ATX {:} holds relative PCVs; cannot handle this type of information!'.format(fn), file=sys.stderr)
                line = fin.readline()
            self.filename = fn


    def __init__(self, fn): self.parse_header(fn)

    def goto_antenna(self, antenna):
        fin = open(self.filename, 'r')
        antenna_found = False
        line = fin.readline()
        while not antenna_found and line:
            while line and line.strip() != "START OF ANTENNA":
                line = fin.readline()
# found new antenna block;
            line = fin.readline()
            if line[0:20].strip() == antenna:
                antenna_found = True
                break
            else:
                line = fin.readline()
        if antenna_found: return fin
        fin.close()
        return None

    def get_noazi(self, antenna, freq_list):
        if len(antenna) < 20: antenna = "{:16s}NONE".format(antenna)
        freqs_collected = []
        fin = self.goto_antenna(antenna)
        if fin is None:
            print("ERROR. Failed locating antenna \'{:}\' in atx file {:}".format(antenna, self.filename), file=sys.stderr)
            return None

        pcv = ReceiverAntennaPcv(antenna)
        line = fin.readline()
        while sorted(freqs_collected) != sorted(freq_list):
            while line and line[60:].rstrip() != "START OF FREQUENCY":
                line = fin.readline()
                if line[60:].rstrip() == "ZEN1 / ZEN2 / DZEN":
                    z1, z2, dz = [ float(x) for x in line[0:20].split() ]
                elif line[60:].rstrip() == "END OF ANTENNA":
                    print("ERROR. Failed locating requested frequencies for antenna \'{:}\'. Atx file is {:}".format(antenna, self.filename), file=sys.stderr)
                    fin.close()
                    return None
# should be at the start of a new frequency
            freq = line[3:6]
            if freq in freq_list:
                line = fin.readline()
                assert line[60:].rstrip() == "NORTH / EAST / UP"
                n, e, u = [ float(line[i*10:i*10+10])*1e-3 for i in range(3) ]
                line = fin.readline()
                assert line.startswith("   NOAZI")
                ln = line[8:]
                vals = [ float(ln[i*8:i*8+8])*1e-3 for i in range(int((z2-z1)/dz)+1) ]
                freqs_collected.append(freq)
                pcv.add_freq(freq, [n,e,u], [z1,z2,dz], vals)
            line = fin.readline()
        fin.close()
        return pcv
