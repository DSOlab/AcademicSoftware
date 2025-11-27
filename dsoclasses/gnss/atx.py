import datetime
import attotime
import math
import sys
import numpy as np


# System Frequencies as in
# https://files.igs.org/pub/data/format/antex14.txt
"""
 GPS:     'G01' - L1             
          'G02' - L2             
          'G05' - L5             
 GLONASS: 'R01' - G1             
          'R02' - G2             
 Galileo: 'E01' - E1             
          'E05' - E5a            
          'E07' - E5b            
          'E08' - E5 (E5a+E5b)   
          'E06' - E6             
 Compass: 'C01' - E1             
          'C02' - E2             
          'C07' - E5b            
          'C06' - E6             
 QZSS:    'J01' - L1             
          'J02' - L2             
          'J05' - L5             
          'J06' - LEX            
 SBAS:    'S01' - L1             
          'S05' - L5
"""
_ATX_FREQS = {
    'G': {'L1': 'G01', 'L2': 'G02', 'L5': 'G04'},
    'R': {'G1': 'R01', 'G2': 'R02', 'L1': 'R01', 'L2': 'R02'},
    'E': {'E1': 'E01', 'E5a': 'E05', 'E5b': 'E07', 'E5(E5a+E5b)': 'E08', 'E6': 'E06'},
    'C': {'E1': 'C01', 'E2': 'C02', 'E5b': 'C07', 'E6': 'C06'},
    'J': {'L1': 'J01', 'L2': 'J02', 'L5': 'J05', 'LEX': 'J06'},
    'S': {'L1': 'S01', 'L5': 'S05'}
}
def _to_atx_freq(sys, freq):
    fstr = None
    try:
        fstr = _ATX_FREQS[sys.upper()[0]][freq]
    except:
        pass
    if not fstr:
        raise RuntimeError(f'Failed to match given SYS/FREQ={sys}{freq} to a valid atx frequency')
    return fstr

class ReceiverAntennaPcv:
    def __init__(self, antenna):
        self.antenna = antenna
        self.pcv = {}
    def add_freq(self, freq, neu, z1z2dz, vals):
        d = {'neu': neu, 'z1z2dz': z1z2dz, 'pcv': vals}
        self.pcv[freq] = d
    
    # def pco(self, freq): return self.pcv[freq]['neu']
    # def pco(self, sys, freq): return self.pco(_to_atx_freq(sys, freq))
    def pco(self, *args):
        if len(args) == 1:
            (freq,) = args
            return self.pcv[freq]['neu']
        elif len(args) == 2:
            sys, freq = args
            return self.pcv[_to_atx_freq(sys, freq)]['neu']
        else:
            raise RuntimeError("ReceiverAntennaPcv.pco expects 1 or 2 arguments")

    def _pcv_hgt(self, freq, el): 
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
    def pcv_hgt(self, *args):
        if len(args) == 2:
            freq, el = args
            return self._pcv_hgt(freq, el)
        elif len(args) == 3:
            sys, freq, el = args
            return self._pcv_hgt(_to_atx_freq(sys, freq), el)
        else:
            raise RuntimeError("ReceiverAntennaPcv.pcv_hgt expects 2 or 3 arguments")

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

    def _get_noazi(self, antenna, freq_list):
        if len(antenna) < 20: antenna = "{:16s}NONE".format(antenna)
        freqs_collected = []
        fin = self.goto_antenna(antenna)
        if fin is None:
            raise RuntimeError("Failed locating antenna \'{antenna}\' in atx file {self.filename}")

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
    
    def get_noazi(self, antenna, freq_list):
        if isinstance(freq_list, list):
            return  self._get_noazi(antenna, freq_list)
        elif isinstance(freq_list, dict):
            # example: {"E": ['FREQ1', 'FREQ2'], 'G': [...], ...}
            freqs = []
            for sys, fs in freq_list.items():
                if isinstance(fs, list):
                    for f in fs:
                        freqs.append(_to_atx_freq(sys, f))
                else:
                    freqs.append(_to_atx_freq(sys, fs))
            return self._get_noazi(antenna, freqs)
        else:
            self._get_noazi(antenna, [freq_list])

