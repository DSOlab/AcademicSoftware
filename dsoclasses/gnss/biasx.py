import datetime

def _bias_sinex_parse_date(dstr):
    year_str, doy_str, sod_str = dstr.split(':')
    year = int(year_str)
    doy  = int(doy_str)
    sod  = int(sod_str)  # seconds of day
    dt = datetime.datetime(year, 1, 1) + datetime.timedelta(days=doy-1, seconds=sod)
    return dt

class Bsinex:
    
    def __init__(self, fn):
        self._fn = fn
        self._sat_clock_ref_obs = {}
        self._sat_bias = {}

    def parse_description_block(self):
        with open(self._fn, 'r', encoding="utf-8", errors="ignore") as fin:
            line = fin.readline()
            while line and line.strip() != '+BIAS/DESCRIPTION':
                line = fin.readline()

            while line and line.strip() != '-BIAS/DESCRIPTION':
                if line.lstrip().startswith('BIAS_MODEL'):
                    self._bias_model = line[41:].strip()
                elif line.lstrip().startswith('APC_MODEL'):
                    self._apc_model = line[41:].strip()
                elif line.lstrip().startswith('TIME_SYSTEM'):
                    self._time_system = line[41:].strip()
                elif line.lstrip().startswith('SATELLITE_CLOCK_REFERENCE_OBSERVABLES'):
                    info = line[41:].split()
                    # at least the satellite system should follow, e.g.
                    # G C1W C2W, or
                    # R
                    assert(len(info)>0)
                    system = info[0]
                    assert system in ['GREC']
                    refs = info[1:] if len(info)>0 else []
                    self._sat_clock_ref_obs[system] = refs
                else:
                    pass
                line  = fin.readline()
        return

    def sat_bias(self, prn):
        if self._sat_bias == {}: self._sat_bias = self._parse_sat_bias()
        return self._sat_bias[prn]

    def _parse_sat_bias(self):
        """
        +BIAS/SOLUTION
        *BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___ __ESTIMATED_SLOPE____ _STD_DEV___
         OSB  G063 G01           C1C       2024:001:00000 2024:002:00000 ns                  9.5274      0.0168
         OSB  E209 E09           C5X       2024:001:00000 2024:002:00000 ns                  2.9844      0.0158
         OSB  E224 E10           C1C       2024:001:00000 2024:002:00000 ns                  0.1316      0.0000
        """
        if self._sat_bias != {}: return self._sat_bias;
        _sat_bias = {}
        with open(self._fn, 'r', encoding="utf-8", errors="ignore") as fin:
            line = fin.readline()
            while line and line.strip() != '+BIAS/SOLUTION':
                line = fin.readline()
            
            while line and line.strip() != '-BIAS/SOLUTION':
                fields = line.split()
                if (not line.startswith('*')) and line[15:24].strip()=='' and fields[0] == 'OSB' and fields[1][0] in 'GREC' and len(fields[1])==4 and fields[2][0] in 'GREC' and len(fields[2])==3:
                    prn  = fields[2]
                    obs1  = fields[3]
                    start = _bias_sinex_parse_date(fields[4])
                    stop  = _bias_sinex_parse_date(fields[5])
                    units = fields[6]
                    val   = float(fields[7])
                    if prn not in _sat_bias:
                        _sat_bias[prn] = {'from': start, 'to': stop, 'units': units, obs1: val}
                    else:
                        assert units == _sat_bias[prn]['units']
                        assert (obs1 not in _sat_bias[prn])
                        assert (start == _sat_bias[prn]['from'] and stop == _sat_bias[prn]['to'])
                        _sat_bias[prn][obs1] = val
                line = fin.readline()
        return _sat_bias
