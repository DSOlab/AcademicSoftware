import datetime

def parse_sinex_date(dstr):
    """ Parse a ate string geiven in the format: YYYY:DDD:SSSSS, where :
        * YYYY is the year,
        * DDD is the day of year, and
        * SSSSS is the seconds of day (integral)
    """
    l = dstr.strip().split(":")
    if len(l) != 3:
        print("Error. Failed resolving SINEX date: {:}".format(dstr.strip()), file=sys.stderr)
        raise RuntimeError
    date = datetime.datetime.strptime(" ".join([l[0],l[1]]), "%Y %j")
    time = datetime.timedelta(seconds=int(l[2]))
    return date + time

class SinexBias:
    
    def read_header(self, fn):
        with open(fn, 'r') as fin:
            line = fin.readline()
            assert line[0:5] == "%=BIA"
            self.vesion = float(line[6:6+4])
            self.start = parse_sinex_date(line[34:34+14])
            self.stop  = parse_sinex_date(line[49:49+14])
            self.bias_mode = line[64]
            self.num_estimates = int(line[66:])

    def __init__(self, fn):
        self.read_header(fn)
        self.fn = fn

    def prn_matches(self, prn, prn_list):
        if prn_list == []: return True
        if prn[0] in prn_list: return True
        return True if prn in prn_list else False
    def station_matches(self, station, station_list):
        if station_list is None: return (station == '')
        return True if (station_list == [] or station.upper() in station_list) else False
    def obs_matches(self, obs, obs_list):
        if obs_list == []: return True
        # ablosute obervable (OBS2 is missing, should be None)
        if obs[1] == '': obs = obs[0]
        return True if obs in obs_list else False

    def get(self, bias_types, prn_list=[], station_list=[], obs_list=[]):
        """ An empty (input) list denotes everything. E.g., if prn_list=[], 
            then all PRN's will be considered a match, even if they only 
            designate a GNSS (e.g. 'G', or 'E').
            
            If we don;t want to consider stations (i.e. only consider lines 
            with the field 'STATION' missing), then station_list should be 
            set to None.
        """
        out = []
        with open(self.fn, 'r') as fin:
            line = fin.readline()
            while line and line.strip() != '+BIAS/SOLUTION': line = fin.readline()
            while line and line.strip() != '-BIAS/SOLUTION':
                line = fin.readline()
#*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___
#      6    11  15        25   30   35             50             65   70                    92
                bias = line[0:6].strip()
                if bias in bias_types:
                    prn = line[11:15].strip()
                    if self.prn_matches(prn, prn_list):
                        station = line[15:25].strip()
                        if self.station_matches(station, station_list):
                            obs = (line[25:30].strip(), line[30:35].strip())

                            if self.obs_matches(obs, obs_list):
                                bias_start = parse_sinex_date(line[35:50].strip())
                                bias_end   = parse_sinex_date(line[50:65].strip())
                                assert line[65:70].strip() == 'ns'
                                bias_value = float(line[70:92].strip())
                                bias_std   = float(line[92:].strip())
                                out.append({'bias': bias, 'prn': prn, 'station': station, 'obs': obs, 'start': bias_start, 'end': bias_end, 'value': bias_value, 'std': bias_std})
        return out
