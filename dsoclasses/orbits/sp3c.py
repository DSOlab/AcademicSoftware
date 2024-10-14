#! /usr/bin/python

import sys
from datetime import datetime

def approx_gast(t):
    sec_day = (t - t.replace(hour=0, minute=0, second=0, microsecond=0)).\
                total_seconds()
    omega_earth = 7.292115e-5 ## [rad/s]
    return omega_earth * sec_day

def R3(theta):
    r = R.from_euler('z', theta, degrees=False).as_matrix()

class Sp3:

    def get_satellite(self, satid, toSI=True):
        """ Extract position and velocity (if present) info from an Sp3 
            instance.
            
            The function can only extract infor for one satellite. For more
            satellites, the function needs to be called again (one time for 
            each satellite).

            The info are returned as a dictionary, where the keys are the 
            time stampes (i.e. epochs) and the values are dictionaries of 
            info collected. E.g.

        """

        if satid not in self.sat_ids:
            print('Warning! Satellite {:} not in sp3 file'.format(satid))
            return {}

        data = {}
        with open(self.fn, 'r') as fin:
            fin.seek(self.end_of_head)
            line = fin.readline()
            while line:
# should be a header line or a line marking the end of the file (i.e. EOF)
                if line.startswith('EOF'): break
                if not line.startswith('*  '):
                    print('Error. Expected header line, found: {:}'.format(line.strip()))
                    raise RuntimeError
                t = datetime.strptime(line[3:29], "%Y %m %d %H %M %S.%f")
                epoch_entries = {}
# keep on reading lines that can be 'P', 'V'. there can also be corellation 
# lines, starting eith 'E[PV]' but these are just skipped
                line = fin.readline()
                while line[0] in ['P', 'V', 'E']:
                    if line[0] == 'P':
# parse position and clock info
                        csatid = line[1:4]
                        if csatid == satid:
                            pscale = 1e3 if toSI == True else 1e0
                            cscale = 1e6 if toSI == True else 1e0
                            x, y, z, clk = [float(x) for x in line[4:62].split()]
                            epoch_entries['x'] = x*pscale; # default is km
                            epoch_entries['y'] = y*pscale; # default is km
                            epoch_entries['z'] = z*pscale; # default is km
                            epoch_entries['c'] = clk*cscale; # default is microsec
                            if len(line) > 61:
                                xsdev, ysdev, zsdev, csdev = [int(x) for x in line[61:74].split()]
                                epoch_entries['xsdev'] = xsdev; # 
                                epoch_entries['ysdev'] = ysdev; # 
                                epoch_entries['zsdev'] = zsdev; # 
                                epoch_entries['csdev'] = csdev; # 
                            if len(line) > 74:
                                events = (line[74]+line[75]+line[78]+line[79]).strip()
                                epoch_entries['events'] = events # string

                    elif line[0] == 'V':
# parse velocity and clock info
                        csatid = line[1:4]
                        if csatid == satid:
                            vscale = 1e-1 if toSI == True else 1e0
                            cscale = 1e-6 if toSI == True else 1e0
                            vx, vy, vz, clk_rate = [float(x) for x in line[4:62].split()]
                            epoch_entries['vx'] = vx*vscale; # default is dm/sec
                            epoch_entries['vy'] = vy*vscale; # default is dm/sec
                            epoch_entries['vz'] = vz*vscale; # default is dm/sec
                            epoch_entries['cr'] = clk_rate*cscale; # 10**-4 microsec/sec
                            if len(line)>61:
                                xvsdev, yvsdev, zvsdev, crsdev = [int(x) for x in line[61:74].split()]
                                epoch_entries['xvsdev'] = xvsdev; # 
                                epoch_entries['yvsdev'] = yvsdev; # 
                                epoch_entries['zvsdev'] = zvsdev; # 
                                epoch_entries['crsdev'] = crsdev; # 
# if the line starts with a 'E' it could be EOF
                    elif line.startswith('EOF'): break
                    else:
# neither position nor velocity line, read next
                        pass
# read next line
                    line = fin.readline()
# we are through with this block; store the values collected, using the epoch
# as key
                data[t] = epoch_entries
            return data


    # Create an Sp3 instance from a filename
    def __init__(self, fn):
        self.fn = fn
        with open(fn, 'r') as fin:
# line #1
            line = fin.readline()
            self.pos_vel    =  line[2]
            self.t_start = datetime.strptime(line[3:29], "%Y %m %d %H %M %S.%f")
            self.num_epochs = int(line[32:40])
            self.crd_system = line[46:52].strip()
            self.orb_type   = line[52:56].strip()
            self.agency     = line[56:].strip()
            
# line #2 (ignore)
            line = fin.readline()
            
# lines #[3, 7], all start with '+ ' and contain satellites included in the 
# Sp3-c.
# Note that since sp3-d, the lines starting with '+ ' can be more than five.
# So here, we are going to read as many lines as are needed, i.e. all lines 
# starting with '+ ', but they should be more than five.
            self.sat_ids = []
            num_lines_sat_id = 0
            line = fin.readline()
            while line.startswith('+ '):
                num_lines_sat_id += 1
# satellite id's start at column 9 and takeup 3 characters (each)
                self.sat_ids += [sat for sat in [line[i:i+3].strip() \
                    for i in range(9, len(line), 3)] if (sat != "" and sat != "0")]
                line = fin.readline()
            if num_lines_sat_id < 5:
                print('ERROR Expected a minimum of five lines to start with \'+ \' but found {:}'.format(num_lines_sat_id))
                raise RuntimeError
                
# lines #[8, 12], all start with '++' and contain accuracies, for satellites 
# included in the sp3-c.
# Note that since sp3-d:
# a. this block can start at a later line than 8,
# b. this block can have nore than five lines, however the number of lines
#    should exactly match the num_lines_sat_ids
# The first line to examine (already buffered), should be the first line of 
# this block.
            self.sat_acc = []
            for j in range(0, num_lines_sat_id-1):
                if not line.startswith('++'):
                    print('ERROR Expected line to start with \'++\' but found {:}'.format(line.strip()))
                    raise RuntimeError
# satellite accuracies start at column 9 and takeup 3 characters (each)
                self.sat_acc += [int(acc) for acc in [line[i:i+3].strip() \
                    for i in range(9, len(line), 3)] if acc != ""]
                line = fin.readline()

# line #13 -- we only case about the time system. For sp3-d, this can be a
# latter line
            self.time_sys = fin.readline()[9:13]

# line #14 has no interesting info. For sp3-d, this can be a latter line
            line = fin.readline()
            
# line #15 has base values for accuracies. For sp3-d, this can be a latter line
            line = fin.readline()
            self.base_posvel = float(line[3:14])
            self.base_clkrate= float(line[14:27])

# lines #[16, 22] are skiped, no info.
            for j in range(16, 23): line = fin.readline()
# if we are reading an sp3d file, there may be more comment lines. Keep reading 
# until we mee the first line that does not start with '/*'
# Store the current position whithin the file, so we can later skip the 
# header easily
            self.end_of_head = fin.tell()
            while line.startswith('/*'): 
                self.end_of_head = fin.tell()
                line = fin.readline()

