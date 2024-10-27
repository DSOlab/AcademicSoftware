import datetime
import sys
import math

def parse_sinex_date(dstr):
    """ Parse a ate string geiven in the format: YY:DDD:SSSSS, where :
        * YY is the two-digit year,
        * DDD is the day of year, and
        * SSSSS is the seconds of day (integral)
    """
    l = dstr.strip().split(":")
    if len(l) != 3:
        print("Error. Failed resolving SINEX date: {:}".format(dstr.strip()), file=sys.stderr)
        raise RuntimeError
    date = datetime.datetime.strptime(" ".join([l[0],l[1]]), "%y %j")
    time = datetime.timedelta(seconds=int(l[2]))
    return date + time

class Sinex:
    def goto_block(self, fin, block_str):
        line = fin.readline()
        while line and line != block_str and line.strip() != "%ENDSNX":
            if line.strip() == "+" + block_str.strip():
                return
            line = fin.readline()
        print("Warning. Failed finding block: {:} in SINEX".format(block_str), file=sys.stderr)
        raise RuntimeError

    def parse_header_line(self, fn):
        with open(fn, 'r') as fin:
            line = fin.readline()
            if not line.startswith('%=SNX'):
                print("Error. Failed to identify SINEX first line ({:})".format(fn), file=sys.stderr)
                raise RuntimeError
            self.version = float(line[6:10])
            self.agency  = line[11:14].strip()
            self.time_creation = parse_sinex_date(line[14:28])
            self.data_provider = line[28:31].strip() 
            self.data_start = parse_sinex_date(line[31:44])
            self.data_stop  = parse_sinex_date(line[44:57])
            self.technique  = line[58]
            self.num_estimates     = int(line[60:65])
            self.constraint_code   = int(line[66])
            self.solution_contents = line[67:].strip()
        return

    def __init__(self, fn):
        self.parse_header_line(fn)
        self.filename = fn

    def site_id(self, site_list_in=[]):
        site_list_out = []
        with open(self.filename, 'r') as fin:
            self.goto_block(fin, "SITE/ID")
            line = fin.readline()
            while line.strip() != "-SITE/ID":
                dct = {}
                if line[0] != "*":
                    dct["code"]       = line[0:5].strip()
                    dct["point"]      = line[6:8].strip()
                    dct["domes"]      = line[9:19].strip()
                    dct["technique"]  = line[19]
                    dct["description"]= line[21:44].strip()
                    l = line[44:].split()
                    dct["lon"] = math.radians(float(l[0]) + float(l[1])/ 60e0 + float(l[2])/ 60e0 / 60e0)
                    dct["lat"] = math.radians(float(l[3]) + float(l[4])/ 60e0 + float(l[5])/ 60e0 / 60e0)
                    dct["hgt"] = float(l[6])
                    if site_list_in == [] or dct["code"] in site_list_in:
                        site_list_out.append(dct)
                line = fin.readline()
            return site_list_out
            
    def includes_site(self, site): 
        sites = self.site_id([site])
        if len(sites) > 1:
            print("ERROR. Site with code {:} is included more than once in SINEX {:}".format(site, self.filename), file=sys.stderr)
        return bool(len(sites)==1)
    
    def solution_epochs(self, site_list_in=[]):
        site_list_out = []
        with open(self.filename, 'r') as fin:
            self.goto_block(fin, "SOLUTION/EPOCHS")
            line = fin.readline()
            while line.strip() != "-SOLUTION/EPOCHS":
                dct = {}
                if line[0] != "*":
                    dct["code"]       = line[0:5].strip()
                    dct["point"]      = line[6:8].strip()
                    dct["solution_id"]= int(line[9:13].strip())
                    dct["technique"]  = line[14]
                    l = line[16:].split()
                    dct["start"] = parse_sinex_date(l[0])
                    dct["stop"] = parse_sinex_date(l[1])
                    dct["meant"] = parse_sinex_date(l[2])
                    if site_list_in == [] or dct["code"] in site_list_in:
                        site_list_out.append(dct)
                line = fin.readline()
            return site_list_out
    
    def solution_estimate(self, site_list_in=[], append_solution_epoch_info=False):
        site_list_out = []
        with open(self.filename, 'r') as fin:
            self.goto_block(fin, "SOLUTION/ESTIMATE")
            line = fin.readline()
            while line.strip() != "-SOLUTION/ESTIMATE":
                dct = {}
                if line[0] != "*":
                    dct["index"]  = int(line[0:6].strip())
                    dct["parameter_type"] = line[7:13].strip()
                    dct["code"]   = line[14:18].strip()
                    dct["point"]  = line[19:21].strip()
# *psd.snx files, contain no SOLN record, and the field is instead filled 
# with '-----'
                    if line[22:26] == "----":
                        dct["solution_id"] = None
                    else:
                        dct["solution_id"] = int(line[22:26])
                    dct["epoch"] = parse_sinex_date(line[27:39].strip())
                    dct["unit"] = line[40:44].strip()
                    l = line[44:].split()
                    dct["constraint_code"] = l[0]
                    dct["value"] = float(l[1])
                    dct["std_dev"] = float(l[2])
                    if site_list_in == [] or dct["code"] in site_list_in:
                        site_list_out.append(dct)
                line = fin.readline()
            if not append_solution_epoch_info: return site_list_out
# read the SOLUTION/EPOCHS block and append solution specific info. note that 
# some SINEX files do not include a SOLUTION/EPOCHS block, e.g. ITRF SINEX 
# files for PSD
        sta_list_out_extra = []
        try:
            sol_info = self.solution_epochs(site_list_in)
        except:
            print("Warning. No block SOLUTION/EPOCHS found in SINEX file {:}".format(self.filename))
            return site_list_out

        for estimate in site_list_out:
            dct = estimate
            estimate_matched = False
            for info in sol_info:
                if info["code"] == dct["code"] and info["point"] == dct["point"] and info["solution_id"] == dct["solution_id"]:
                    dct["technique"] = info["technique"]
                    dct["start"] = info["start"] 
                    dct["stop"]  = info["stop"]
                    estimate_matched = True
                    break
            if not estimate_matched:
                print("ERROR. Failed matching SOLUTION/EPOCHS for SOLUTION/ESTIMATE!")
                #print("ERROR. SINEX {:} estimate {:}".format(self.filename, estimate))
                raise RuntimeError
            sta_list_out_extra.append(dct)
        return sta_list_out_extra
