import datetime
from dsoclasses.sinex import sinex
from dsoclasses.geodesy import transformations
import numpy as np
import math

class ItrfSiteMotion:

    class LinearTerm:
        def __init__(self, t0, x0, vx, valid_from=datetime.datetime.min, valid_to=datetime.datetime.max):
            self.t0 = t0
            self.x0 = x0
            self.vx = vx
            self.valid_from = valid_from
            self.valid_to = valid_to
        def value(self, t):
            if t<self.valid_from or t>self.valid_to:
                raise RuntimeError
            dt = (t-self.t0).days / 365.25
            return self.x0 + dt * self.vx
    
    class PsdTerm:
        def __init__(self, amp, tau, terq, logexp):
            self.amp = amp
            self.tau = tau
            self.terq = terq
            self.logexp = logexp
        def value(self, t):
            if t >= self.terq:
                arg = (t-self.terq).days / 365.25 / self.tau
                if self.logexp == "log":
                    return self.amp * math.log10(1. + arg)
                else:
                    return self.amp * (1. - math.exp(-arg))

    def allow_extrapolation_in_linear_terms(self):
        for i in range(1,len(self.xlt)):
            assert self.xlt[i-1].valid_from < self.xlt[i].valid_from
            assert self.xlt[i-1].valid_to < self.xlt[i].valid_to
# the following should be valid, but i have seen circumstances where it isn't
# E.g. ITRF2020, site ARFB, file ITRF2020-IDS-TRF.SSC
            # assert self.xlt[i-1].valid_to <= self.xlt[i].valid_from
# so, relax the check a little bit
            assert((self.xlt[i].valid_from >= self.xlt[i-1].valid_to) or ((self.xlt[i-1].valid_to - self.xlt[i].valid_from).days < 2))
        self.xlt[0].valid_from = datetime.datetime.min
        self.xlt[-1].valid_to = datetime.datetime.max
        for i in range(1,len(self.ylt)):
            assert self.ylt[i-1].valid_from < self.ylt[i].valid_from
            assert self.ylt[i-1].valid_to < self.ylt[i].valid_to
            # assert self.ylt[i-1].valid_to <= self.ylt[i].valid_from
            assert((self.ylt[i].valid_from >= self.ylt[i-1].valid_to) or ((self.ylt[i-1].valid_to - self.ylt[i].valid_from).days < 2))
        self.ylt[0].valid_from = datetime.datetime.min
        self.ylt[-1].valid_to = datetime.datetime.max
        for i in range(1,len(self.zlt)):
            assert self.zlt[i-1].valid_from < self.zlt[i].valid_from
            assert self.zlt[i-1].valid_to < self.zlt[i].valid_to
            # assert self.zlt[i-1].valid_to <= self.zlt[i].valid_from
            assert((self.zlt[i].valid_from >= self.zlt[i-1].valid_to) or ((self.zlt[i-1].valid_to - self.zlt[i].valid_from).days < 2))
        self.zlt[0].valid_from = datetime.datetime.min
        self.zlt[-1].valid_to = datetime.datetime.max
        return
    
    def enu2xyz(self, x0, y0, z0, de, dn, du):
        lat, lon, hgt = transformations.car2ell(x0, y0, z0)
        R = transformations.geodetic2lvlh(lat, lon)
        return R @ np.array((de, dn, du))

    def at(self, t):
        x = 0e0; y = 0e0; z = 0e0;
        de = 0e0; dn = 0e0; du = 0e0;
        for xl in self.xlt:
            try: x += xl.value(t)
            except: pass
        for yl in self.ylt:
            try: y += yl.value(t)
            except: pass
        for zl in self.zlt:
            try: z += zl.value(t)
            except: pass
        if x == .0 or y == .0 or z == .0:
            print("ERROR. Failed to find suitable interval for linear term at {:}".format(t))
            raise RuntimeError
        for ep in self.epsd:
            de += ep.value(t)
        for np in self.npsd:
            dn += np.value(t)
        for up in self.upsd:
            du += up.value(t)
        dcar = self.enu2xyz(x, y, z, de, dn, du)
        return x + dcar[0], y + dcar[1], z + dcar[2]
    
    def get_psd_terms(self, estimates):
        epsd_terms = []
        npsd_terms = []
        upsd_terms = []
        for ctype in ['E', 'N', 'U']:
            site = estimates[0]['code']
            for estimate in estimates:
                if estimate['parameter_type'] == 'ALOG_' + ctype:
                    ne = self.PsdTerm(estimate['value'], 0e0, estimate['epoch'], 'log')
                elif estimate['parameter_type'] == 'TLOG_' + ctype:
                    assert ne.terq == estimate['epoch']
                    ne.tau = estimate['value']
                    if   ctype == 'E': epsd_terms.append(ne)
                    elif ctype == 'N': npsd_terms.append(ne)
                    else             : upsd_terms.append(ne)
        return epsd_terms, npsd_terms, upsd_terms

    def get_linear_terms(self, estimates):
# get (max) number of individual solutions
        num_solns = 0
        for estimate in estimates:
            if estimate['parameter_type'] in ['STAX', 'STAY', 'STAZ']:
                num_solns = max(num_solns, estimate['solution_id'])
# to be returned; a list of piece-wise linear models per component
        xterms=[]; yterms=[]; zterms=[];
# gert parameter of a givel dolution id (from list of estimates)
        def get_parameter(ptype, solnid):
            for estimate in estimates:
                if estimate['parameter_type'] == ptype and estimate['solution_id'] == solnid:
                    return estimate
            return None
# arrange pice-wise models per component
        for ctype in ['X', 'Y', 'Z']:
            for soln in range(1,num_solns+1):
                a = get_parameter('STA' + ctype, soln)
                assert a['unit'] == "m"
                b = get_parameter('VEL' + ctype, soln)
                assert b['unit'] == "m/y"
                for p in ['code', 'start', 'stop', 'point', 'epoch']:
                    assert a[p] == b[p]
                if ctype == 'X': xterms.append(self.LinearTerm(a['epoch'], a['value'],b['value'],a['start'],a['stop']))
                if ctype == 'Y': yterms.append(self.LinearTerm(a['epoch'], a['value'],b['value'],a['start'],a['stop']))
                if ctype == 'Z': zterms.append(self.LinearTerm(a['epoch'], a['value'],b['value'],a['start'],a['stop']))
        return xterms, yterms, zterms

    def __init__(self, site, *args):
# collect here all (individual) estimates, from all SINEX files
        all_estimates = []
        for fn in args:
# create a SINEX instance
            snx = sinex.Sinex(fn)
# check if site is included in the SINEX instance
            if not snx.includes_site(site):
                print("No entries for site {:} in SINEX file {:}".format(site, fn))
            else:
# collect (all) estimates in SINEX for the site
                all_estimates += snx.solution_estimate([site], True)
# given estimates, construct LinearTerm's for X-, Y- and Z- components
        self.xlt, self.ylt, self.zlt = self.get_linear_terms(all_estimates)
        self.epsd, self.npsd, self.upsd = self.get_psd_terms(all_estimates)
        self.allow_extrapolation_in_linear_terms()
