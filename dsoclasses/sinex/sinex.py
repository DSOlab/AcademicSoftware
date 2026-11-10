import datetime
import sys
import math
import numpy as np
from dsoclasses.geodesy.transformations import geodetic2lvlh, car2ell
from dsoclasses.sinex.dpod import dpod_freq_corr
from dsoclasses.time.calmjd import cal2fmjd
from dsoclasses.time.pyattotime import to_attoseconds, at2pt


def parse_sinex_date(text):
    """Read SINEX ``YY:DDD:SSSSS`` into a naive calendar time.

    This teaching reader expects finite, explicitly dated solution epochs.
    Unknown/open dates such as ``00:000:00000`` are not silently turned into
    a real date; handling open intervals belongs to a more general parser.
    """
    try:
        year_s, doy_s, second_s = text.strip().split(':')
        if len(year_s) != 2:
            raise ValueError('Expected 2-digit SINEX year')
        yy, doy, seconds = int(year_s), int(doy_s), int(second_s)
        year = 2000 + yy if yy < 80 else 1900 + yy
        days = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
        if not (1 <= doy <= days) or not (0 <= seconds <= 86400):
            raise ValueError('Invalid SINEX day-of-year or seconds-of-day')
        return datetime.datetime(year, 1, 1) + datetime.timedelta(days=doy-1, seconds=seconds)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'Invalid finite SINEX epoch {text!r}') from exc


def _parse_interval_date(value, *, start):
    """Return a finite epoch, or an open SINEX interval endpoint.

    Unknown boundaries (00:000:00000) are allowed only for a solution's
    applicability, never for its reference-coordinate epoch.
    """
    if value.strip() == '00:000:00000':
        return datetime.datetime.min if start else datetime.datetime.max
    return parse_sinex_date(value)


class Sinex:
    def goto_block(self, fin, block_str):
        line = fin.readline()
        while line and line != block_str and line.strip() != "%ENDSNX":
            if line.strip() == "+" + block_str.strip():
                return
            line = fin.readline()
        print(
            "Warning. Failed finding block: {:} in SINEX".format(block_str),
            file=sys.stderr,
        )
        raise RuntimeError

    def parse_header_line(self, fn):
        with open(fn, "r") as fin:
            line = fin.readline()
            if not line.startswith("%=SNX"):
                print(
                    "Error. Failed to identify SINEX first line ({:})".format(fn),
                    file=sys.stderr,
                )
                raise RuntimeError
            self.version = float(line[6:10])
            self.agency = line[11:14].strip()
            self.time_creation = parse_sinex_date(line[14:28])
            self.data_provider = line[28:31].strip()
            self.data_start = _parse_interval_date(line[31:44], start=True)
            self.data_stop = _parse_interval_date(line[44:57], start=False)
            self.technique = line[58]
            self.num_estimates = int(line[60:65])
            self.constraint_code = int(line[66])
            self.solution_contents = line[67:].strip()
        return

    def __init__(self, fn):
        self.parse_header_line(fn)
        self.filename = fn

    def site_id(self, site_list_in=[]):
        site_list_out = []
        with open(self.filename, "r") as fin:
            self.goto_block(fin, "SITE/ID")
            line = fin.readline()
            while line.strip() != "-SITE/ID":
                if not line or line.strip() == "%ENDSNX":
                    raise ValueError("Truncated SINEX block SITE/ID: missing -SITE/ID")
                dct = {}
                if line[0] != "*":
                    dct["code"] = line[0:5].strip()
                    dct["point"] = line[6:8].strip()
                    dct["domes"] = line[9:19].strip()
                    dct["technique"] = line[19]
                    dct["description"] = line[21:44].strip()
                    l = line[44:].split()
                    # A negative degree value makes *all* DMS components
                    # negative; do not accidentally add positive arcminutes.
                    def dms(deg, minute, second):
                        degrees = float(deg)
                        angle = abs(degrees) + abs(float(minute))/60 + abs(float(second))/3600
                        return math.radians(math.copysign(angle, degrees))
                    dct["lon"] = dms(*l[:3])
                    dct["lat"] = dms(*l[3:6])
                    dct["hgt"] = float(l[6])
                    if site_list_in == [] or dct["code"] in site_list_in:
                        site_list_out.append(dct)
                line = fin.readline()
            return site_list_out

    def includes_site(self, site):
        sites = self.site_id([site])
        if len(sites) > 1:
            print(
                "ERROR. Site with code {:} is included more than once in SINEX {:}".format(
                    site, self.filename
                ),
                file=sys.stderr,
            )
        return bool(len(sites) == 1)

    def solution_epochs(self, site_list_in=[]):
        site_list_out = []
        with open(self.filename, "r") as fin:
            self.goto_block(fin, "SOLUTION/EPOCHS")
            line = fin.readline()
            while line.strip() != "-SOLUTION/EPOCHS":
                if not line or line.strip() == "%ENDSNX":
                    raise ValueError("Truncated SINEX block SOLUTION/EPOCHS: missing -SOLUTION/EPOCHS")
                dct = {}
                if line[0] != "*":
                    dct["code"] = line[0:5].strip()
                    dct["point"] = line[6:8].strip()
                    dct["solution_id"] = int(line[9:13].strip())
                    dct["technique"] = line[14]
                    l = line[16:].split()
                    dct["start"] = _parse_interval_date(l[0], start=True)
                    dct["stop"] = _parse_interval_date(l[1], start=False)
                    dct["meant"] = parse_sinex_date(l[2])
                    if site_list_in == [] or dct["code"] in site_list_in:
                        site_list_out.append(dct)
                line = fin.readline()
        return site_list_out

    def eccentricity(self, site_list_in=[]):
        site_list_out = []
        with open(self.filename, "r") as fin:
            self.goto_block(fin, "SITE/ECCENTRICITY")
            line = fin.readline()
            while line.strip() != "-SITE/ECCENTRICITY":
                if not line or line.strip() == "%ENDSNX":
                    raise ValueError("Truncated SINEX block SITE/ECCENTRICITY: missing -SITE/ECCENTRICITY")
                dct = {}
                if line[0] != "*":
                    dct["code"] = line[0:5].strip()
                    dct["point"] = line[6:8].strip()
                    dct["solution_id"] = int(line[9:13].strip())
                    dct["technique"] = line[14]
                    l = line[16:].split()
                    dct["start"] = _parse_interval_date(l[0], start=True)
                    dct["stop"] = _parse_interval_date(l[1], start=False)
                    dct["axe"] = l[2]
                    dct["up"], dct["north"], dct["east"] = [float(d) for d in l[3:]]
                    if site_list_in == [] or dct["code"] in site_list_in:
                        site_list_out.append(dct)
                line = fin.readline()
        return site_list_out

    def solution_estimate(self, site_list_in=[], append_solution_epoch_info=False):
        site_list_out = []
        with open(self.filename, "r") as fin:
            self.goto_block(fin, "SOLUTION/ESTIMATE")
            line = fin.readline()
            while line.strip() != "-SOLUTION/ESTIMATE":
                if not line or line.strip() == "%ENDSNX":
                    raise ValueError("Truncated SINEX block SOLUTION/ESTIMATE: missing -SOLUTION/ESTIMATE")
                dct = {}
                if line[0] != "*":
                    dct["index"] = int(line[0:6].strip())
                    dct["parameter_type"] = line[7:13].strip()
                    dct["code"] = line[14:18].strip()
                    dct["point"] = line[19:21].strip()
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
            if not append_solution_epoch_info:
                return site_list_out
        # read the SOLUTION/EPOCHS block and append solution specific info. note that
        # some SINEX files do not include a SOLUTION/EPOCHS block, e.g. ITRF SINEX
        # files for PSD
        sta_list_out_extra = []
        try:
            sol_info = self.solution_epochs(site_list_in)
        except RuntimeError as exc:
            # A static single-epoch/PSD SINEX may omit this block, but
            # coordinate extrapolation requires its validity intervals.
            raise LookupError('SINEX lacks SOLUTION/EPOCHS; cannot select a '
                              'valid station coordinate solution') from exc

        for estimate in site_list_out:
            dct = estimate
            estimate_matched = False
            for info in sol_info:
                if (
                    info["code"] == dct["code"]
                    and info["point"] == dct["point"]
                    and info["solution_id"] == dct["solution_id"]
                ):
                    dct["technique"] = info["technique"]
                    dct["start"] = info["start"]
                    dct["stop"] = info["stop"]
                    estimate_matched = True
                    break
            if not estimate_matched:
                raise LookupError(
                    f'SINEX SOLUTION/ESTIMATE has no matching SOLUTION/EPOCHS '
                    f'for station {dct["code"]}, solution {dct["solution_id"]}'
                )
            sta_list_out_extra.append(dct)
        return sta_list_out_extra


def extract_sinex_coordinates(sinex_fn, site_list, epoch, apply_eccentricity=False, dpod_freq_corr_fn=None):
    """Propagate SINEX station coordinates to ``epoch`` using station velocity.

    Arguments
    ---------
    sinex_fn: file path to a station coordinate/velocity SINEX product.
    site_list: iterable of four-character station identifiers.
    epoch: datetime, attodatetime or other supported precise epoch label.
    apply_eccentricity: if true, add station SITE/ECCENTRICITY to XYZ.
    dpod_freq_corr_fn: optional DORIS periodic-displacement table.

    Returns
    -------
    ``{site: {'X': m, 'Y': m, 'Z': m, 'soln': solution_id}}``.

    Important conventions
    ---------------------
    Positions refer to the stated reference frame and reference epoch of the
    SINEX file. Velocities are metres per Julian year (365.25 x 86400 s).
    This method does not transform between terrestrial reference frames or
    account for nonlinear motion, earthquakes, or antenna phase centres.

    All XYZ and VELXYZ values must belong to the *same point/solution* and
    have a complete validity interval. If this cannot be verified, fail with
    an informative exception rather than return partial station coordinates.
    """
    if isinstance(site_list, str):
        site_list = [site_list]
    site_list = list(site_list)
    if not site_list:
        raise ValueError('site_list must contain at least one station code')
    instant = to_attoseconds(epoch)
    display_epoch = at2pt(epoch)

    snx = Sinex(sinex_fn)
    # SOLUTION/EPOCHS supplies the validity and solution number. Estimate
    # rows alone do not identify the complete time validity unambiguously.
    entries = snx.solution_estimate(site_list, append_solution_epoch_info=True)
    result = {}
    for site in site_list:
        valid = []
        for entry in entries:
            if entry['code'] != site:
                continue
            if 'start' not in entry or 'stop' not in entry:
                raise LookupError(f'SINEX {site}: missing SOLUTION/EPOCHS interval')
            if to_attoseconds(entry['start']) <= instant <= to_attoseconds(entry['stop']):
                valid.append(entry)
        if not valid:
            raise LookupError(f'No valid SINEX coordinate solution for {site} at {display_epoch}')

        # Group by monument point and solution number; never accidentally
        # combine the XYZ of one solution with velocities of another.
        solutions = {}
        for item in valid:
            key = (item['point'], item['solution_id'])
            parameters = solutions.setdefault(key, {})
            name = item['parameter_type']
            if name in parameters:
                raise ValueError(f'Duplicate SINEX parameter {name} for {site}, solution {key}')
            parameters[name] = item
        if len(solutions) != 1:
            raise LookupError(f'SINEX {site}: ambiguous overlapping point/solution records')
        required = {'STA'+axis for axis in 'XYZ'} | {'VEL'+axis for axis in 'XYZ'}
        matches = [(key, p) for key, p in solutions.items() if required <= p.keys()]
        if len(matches) != 1:
            raise LookupError(f'Expected one complete XYZ/velocity solution for {site} '
                              f'at {display_epoch}, found {len(matches)}')
        (point, soln), parameters = matches[0]
        ref_epochs = {to_attoseconds(parameters['STA'+axis]['epoch']) for axis in 'XYZ'}
        if len(ref_epochs) != 1:
            raise ValueError(f'SINEX {site} XYZ components have different reference epochs')
        years_since_ref = (instant - ref_epochs.pop()) / (365.25 * 86400 * 10**18)
        coords = {}
        for axis in 'XYZ':
            position = parameters['STA'+axis]
            velocity = parameters['VEL'+axis]
            if position['unit'] != 'm' or velocity['unit'] != 'm/y':
                raise ValueError(f'Unsupported SINEX coordinate/velocity units for {site} {axis}')
            xyz = float(position['value']) + float(velocity['value']) * years_since_ref
            if not np.isfinite(xyz):
                raise ValueError(f'Non-finite propagated coordinate for {site} {axis}')
            coords[axis] = xyz
        coords['soln'] = soln
        result[site] = coords

    if apply_eccentricity:
        eccentricities = snx.eccentricity(site_list)
        for site, crd in result.items():
            matches = [r for r in eccentricities if r['code'] == site
                       and r['solution_id'] == crd['soln']
                       and to_attoseconds(r['start']) <= instant <= to_attoseconds(r['stop'])]
            if len(matches) != 1:
                raise LookupError(f'Expected one valid SITE/ECCENTRICITY for {site}, '
                                  f'found {len(matches)}')
            entry = matches[0]
            if entry['axe'] != 'UNE':
                raise ValueError(f'Unsupported SINEX eccentricity axes {entry["axe"]!r} for {site}; need UNE')
            # The SINEX eccentricity is specified as U,N,E (metres). The
            # geodetic2lvlh matrix has E,N,U columns, so reorder explicitly.
            local_enu = np.array([entry['east'], entry['north'], entry['up']])
            lat, lon, _ = car2ell(crd['X'], crd['Y'], crd['Z'])
            vec = np.array([crd['X'], crd['Y'], crd['Z']]) + geodetic2lvlh(lat, lon) @ local_enu
            crd.update(zip('XYZ', vec.tolist()))

    if dpod_freq_corr_fn:
        corrections = dpod_freq_corr(dpod_freq_corr_fn, site_list)
        # This optional DORIS deformation model is not used in the GNSS
        # tutorial; retain its old interface with corrected fractional-DOY.
        from dsoclasses.time.calmjd import cal2fmjd
        fmjd = cal2fmjd(display_epoch)
        fdoy = display_epoch.timetuple().tm_yday + (fmjd - math.floor(fmjd))
        for site, crd in result.items():
            for item in corrections.get(site, ()):
                if int(item['soln']) != int(crd['soln']):
                    continue
                w = 2 * math.pi * fdoy / item['freq']
                for axis in 'XYZ':
                    k = axis.lower()
                    crd[axis] += (math.cos(w) * item['xyz'][k+'cos']
                                  + math.sin(w) * item['xyz'][k+'sin']) * 1e-3
    return result
