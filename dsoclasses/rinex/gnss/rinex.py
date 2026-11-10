"""Small, strict RINEX 3 observation reader used by the GNSS notebooks.

This reader intentionally targets observation data, not navigation data. It
preserves the original ``GnssRinex`` / ``DataBlock`` / ``fetch`` API. Incomplete
measurements are omitted, while malformed file structure raises a contextual
exception instead of silently shifting the observation columns.

RINEX event epochs (flags 2--6) contain special records, not normal satellite
observations. They are skipped to avoid interpreting those records as ranges.
"""
from __future__ import annotations

from datetime import datetime
from decimal import Decimal, InvalidOperation
import math
from pathlib import Path

from dsoclasses.time.pyattotime import (
    from_attoseconds, seconds_to_attoseconds, to_attoseconds,
)


_SYSTEMS = {
    'gps': 'G', 'g': 'G', 'glonass': 'R', 'glo': 'R', 'r': 'R',
    'galileo': 'E', 'gal': 'E', 'e': 'E', 'beidou': 'C', 'bds': 'C',
    'c': 'C', 'qzss': 'J', 'j': 'J', 'irnss': 'I', 'navic': 'I',
    'navic/irnss': 'I', 'i': 'I', 'sbas': 'S', 's': 'S',
}


def sysstr2sysid(sat_sys):
    """Normalize a RINEX constellation name/code to its letter identifier."""
    if not isinstance(sat_sys, str):
        raise TypeError('Satellite system must be a string')
    try:
        return _SYSTEMS[sat_sys.strip().lower()]
    except KeyError as exc:
        raise ValueError(f'Unknown GNSS satellite system: {sat_sys!r}') from exc


def fetch(observables, *preferred_codes):
    """Find the first *valid, finite* observation in preference order.

    A RINEX record may contain a code with an empty/non-finite ``value``.
    Presence of the dictionary key alone does not mean an observation exists.
    """
    for code in preferred_codes:
        item = observables.get(code)
        if not isinstance(item, dict):
            continue
        value = item.get('value')
        if value is not None:
            try:
                if math.isfinite(float(value)):
                    return item, code
            except (ValueError, TypeError):
                pass
    raise LookupError(f'No valid observable found among {preferred_codes!r}')


def fetchv(observables, *preferred_codes):
    """As :func:`fetch`, but return only the observation's numeric value."""
    item, _ = fetch(observables, *preferred_codes)
    return float(item['value'])


class GnssRinex:
    """Repeatedly iterable RINEX 3 observation file reader.

    Attributes in the header are accessible by their historical names, e.g.
    ``marker_name``, ``antenna_type`` and ``time_first_obs``. Times are kept
    as printed in the file's time system: GPS/UTC conversion is **not** done.
    """

    @staticmethod
    def resolve_date(date_text):
        """Parse a RINEX epoch, retaining fractional seconds up to attoseconds."""
        parts = date_text.split()
        if len(parts) < 6:
            raise ValueError(f'Invalid RINEX epoch: {date_text!r}')
        try:
            year, month, day, hour, minute = (int(v) for v in parts[:5])
            seconds = Decimal(parts[5])
            if not seconds.is_finite() or seconds < 0 or seconds >= 60:
                raise ValueError('Invalid seconds-of-minute')
            whole = int(seconds)
            fractional_as = seconds_to_attoseconds(seconds - whole)
            base = datetime(year, month, day, hour, minute, whole)
        except (ValueError, InvalidOperation, OverflowError) as exc:
            raise ValueError(f'Invalid RINEX epoch: {date_text!r}') from exc
        # Keep the complete fractional seconds printed in RINEX. In
        # particular, do not truncate to nanoseconds or native microseconds.
        # Both SP3 interpolation and RINEX parsing use the same time origin.
        return from_attoseconds(to_attoseconds(base) + fractional_as)

    def __init__(self, filename):
        self.filename = str(Path(filename))
        self._stream = None
        self._pending_line = None
        self.parse_header(self.filename)

    def parse_header(self, filename):
        """Validate the fixed-column RINEX header and observation-code lists."""
        self.obscodes = {}
        self.glo_slots = {}
        self.marker_name = None
        self.antenna_type = None
        self.time_first_obs = None
        self.time_sys = None
        self.xapprox = self.yapprox = self.zapprox = None
        with open(filename, 'r', encoding='ascii', errors='replace') as stream:
            line = stream.readline()
            if not line:
                raise ValueError(f'Empty RINEX file: {filename}')
            if line[60:80].strip() != 'RINEX VERSION / TYPE':
                raise ValueError('Missing RINEX VERSION / TYPE header record')
            try:
                self.version = float(line[:9])
            except ValueError as exc:
                raise ValueError('Invalid RINEX version') from exc
            if not (3 <= self.version < 4):
                raise ValueError(f'RINEX version {self.version} not supported; expected 3.x')
            self.type = line[20:21]
            self.system = line[40:41]
            if self.type != 'O':
                raise ValueError('Only RINEX 3 observation files are supported')
            for line_no, line in enumerate(stream, start=2):
                label = line[60:80].strip()
                if label == 'END OF HEADER':
                    break
                if label == 'MARKER NAME':
                    self.marker_name = line[:60].strip()
                elif label == 'MARKER NUMBER':
                    self.marker_number = line[:60].strip()
                elif label == 'MARKER TYPE':
                    self.marker_type = line[:60].strip()
                elif label == 'OBSERVER / AGENCY':
                    self.observer = line[:20].strip()
                    self.agency = line[20:60].strip()
                elif label == 'REC # / TYPE / VERS':
                    self.receiver_number = line[:20].strip()
                    self.receiver_type = line[20:40].strip()
                    self.receiver_version = line[40:60].strip()
                elif label == 'ANT # / TYPE':
                    self.antenna_number = line[:20].strip()
                    self.antenna_type = line[20:40].strip()
                elif label == 'APPROX POSITION XYZ':
                    values = line[:60].split()
                    if len(values) != 3:
                        raise ValueError(f'Invalid approximate XYZ at header line {line_no}')
                    self.xapprox, self.yapprox, self.zapprox = map(float, values)
                elif label == 'ANTENNA: DELTA H/E/N':
                    self.dh, self.de, self.dn = map(float, line[:60].split())
                elif label == 'SYS / # / OBS TYPES':
                    system = line[:1].strip()
                    if not system:
                        raise ValueError(f'Unexpected observation-code continuation at line {line_no}')
                    if system in self.obscodes:
                        raise ValueError(f'Duplicate observation-code declaration: {system}')
                    try:
                        count = int(line[3:6])
                    except ValueError as exc:
                        raise ValueError(f'Invalid observation count at line {line_no}') from exc
                    codes = line[7:60].split()
                    while len(codes) < count:
                        continuation = next(stream, None)
                        if continuation is None or continuation[60:80].strip() != 'SYS / # / OBS TYPES':
                            raise ValueError(f'Truncated observation codes for {system}')
                        if continuation[:1].strip() not in ('', system):
                            raise ValueError(f'Unexpected system in observation-code continuation: {system}')
                        codes.extend(continuation[7:60].split())
                    if count == 0 or len(codes) != count or len(set(codes)) != count:
                        raise ValueError(f'Invalid observation-code list for {system}')
                    self.obscodes[system] = codes
                elif label == 'TIME OF FIRST OBS':
                    self.time_first_obs = self.resolve_date(line[:44])
                    self.time_sys = line[44:60].strip()
                elif label == 'TIME OF LAST OBS':
                    self.time_last_obs = self.resolve_date(line[:44])
                    self.time_last_obs_sys = line[44:60].strip()
                elif label == 'GLONASS SLOT / FRQ #':
                    # The first line contains a satellite count; continuation
                    # records may omit it. Merge all pairs as found.
                    fields = line[4:60].split()
                    if len(fields) % 2:
                        raise ValueError(f'Incomplete GLONASS slot pair at {line_no}')
                    for sat, slot in zip(fields[::2], fields[1::2]):
                        self.glo_slots[sat] = int(slot)
            else:
                raise ValueError('RINEX header is missing END OF HEADER')
        if not self.obscodes:
            raise ValueError('No observation types were declared in the RINEX header')

    def glonass_slot(self, prn):
        try:
            return self.glo_slots[prn]
        except KeyError as exc:
            raise KeyError(f'No GLONASS frequency slot for {prn}') from exc

    def approx_cartesian(self):
        if self.xapprox is None:
            raise ValueError('RINEX header has no approximate receiver position')
        return [self.xapprox, self.yapprox, self.zapprox]

    class DataBlock:
        """One ordinary observation epoch, preserving notebook API."""
        def __init__(self, values):
            self.dct = values

        def t(self):
            return self.dct['epoch']

        def flag(self):
            return self.dct['flag']

        def nsats(self):
            return self.dct['num_sats']

        def satellite(self, satid):
            return self.dct[satid]

        def filter_satellite_system(self, sat_sys, include_block_info=True):
            system = sysstr2sysid(sat_sys)
            result = {key: value for key, value in self if key[0] == system}
            if include_block_info:
                result = {'epoch': self.t(), 'flag': self.flag(),
                          'num_sats': len(result), **result}
            return GnssRinex.DataBlock(result)

        def __iter__(self):
            return ((key, value) for key, value in self.dct.items()
                    if key not in ('epoch', 'flag', 'num_sats'))

    def close(self):
        """Close the active observation stream; safe to call repeatedly."""
        if self._stream is not None:
            self._stream.close()
            self._stream = None
        self._pending_line = None

    def __iter__(self):
        self.close()  # The notebooks deliberately iterate over the file often.
        self._stream = open(self.filename, 'r', encoding='ascii', errors='replace')
        for line in self._stream:
            if line[60:80].strip() == 'END OF HEADER':
                return self
        self.close()
        raise ValueError('RINEX header is missing END OF HEADER')

    def _readline(self):
        if self._pending_line is not None:
            line, self._pending_line = self._pending_line, None
            return line
        return self._stream.readline()

    def _satellite_line(self, satid, first_line):
        """Read one satellite's fixed-width observation fields.

        RINEX 3 commonly places all fields on one long line, but also allows
        continuation records. Missing *trailing* blank fields are permitted.
        """
        codes = self.obscodes.get(satid[:1])
        if codes is None:
            raise ValueError(f'Satellite {satid} has no declared observation types')
        payload = first_line[3:].rstrip('\r\n')
        width = 16 * len(codes)
        while len(payload) < width:
            line = self._stream.readline()
            if not line:
                break
            if line.startswith('   '):
                payload += line[3:].rstrip('\r\n')
            else:
                self._pending_line = line
                break
        observations = {}
        for i, code in enumerate(codes):
            field = payload[16 * i:16 * (i + 1)].ljust(16)
            if not field[:14].strip():
                continue
            try:
                value = float(field[:14])
                if not math.isfinite(value):
                    continue
                lli = int(field[14]) if field[14].strip() else None
                ssi = int(field[15]) if field[15].strip() else None
            except ValueError as exc:
                raise ValueError(f'Invalid {satid} observation {code}: {field!r}') from exc
            observations[code] = {'value': value, 'lli': lli, 'ssi': ssi}
        return observations

    def __next__(self):
        if self._stream is None:
            raise StopIteration
        while True:
            line = self._readline()
            if not line:
                self.close()
                raise StopIteration
            if not line.startswith('>'):
                self.close()
                raise ValueError(f'Expected RINEX epoch record, got {line[:30]!r}')
            try:
                epoch = self.resolve_date(line[1:29])
                flag = int(line[29:32])
                count = int(line[32:35])
                if count < 0:
                    raise ValueError('Negative observation/special-record count')
            except ValueError as exc:
                self.close()
                raise ValueError(f'Malformed RINEX epoch record: {line[:40]!r}') from exc
            if flag > 1:
                # Event and cycle-slip records are not GNSS measurements.
                for _ in range(count):
                    if not self._stream.readline():
                        self.close()
                        raise ValueError(f'Truncated RINEX special-event block at {epoch}')
                continue
            result = {'epoch': epoch, 'flag': flag, 'num_sats': count}
            for _ in range(count):
                row = self._readline()
                if len(row) < 3 or row.startswith('>'):
                    self.close()
                    raise ValueError(f'Truncated observation block at {epoch}')
                satellite = row[:3]
                if satellite in result:
                    self.close()
                    raise ValueError(f'Duplicate satellite {satellite} at {epoch}')
                result[satellite] = self._satellite_line(satellite, row)
            return self.DataBlock(result)

    def __del__(self):
        self.close()
