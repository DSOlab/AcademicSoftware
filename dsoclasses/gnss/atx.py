"""Receiver antenna phase-centre corrections from ANTEX 1.4 files.

All public PCO/PCV values are returned in **metres**. ANTEX stores them in
millimetres. Zenith angles in ANTEX are degrees; elevations in this API are
radians. The helper handles NOAZI (azimuth-independent) patterns only.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

# ANTEX frequency identifiers differ from PRNs. In particular GPS L5 is G05.
_ATX_FREQS = {
    'G': {'L1': 'G01', 'L2': 'G02', 'L5': 'G05'},
    'R': {'G1': 'R01', 'G2': 'R02', 'G3': 'R03', 'L1': 'R01', 'L2': 'R02'},
    'E': {'E1': 'E01', 'E5a': 'E05', 'E5b': 'E07', 'E5(E5a+E5b)': 'E08', 'E6': 'E06'},
    'C': {'B1C': 'C01', 'B1I': 'C02', 'B2a': 'C05', 'B2b': 'C07', 'B3I': 'C06', 'E1': 'C01', 'E2': 'C02', 'E5b': 'C07', 'E6': 'C06'},
    'J': {'L1': 'J01', 'L2': 'J02', 'L5': 'J05', 'LEX': 'J06'},
    'S': {'L1': 'S01', 'L5': 'S05'},
    'I': {'L5': 'I05', 'S': 'I09'},
}


def _to_atx_freq(system, frequency):
    try:
        return _ATX_FREQS[system.upper()[0]][frequency]
    except (IndexError, KeyError, AttributeError) as exc:
        raise ValueError(f'Unknown ANTEX frequency {system!r} / {frequency!r}') from exc


class ReceiverAntennaPcv:
    """Antenna NEU offsets and NOAZI variation curves in SI metres."""
    def __init__(self, antenna):
        self.antenna = antenna
        self.pcv = {}

    def add_freq(self, freq, neu, z1z2dz, vals):
        """Add one frequency; caller passes offsets/variations in metres."""
        offset = np.asarray(neu, dtype=float)
        angles = np.asarray(z1z2dz, dtype=float)
        variation = np.asarray(vals, dtype=float)
        if offset.shape != (3,) or angles.shape != (3,) or variation.ndim != 1:
            raise ValueError(f'Invalid antenna calibration dimensions for {freq}')
        z1, z2, dz = angles
        if not np.all(np.isfinite(offset)) or not np.all(np.isfinite(variation)):
            raise ValueError(f'Non-finite antenna calibration for {freq}')
        if dz <= 0 or z2 < z1 or z1 < 0 or z2 > 90 or not np.all(np.isfinite(angles)):
            raise ValueError(f'Invalid zenith-angle grid for {freq}')
        expected = int(round((z2 - z1) / dz)) + 1
        if expected != len(variation) or not np.isclose(z1 + (expected - 1) * dz, z2):
            raise ValueError(f'Wrong NOAZI sample count for {freq}: expected {expected}, got {len(variation)}')
        if freq in self.pcv:
            raise ValueError(f'Duplicate ANTEX calibration for {freq}')
        self.pcv[freq] = {'neu': offset, 'z1z2dz': angles, 'pcv': variation}

    @staticmethod
    def _resolve_freq(args):
        if len(args) == 1:
            return args[0]
        if len(args) == 2:
            return _to_atx_freq(*args)
        raise TypeError('Frequency should be a code or a system/frequency pair')

    def pco(self, *args):
        """Return phase-centre offset as [north, east, up] in metres."""
        return self.pcv[self._resolve_freq(args)]['neu'].copy()

    def pcv_hgt(self, *args):
        """Interpolate phase-centre variation [m] at elevation [rad].

        Exact grid endpoints are valid, including the zenith itself. We do
        not silently extrapolate outside the ANTEX zenith-angle coverage.
        """
        if len(args) == 2:
            freq, elevation = args
        elif len(args) == 3:
            system, frequency, elevation = args
            freq = _to_atx_freq(system, frequency)
        else:
            raise TypeError('pcv_hgt expects (freq, elevation) or (system, freq, elevation)')
        record = self.pcv[freq]
        if not math.isfinite(float(elevation)):
            raise ValueError('Elevation must be finite')
        za = 90.0 - math.degrees(elevation)
        z1, z2, dz = record['z1z2dz']
        if za < z1 - 1e-9 or za > z2 + 1e-9:
            raise ValueError(f'Elevation outside ANTEX coverage for {freq} (zenith {z1}..{z2} deg)')
        knots = z1 + np.arange(len(record['pcv'])) * dz
        return float(np.interp(np.clip(za, z1, z2), knots, record['pcv']))


    def range_correction_m(self, freq, azimuth, elevation):
        """Receiver antenna correction from phase centre to ARP [metres].

        ``azimuth`` is clockwise from north, ``elevation`` above horizon;
        both are in radians. The receiver PCO vector is ANTEX North/East/Up.
        Projecting onto the receiver-to-satellite line of sight gives the
        geometrical difference between distances measured at the ARP and
        at the frequency-dependent antenna phase centre.

        The returned ``dot(PCO, LOS) + PCV`` follows the sign convention used
        by the course pseudorange notebook: add it to a measured code range
        to refer that range to the antenna reference point (ARP). It does
        NOT include satellite-antenna offsets, nor azimuth-dependent PCV.
        """
        if not all(math.isfinite(float(x)) for x in (azimuth, elevation)):
            raise ValueError('Azimuth and elevation must be finite radians')
        if not 0 <= elevation <= math.pi/2:
            raise ValueError('Antenna correction requires elevation in 0..90 degrees')
        neu_los = np.array([math.cos(elevation) * math.cos(azimuth),
                            math.cos(elevation) * math.sin(azimuth),
                            math.sin(elevation)])
        return float(np.dot(self.pco(freq), neu_los) + self.pcv_hgt(freq, elevation))


class Atx:
    """Read receiver NOAZI calibration blocks from an ANTEX file."""
    def __init__(self, filename):
        self.filename = str(Path(filename))
        with open(self.filename, encoding='ascii', errors='replace') as stream:
            first = stream.readline()
            if first[60:80].strip() != 'ANTEX VERSION / SYST':
                raise ValueError('Not an ANTEX file (missing version header)')
            self.version = float(first[:8])
            self.sat_sys = first[20:21]
            self.variation_type = None
            for line in stream:
                label = line[60:80].strip()
                if label == 'PCV TYPE / REFANT':
                    self.variation_type = line[:1]
                elif label == 'END OF HEADER':
                    break
            else:
                raise ValueError('ANTEX header lacks END OF HEADER')
            if self.variation_type != 'A':
                raise ValueError('Only absolute antenna calibrations (PCV type A) are supported')

    def goto_antenna(self, antenna):
        """Compatibility method: return a stream at the requested antenna.

        The *caller* must close the returned stream. ``get_noazi`` below uses
        its own context manager and does not leave files open.
        """
        stream = open(self.filename, encoding='ascii', errors='replace')
        for line in stream:
            if line[60:80].strip() == 'START OF ANTENNA':
                antenna_line = stream.readline()
                if self._antenna_matches(antenna_line[:20], antenna):
                    return stream
        stream.close()
        return None

    @staticmethod
    def _antenna_matches(file_name, requested):
        name = file_name.strip()
        requested = requested.strip()
        return name == requested or (name == f'{requested:<16}NONE'.strip()
                                     and len(requested) <= 16)

    @staticmethod
    def _read_mm_values(line):
        """Decode an entire ANTEX NOAZI row, returning metres.

        ANTEX labels on control/header records occupy columns 61--80, but
        **NOAZI is a data record**, not an 80-column header record. Its values
        begin in column 9 and may continue far beyond column 80 (e.g. 19
        samples for zenith angles from 0 to 90 degrees in steps of 5 degrees).
        Slicing at column 60 truncates a valid F8.2 field and can leave a
        lone minus sign, causing the original ``float('-')`` exception.

        The specification uses 8-character F8.2 fields. We decode those first
        so adjacent negative numbers remain separable. For small hand-written
        examples that omit fixed-width padding, we also accept whitespace-
        separated values. Either way, missing/invalid and non-finite values
        must be rejected rather than silently inserted into a calibration.
        """
        payload = line[8:].rstrip('\r\n')
        if not payload.strip():
            raise ValueError('ANTEX NOAZI row contains no calibration samples')

        columns = [payload[i:i + 8].strip()
                   for i in range(0, len(payload), 8)]
        try:
            mm_values = [float(value) for value in columns if value]
        except ValueError:
            # Tolerate simple, nonstandard test/export files while still
            # rejecting malformed numeric tokens in the fallback below.
            try:
                mm_values = [float(value) for value in payload.split()]
            except ValueError as exc:
                raise ValueError(
                    'ANTEX NOAZI contains an invalid millimetre value '
                    f'in {payload!r}'
                ) from exc

        if not mm_values or not all(math.isfinite(value) for value in mm_values):
            raise ValueError('ANTEX NOAZI contains missing or non-finite values')
        return [value * 1e-3 for value in mm_values]

    def _get_noazi(self, antenna, frequencies):
        requested = list(dict.fromkeys(frequencies))
        if not requested:
            raise ValueError('At least one frequency must be requested')
        stream = self.goto_antenna(antenna)
        if stream is None:
            raise LookupError(f'Antenna {antenna!r} not found in {self.filename}')

        result = ReceiverAntennaPcv(antenna)
        zenith_grid = None
        active_freq = None
        neu = None
        values = []
        try:
            for line in stream:
                label = line[60:80].strip()
                if label == 'END OF ANTENNA':
                    break
                if label == 'ZEN1 / ZEN2 / DZEN':
                    zenith_grid = [float(v) for v in line[:60].split()]
                elif label == 'START OF FREQUENCY':
                    active_freq = line[3:6]
                    neu = None
                    values = []
                elif label == 'NORTH / EAST / UP' and active_freq in requested:
                    neu = [float(line[i:i+10]) * 1e-3 for i in (0, 10, 20)]
                elif label == 'END OF FREQUENCY':
                    if active_freq in requested:
                        if zenith_grid is None or neu is None:
                            raise ValueError(f'Incomplete calibration for {active_freq} on {antenna}')
                        result.add_freq(active_freq, neu, zenith_grid, values)
                    active_freq = None
                elif active_freq in requested and line.startswith('   NOAZI'):
                    # NOAZI is a potentially very long data line; the parser
                    # must not apply the 60-character label-field limit.
                    values.extend(self._read_mm_values(line))
                elif (active_freq in requested and neu is not None and values
                      and zenith_grid is not None and not line[:8].strip()
                      and len(values) < int(round((zenith_grid[1] - zenith_grid[0])
                                                       / zenith_grid[2])) + 1):
                    # Normally all NOAZI samples occupy one line. Accept
                    # continuation lines from nonstandard writers only until
                    # the declared grid is complete. Unlike azimuth-dependent
                    # PCV rows, continuation lines have blank columns 1..8.
                    # A long continuation can itself exceed column 60.
                    values.extend(self._read_mm_values(line))
            else:
                raise ValueError(f'ANTEX antenna {antenna!r} lacks END OF ANTENNA')
        finally:
            stream.close()
        missing = [code for code in requested if code not in result.pcv]
        if missing:
            raise LookupError(f'No NOAZI calibration for antenna {antenna!r}, frequencies {missing}')
        return result

    def get_noazi(self, antenna, freq_list):
        """Return the selected NOAZI curves; accept string, list, or map."""
        if isinstance(freq_list, str):
            freqs = [freq_list]
        elif isinstance(freq_list, dict):
            freqs = []
            for system, bands in freq_list.items():
                if isinstance(bands, str):
                    bands = [bands]
                freqs.extend(_to_atx_freq(system, band) for band in bands)
        elif isinstance(freq_list, (list, tuple)):
            freqs = list(freq_list)
        else:
            raise TypeError('freq_list must be a string, sequence or constellation map')
        return self._get_noazi(antenna, freqs)
