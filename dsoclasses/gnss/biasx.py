"""A small, strict reader for satellite *code* OSBs in Bias-SINEX 1.00.

Scope
-----
* Satellite-only OSB rows (no receiver/station-specific biases).
* Code observables, e.g. ``C1C``, whose biases are supplied in nanoseconds.
* GPS, GLONASS, Galileo, BeiDou, QZSS, SBAS, and NavIC PRNs are handled by
  their three-character RINEX identifiers. No constellation-specific sign or
  signal-frequency assumptions are made.
* Phase OSBs, DSBs, and estimates with slopes are **not** applied. This is
  intentional: the teaching positioning exercise uses code measurements.

A key detail: validity intervals belong to *each observable*, not to a whole
satellite. For example, C1C can have a different interval from C2W. Selecting
biases by a satellite-wide interval silently discards otherwise valid data.

The input and output use the time scale stated by the Bias-SINEX product;
no implicit GPS/UTC/other constellation time conversion is performed.
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path
import math

from dsoclasses.time.pyattotime import to_attoseconds

_SPEED_OF_LIGHT = 299_792_458.0  # Exact SI definition [m/s].
_ALLOWED_SYSTEMS = frozenset('GRECJIS')


def _bias_sinex_parse_date(text):
    """Parse YYYY:DDD:SSSSS; refuse undefined/open-ended timestamps.

    Open intervals are allowed by the full specification, but would require
    additional conventions for bias slopes. This deliberately small reader
    makes that unsupported situation explicit rather than guessing.
    """
    try:
        year_s, doy_s, sod_s = text.split(':')
        year, doy, sod = int(year_s), int(doy_s), int(sod_s)
        if len(year_s) != 4:
            raise ValueError('Expected a four-digit year')
        if year == 0 and doy == 0 and sod == 0:
            raise ValueError('Undefined/open Bias-SINEX epoch is not supported')
        if not (1 <= year <= 9999):
            raise ValueError('Year outside datetime range')
        days = 366 if _is_leap(year) else 365
        if not 1 <= doy <= days or not 0 <= sod < 86400:
            raise ValueError('Invalid day-of-year or seconds-of-day')
        return datetime(year, 1, 1) + timedelta(days=doy - 1, seconds=sod)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'Invalid Bias-SINEX epoch {text!r}') from exc


def _is_leap(year):
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


class Bsinex:
    """Read a Bias-SINEX file and query satellite code OSBs.

    Typical student usage::

        osb = Bsinex('daily.BIA')
        bias_ns = osb.sat_bias('G21', epoch)['C1W']
        bias_m = osb.code_bias_m('G21', 'C1W', epoch)
        corrected_code_m = measured_code_m - bias_m

    Bias sign convention: ``measured = ideal + OSB``. Subtraction is the
    convention used in the course notebooks; users must also check their
    precise clock product's signal reference for consistency.
    """

    def __init__(self, filename):
        self._fn = str(Path(filename))
        self._sat_clock_ref_obs = {}
        self._records = None

    def _read_block(self, name):
        """Yield the contents of one +NAME/-NAME block; reject truncation."""
        with open(self._fn, encoding='ascii', errors='replace') as stream:
            for line in stream:
                if line.strip() == f'+{name}':
                    break
            else:
                raise ValueError(f'Bias-SINEX lacks +{name}')
            for line in stream:
                if line.strip() == f'-{name}':
                    return
                yield line
        raise ValueError(f'Bias-SINEX lacks terminating -{name}')

    def parse_description_block(self):
        """Read the optional metadata used to interpret the bias product."""
        for line in self._read_block('BIAS/DESCRIPTION'):
            field = line.lstrip()
            if field.startswith('SATELLITE_CLOCK_REFERENCE_OBSERVABLES'):
                info = line[41:].split()
                if not info or info[0] not in _ALLOWED_SYSTEMS:
                    raise ValueError(f'Invalid satellite clock reference: {line.rstrip()}')
                self._sat_clock_ref_obs[info[0]] = info[1:]
            elif field.startswith('BIAS_MODEL'):
                self._bias_model = line[41:].strip()
            elif field.startswith('APC_MODEL'):
                self._apc_model = line[41:].strip()
            elif field.startswith('TIME_SYSTEM'):
                self._time_system = line[41:].strip()

    def _parse_sat_bias(self):
        # One satellite -> one observable -> a list of (start, stop, bias_ns).
        # Do not group by common satellite interval: C1C/C2W can have different
        # validity windows even within the same file.
        result = defaultdict(lambda: defaultdict(list))
        for line_no, line in enumerate(self._read_block('BIAS/SOLUTION'), start=1):
            if not line.strip() or line.lstrip().startswith(('*', '%')):
                continue
            tokens = line.split()
            if not tokens or tokens[0] != 'OSB':
                continue  # Ignore DSB, ISB, and other supported file records.
            # BIAS/SOLUTION uses fixed columns for an optional station field.
            # Satellite-only records must leave columns 16 through 24 blank.
            if line[15:24].strip():
                continue
            if len(tokens) < 8:
                raise ValueError(f'Malformed satellite OSB at BIAS/SOLUTION line {line_no}')
            _, svn, prn, observable, start, stop, unit, amount = tokens[:8]
            if (len(prn) != 3 or prn[0] not in _ALLOWED_SYSTEMS
                    or not prn[1:].isdigit() or len(svn) != 4):
                continue  # Not an individual identified GNSS satellite.
            if not observable.startswith('C'):
                continue  # Phase OSBs in cycles must not crash a code-only reader.
            if len(observable) != 3:
                raise ValueError(f'Invalid code observable {observable!r} for {prn}')
            if unit.lower() != 'ns':
                raise ValueError(f'Code OSB {prn} {observable}: expected ns, got {unit!r}')
            t0, t1 = _bias_sinex_parse_date(start), _bias_sinex_parse_date(stop)
            if t0 >= t1:
                raise ValueError(f'Invalid interval for {prn} {observable}: {start}..{stop}')
            try:
                bias = float(amount.replace('D', 'E').replace('d', 'E'))
            except ValueError as exc:
                raise ValueError(f'Non-numeric bias for {prn} {observable}') from exc
            if not math.isfinite(bias):
                raise ValueError(f'Non-finite bias for {prn} {observable}')
            result[prn][observable].append((t0, t1, bias))

        if not result:
            raise ValueError('No satellite code OSBs found in Bias-SINEX file')
        # Reject overlapping intervals for the same observable. Selecting a
        # value from two simultaneously applicable estimates is ambiguous.
        for prn, codes in result.items():
            for code, intervals in codes.items():
                intervals.sort(key=lambda item: item[0])
                for prev, current in zip(intervals, intervals[1:]):
                    if current[0] < prev[1]:
                        raise ValueError(f'Overlapping OSB intervals for {prn} {code}')
        return {sat: dict(codes) for sat, codes in result.items()}

    def sat_bias(self, prn, epoch=None):
        """Return available satellite code OSBs in nanoseconds.

        ``epoch`` selects each observable's record independently. If omitted,
        every observable must have exactly one interval, and all intervals
        must overlap; otherwise the request would be ambiguous.

        Missing observables are absent from the returned dictionary: indexing
        them raises KeyError, so bad measurements cannot be silently corrected.
        Results also contain ``from``, ``to`` and ``units`` metadata, retained
        for compatibility with existing notebooks.
        """
        if self._records is None:
            self._records = self._parse_sat_bias()
        codes = self._records.get(prn)
        if not codes:
            raise KeyError(f'No satellite code OSBs for {prn!r}')
        query = None if epoch is None else to_attoseconds(epoch)
        selected = {}
        for observable, intervals in codes.items():
            if query is None:
                if len(intervals) != 1:
                    raise ValueError(f'{prn} has multiple OSB intervals; supply epoch')
                picked = intervals[0]
            else:
                matches = [rec for rec in intervals
                           if to_attoseconds(rec[0]) <= query < to_attoseconds(rec[1])]
                if not matches:
                    continue  # Some signals may not be available at this epoch.
                picked = matches[0]  # Overlaps were rejected when loading.
            selected[observable] = picked
        if not selected:
            raise ValueError(f'Expected one OSB interval for {prn} at {epoch}, found 0')
        valid_from = max(rec[0] for rec in selected.values())
        valid_to = min(rec[1] for rec in selected.values())
        if valid_from >= valid_to:
            raise ValueError(f'{prn} has disjoint OSB intervals; supply epoch')
        return {'from': valid_from, 'to': valid_to, 'units': 'ns',
                **{obs: rec[2] for obs, rec in selected.items()}}

    def code_bias_m(self, satellite, observable, epoch):
        """Get the code OSB in **metres**, ready to subtract from pseudorange.

        Use the same epoch and time system as the observations. The underlying
        satellite-clock reference convention still needs to match the bias
        product; that cannot be inferred from a numerical bias alone.
        """
        if not isinstance(observable, str) or not observable.startswith('C'):
            raise ValueError('code_bias_m expects a RINEX code observable such as C1C')
        try:
            return self.sat_bias(satellite, epoch)[observable] * 1e-9 * _SPEED_OF_LIGHT
        except KeyError as exc:
            raise LookupError(f'No {observable} code OSB for {satellite} at {epoch}') from exc
