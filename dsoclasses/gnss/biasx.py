"""Read satellite observable-specific signal biases from Bias-SINEX.

OSB estimates are expressed in nanoseconds in CODE/IGS files. The public
``sat_bias`` method returns nanoseconds for compatibility with the notebooks;
conversion to metres is ``bias_ns * 1e-9 * speed_of_light``.
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path
import math

from dsoclasses.time.pyattotime import to_attoseconds


def _bias_sinex_parse_date(text):
    """Parse year:day-of-year:second-of-day (including open-ended 000)."""
    year, doy, sod = map(int, text.split(':'))
    if year < 100:
        year += 2000 if year < 80 else 1900
    if not (1 <= doy <= (366 if _is_leap(year) else 365)) or not (0 <= sod <= 86400):
        raise ValueError(f'Invalid Bias-SINEX epoch: {text!r}')
    return datetime(year, 1, 1) + timedelta(days=doy - 1, seconds=sod)


def _is_leap(year):
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


class Bsinex:
    """Satellite OSB reader with explicit validity-period selection.

    When a satellite has exactly one validity interval, ``sat_bias('G21')``
    behaves as before. For multiple intervals, pass ``epoch`` to select the
    applicable solution and prevent silently using the wrong daily bias.
    """
    def __init__(self, filename):
        self._fn = str(Path(filename))
        self._sat_clock_ref_obs = {}
        self._records = None

    def _read_block(self, name):
        """Read a named +BLOCK/-BLOCK pair, rejecting truncated input."""
        start, stop = '+' + name, '-' + name
        with open(self._fn, encoding='ascii', errors='replace') as stream:
            for line in stream:
                if line.strip() == start:
                    break
            else:
                raise ValueError(f'Bias-SINEX lacks {start}')
            for line in stream:
                if line.strip() == stop:
                    return
                yield line
        raise ValueError(f'Bias-SINEX lacks terminating {stop}')

    def parse_description_block(self):
        """Parse clock reference observables; report absent blocks clearly."""
        for line in self._read_block('BIAS/DESCRIPTION'):
            field = line.lstrip()
            if field.startswith('SATELLITE_CLOCK_REFERENCE_OBSERVABLES'):
                info = line[41:].split()
                if not info or info[0] not in 'GRECJIS':
                    raise ValueError(f'Invalid clock reference observables: {line.rstrip()}')
                self._sat_clock_ref_obs[info[0]] = info[1:]
            elif field.startswith('BIAS_MODEL'):
                self._bias_model = line[41:].strip()
            elif field.startswith('APC_MODEL'):
                self._apc_model = line[41:].strip()
            elif field.startswith('TIME_SYSTEM'):
                self._time_system = line[41:].strip()

    def _parse_sat_bias(self):
        per_sat = defaultdict(dict)
        for line in self._read_block('BIAS/SOLUTION'):
            if not line.strip() or line.lstrip().startswith(('*', '%')):
                continue
            values = line.split()
            if not values or values[0] != 'OSB':
                continue
            # Station-specific OSBs are not suitable for the satellite-only
            # correction performed by the notebooks. Check the station field
            # before tokenizing the remaining fixed-width record.
            if line[15:24].strip():
                continue
            if len(values) < 8:
                raise ValueError(f'Malformed OSB line: {line.rstrip()!r}')
            _, svn, prn, obs, start, stop, unit, estimate = values[:8]
            if len(prn) != 3 or prn[0] not in 'GRECJIS' or len(svn) != 4:
                continue
            if unit.lower() != 'ns':
                raise ValueError(f'OSB {prn} {obs} has unsupported unit {unit!r}')
            start_dt = _bias_sinex_parse_date(start)
            end_dt = _bias_sinex_parse_date(stop)
            if end_dt <= start_dt:
                raise ValueError(f'Non-positive OSB interval for {prn} {obs}')
            bias = float(estimate)
            if not math.isfinite(bias):
                raise ValueError(f'Non-finite OSB for {prn} {obs}')
            key = (start_dt, end_dt)
            group = per_sat[prn].setdefault(key, {})
            if obs in group and group[obs] != bias:
                raise ValueError(f'Conflicting duplicate OSB for {prn} {obs} at {key}')
            group[obs] = bias
        if not per_sat:
            raise ValueError('No usable satellite OSB data in Bias-SINEX file')
        return {prn: [{'from': interval[0], 'to': interval[1], 'units': 'ns', **codes}
                       for interval, codes in sorted(groups.items())]
                for prn, groups in per_sat.items()}

    def sat_bias(self, prn, epoch=None):
        """Find satellite code biases valid at ``epoch`` (optional).

        The validity interval is [start, stop), as used in Bias-SINEX. If
        multiple records exist and no epoch is supplied, ambiguity is an
        error. This is preferable to silently taking the first day's bias.
        """
        if self._records is None:
            self._records = self._parse_sat_bias()
        entries = self._records.get(prn)
        if not entries:
            raise KeyError(f'No satellite OSB data for {prn}')
        if epoch is None:
            if len(entries) != 1:
                raise ValueError(f'{prn} has {len(entries)} OSB intervals; supply epoch')
            return entries[0].copy()
        # The Bias-SINEX header stores native datetimes, while the RINEX
        # epoch may be an attodatetime. Compare normalized integer epochs,
        # preserving the [start, stop) validity convention.
        query = to_attoseconds(epoch)
        matches = [record for record in entries
                   if to_attoseconds(record['from']) <= query < to_attoseconds(record['to'])]
        if len(matches) != 1:
            raise ValueError(f'Expected one OSB interval for {prn} at {epoch}, found {len(matches)}')
        return matches[0].copy()
