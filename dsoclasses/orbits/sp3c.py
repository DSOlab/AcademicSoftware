"""Read SP3-c/SP3-d precise-orbit files without assuming a GNSS constellation.

The public ``Sp3`` interface is deliberately small and preserves the historical
methods ``get_system_pos`` and ``get_satellite``.  Internally, one parser handles
both a single LEO satellite and a mixed GNSS constellation.  Velocity records
are optional; the interpolator currently uses positions and clocks only.

SP3 file units (see https://files.igs.org/pub/data/format/sp3c.txt):
    position:       km                  -> SI: metres
    clock:          microseconds        -> SI: seconds
    velocity:       decimetres / second -> SI: metres / second
    clock rate:     1e-4 microsec/sec   -> SI: seconds / second (factor 1e-10)

SP3 time labels are *not* automatically converted between GPS/TAI/UTC/etc.
The file's ``time_sys`` identifies the scale of every epoch.  A Python naive
``datetime`` supplied to the interpolator must represent that same scale.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import re
from typing import Iterator

import numpy as np


# SP3 uses these all-zero / all-nine sentinels for unavailable observations.
_MISSING_CLOCK = 999999.0
_VALID_SAT_ID = re.compile(r"^[A-Za-z][0-9]{2}$")


def _parse_epoch(line: str, *, line_number: int) -> np.datetime64:
    """Parse the epoch's fixed-layout date into an exact nanosecond timestamp.

    SP3-c prints eight decimal places for seconds, more than the six supported
    by ``datetime``.  Convert the fractional part to integer nanoseconds so
    valid SP3 timestamps are not silently rounded to microseconds.
    """
    try:
        # Splitting only the date fields is safe even if their padding varies.
        year, month, day, hour, minute, seconds = line[3:31].split()
        sec, dot, fractional = seconds.partition(".")
        if len(fractional) > 9 or (fractional and not fractional.isdigit()):
            raise ValueError("fractional seconds exceed nanosecond precision")
        if not 0 <= int(sec) < 60:
            # Leap-second labels cannot be represented by numpy.datetime64;
            # do not silently map a leap second onto the following minute.
            raise ValueError("leap-second or invalid second label")
        base = datetime(int(year), int(month), int(day), int(hour), int(minute), int(sec))
        ns = int(fractional.ljust(9, "0")) if dot else 0
        return np.datetime64(base, "ns") + np.timedelta64(ns, "ns")
    except (ValueError, IndexError, OverflowError) as exc:
        raise ValueError(f"Invalid SP3 epoch at line {line_number}: {line.rstrip()!r}") from exc


def _optional_float(field: str, *, name: str, line_number: int) -> float:
    """Read one optional fixed-width number, mapping blanks to NaN."""
    if not field.strip():
        return float("nan")
    try:
        return float(field)
    except ValueError as exc:
        raise ValueError(f"Invalid {name} at SP3 line {line_number}: {field!r}") from exc


def _optional_int(field: str) -> int | None:
    """SP3 accuracy exponents are optional and need not be present in a file."""
    return int(field) if field.strip() else None


class Sp3:
    """Read a position (``#cP``) or position/velocity (``#cV``) SP3 product.

    Parameters
    ----------
    fn:
        Path to an uncompressed SP3-c or SP3-d text file.  The reader accepts
        the standard multi-record header or a shorter, otherwise valid header.

    Notes
    -----
    ``get_system_pos`` produces a mapping from ``numpy.datetime64[ns]`` epochs
    to satellite observations.  This is a deliberate improvement over Python
    ``datetime`` keys, which cannot represent all eight SP3 fractional digits.
    """

    def __init__(self, fn: str | Path):
        self.fn = str(fn)
        self.sat_ids: list[str] = []
        self.sat_acc: list[int | None] = []
        self.time_sys = "UNKNOWN"
        self.base_posvel: float | None = None
        self.base_clkrate: float | None = None
        self.end_of_head = 0
        self._header_lines = 0

        # Scan until the first epoch rather than assuming exactly four comment
        # lines or exactly five +/- records.  SP3-d permits variable lengths.
        with open(self.fn, "r", encoding="ascii") as fin:
            first = fin.readline()
            if not first.startswith(("#c", "#d", "#C", "#D")):
                raise ValueError(f"Not an SP3-c/SP3-d file: {self.fn!r}")
            if len(first) < 40 or first[2:3].upper() not in ("P", "V"):
                raise ValueError("Malformed SP3 first header record")
            self.pos_vel = first[2].upper()
            self.t_start = _parse_epoch("*  " + first[3:31], line_number=1)
            try:
                self.num_epochs = int(first[32:39])
            except ValueError as exc:
                raise ValueError("Invalid number of epochs in SP3 header") from exc
            self.crd_system = first[46:51].strip()
            self.orb_type = first[52:55].strip()
            self.agency = first[56:60].strip()
            self._header_lines = 1

            sat_count: int | None = None
            while True:
                offset = fin.tell()
                line = fin.readline()
                if not line:
                    raise ValueError("SP3 file has no epoch records")
                self._header_lines += 1
                if line.startswith("* "):
                    self.end_of_head = offset
                    self._header_lines -= 1
                    break
                if line.startswith("+ "):
                    if sat_count is None:
                        count_field = line[3:6].strip()
                        sat_count = int(count_field) if count_field else None
                    for start in range(9, min(len(line.rstrip("\r\n")), 60), 3):
                        candidate = line[start : start + 3].strip().upper()
                        if _VALID_SAT_ID.fullmatch(candidate):
                            self.sat_ids.append(candidate)
                elif line.startswith("++"):
                    for start in range(9, min(len(line.rstrip("\r\n")), 60), 3):
                        self.sat_acc.append(_optional_int(line[start : start + 3]))
                elif line.startswith("%c"):
                    if self.time_sys == "UNKNOWN":
                        code = line[9:12].strip().upper()
                        if code and code != "CCC":
                            self.time_sys = code
                elif line.startswith("%f"):
                    if self.base_posvel is None:
                        self.base_posvel = _optional_float(line[3:13], name="position accuracy base", line_number=self._header_lines)
                        self.base_clkrate = _optional_float(line[14:26], name="clock accuracy base", line_number=self._header_lines)
                elif line.startswith(("##", "%i", "/*")):
                    pass  # Metadata not needed to calculate position and clock.
                else:
                    raise ValueError(f"Unexpected SP3 header record: {line.rstrip()!r}")

        # Header + records specify the set of satellites.  An incomplete list
        # is not fatal: some single-satellite files omit placeholder records.
        self.sat_ids = list(dict.fromkeys(self.sat_ids))
        if sat_count is not None and sat_count != len(self.sat_ids):
            raise ValueError(
                f"SP3 header declares {sat_count} satellites but lists {len(self.sat_ids)}"
            )
        self.sat_acc = self.sat_acc[: len(self.sat_ids)]

    @staticmethod
    def _record(line: str, *, line_number: int, toSI: bool) -> tuple[str, dict]:
        """Parse an SP3 P or V record, including optional flags and accuracies.

        Fixed-width slices are essential: negative values can abut adjacent
        fields and whitespace splitting is not guaranteed to work.
        """
        record_type = line[0]
        satellite = line[1:4].strip().upper()
        if not _VALID_SAT_ID.fullmatch(satellite):
            raise ValueError(f"Invalid satellite ID at line {line_number}: {satellite!r}")
        if len(line.rstrip("\r\n")) < 46:
            raise ValueError(f"Truncated SP3 {record_type} record at line {line_number}")

        values = [
            _optional_float(line[i : i + 14], name="orbit component", line_number=line_number)
            for i in (4, 18, 32)
        ]
        clock = _optional_float(line[46:60], name="clock component", line_number=line_number)
        if record_type == "P":
            # The SP3 missing-position indicator is (0, 0, 0).  Zero in a
            # single component is perfectly legitimate; do not reject it.
            if not np.all(np.isfinite(values)) or all(v == 0 for v in values):
                values = [float("nan")] * 3
            if not np.isfinite(clock) or abs(clock) >= _MISSING_CLOCK:
                clock = float("nan")
            if toSI:
                values = [v * 1000.0 for v in values]
                clock *= 1e-6
            # SP3-c flags: clock event (75), clock predicted (76), orbit
            # maneuver (79), orbit predicted (80).  Preserve their letters.
            flags = "".join(line[i] if i < len(line) else " " for i in (74, 75, 78, 79)).strip()
            result = dict(zip(("x", "y", "z"), values))
            result["c"] = clock
            result["events"] = flags
        else:
            if not np.all(np.isfinite(values)):
                values = [float("nan")] * 3
            if not np.isfinite(clock) or abs(clock) >= _MISSING_CLOCK:
                clock = float("nan")
            if toSI:
                values = [v * 0.1 for v in values]
                clock *= 1e-10  # The SP3 clock-rate unit is 10^-4 microsec/s.
            result = dict(zip(("vx", "vy", "vz"), values))
            result["cr"] = clock

        # Read *known* optional accuracy fields by column, without assuming
        # that every SP3 producer supplies them.  Blank fields remain None.
        if len(line.rstrip("\r\n")) >= 69:
            prefix = "" if record_type == "P" else "v"
            for name, i in zip(("xsdev", "ysdev", "zsdev"), (61, 64, 67)):
                result[prefix + name] = _optional_int(line[i : i + 2])
        if len(line.rstrip("\r\n")) >= 73:
            result["csdev" if record_type == "P" else "crsdev"] = _optional_int(line[70:73])
        return satellite, result

    def _iter_epochs(self, *, toSI: bool) -> Iterator[tuple[np.datetime64, dict[str, dict]]]:
        """Yield one epoch at a time, grouping satellites and optional V lines.

        We check ordering, duplicates and the expected epoch count.  Correlation
        records (EP/EV) are recognized but intentionally ignored: they are not
        needed for position/clock interpolation.
        """
        current_time: np.datetime64 | None = None
        observations: dict[str, dict] = {}
        epochs_seen = 0
        found_eof = False

        with open(self.fn, "r", encoding="ascii") as fin:
            fin.seek(self.end_of_head)
            for line_number, line in enumerate(fin, start=self._header_lines + 1):
                if not line.strip():
                    continue
                if line.startswith("EOF"):
                    found_eof = True
                    break
                if line.startswith("* "):
                    new_time = _parse_epoch(line, line_number=line_number)
                    if current_time is not None:
                        if new_time <= current_time:
                            raise ValueError("SP3 epochs must be strictly increasing")
                        yield current_time, observations
                    current_time, observations = new_time, {}
                    epochs_seen += 1
                    continue
                if line.startswith(("EP", "EV")):
                    continue  # Optional correlation records have no P/V satellite ID.
                if line.startswith(("P", "V")):
                    if current_time is None:
                        raise ValueError(f"Satellite record before first epoch at line {line_number}")
                    satellite, record = self._record(line, line_number=line_number, toSI=toSI)
                    if self.sat_ids and satellite not in self.sat_ids:
                        raise ValueError(f"Satellite {satellite} at line {line_number} not listed in header")
                    entry = observations.setdefault(satellite, {})
                    if line[0] == "P":
                        if "x" in entry:
                            raise ValueError(f"Duplicate position record for {satellite} at {current_time}")
                    elif "vx" in entry:
                        raise ValueError(f"Duplicate velocity record for {satellite} at {current_time}")
                    entry.update(record)
                    continue
                raise ValueError(f"Unexpected SP3 record at line {line_number}: {line.rstrip()!r}")

        if not found_eof:
            raise ValueError("SP3 data ended without the required EOF record")
        if current_time is not None:
            yield current_time, observations
        if self.num_epochs != epochs_seen:
            raise ValueError(f"SP3 declares {self.num_epochs} epochs, found {epochs_seen}")

    def get_system_pos(self, satsys, toSI: bool = True) -> dict:
        """Return positions and clocks for the selected satellite systems.

        The nested value is ``[x, y, z, clock, event_flags]``.  All values are
        in SI units by default.  Missing positions/clocks are represented by
        ``NaN`` and can be handled independently by the interpolator.

        ``satsys`` may be ``['G', 'R']``, ``'G'`` or ``'GR'``.  An empty selection
        is an error rather than silently returning an empty interpolator.
        """
        entries = [satsys] if isinstance(satsys, str) else list(satsys)
        if not entries or any(not isinstance(item, str) or not item for item in entries):
            raise ValueError("satsys must contain satellite-system letters such as 'G' or 'L'")
        systems = {letter.upper() for item in entries for letter in item}
        if not systems or any(len(s) != 1 or not s.isalpha() for s in systems):
            raise ValueError("satsys must contain satellite-system letters such as 'G' or 'L'")
        data = {}
        for epoch, satellites in self._iter_epochs(toSI=toSI):
            data[epoch] = {
                sat: [entry["x"], entry["y"], entry["z"], entry["c"], entry["events"]]
                for sat, entry in satellites.items()
                if sat[0] in systems and "x" in entry
            }
        return data

    def get_satellite(self, satid: str, toSI: bool = True) -> dict:
        """Return all available P/V fields for one satellite, indexed by epoch.

        Unlike the old reader, missing epochs are omitted instead of appearing
        as empty dictionaries.  An unknown satellite raises a helpful KeyError.
        """
        satid = satid.upper()
        if self.sat_ids and satid not in self.sat_ids:
            raise KeyError(f"Satellite {satid!r} is absent from SP3 file")
        return {
            epoch: satellites[satid]
            for epoch, satellites in self._iter_epochs(toSI=toSI)
            if satid in satellites and "x" in satellites[satid]
        }
