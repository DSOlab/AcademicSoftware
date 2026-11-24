"""Make a small, station-only SINEX teaching extract without covariance data.

Usage::

    python -m dsoclasses.sinex.slim ITRF2020-IGS-TRF.SNX.gz DYNG_ITRF2020.snx DYNG

The source may be compressed with gzip. The input is read sequentially, so
even a very large SINEX file does not need to be loaded into memory. This
creates a *reduced teaching extract*, not a replacement for the official
full-covariance reference-frame product.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


_BLOCKS = ('FILE/REFERENCE', 'SITE/ID', 'SOLUTION/EPOCHS', 'SOLUTION/ESTIMATE')
_PARAMETERS = {'STAX', 'STAY', 'STAZ', 'VELX', 'VELY', 'VELZ'}


def slim_sinex(source, destination, stations):
    """Retain only requested station XYZ/velocity records and needed blocks.

    Matrix, normal-equation, apriori, history, receiver/antenna and all other
    blocks are discarded. SOLUTION/ESTIMATE indices are renumbered and the
    header parameter count is updated. Existing destination files are never
    overwritten.
    """
    sites = {str(site).strip().upper() for site in stations}
    if not sites or any(len(site) != 4 or not site.isalnum() for site in sites):
        raise ValueError('stations must be four-character station codes')
    source = Path(source)
    destination = Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError('Source and destination must be different files')

    retained = {block: [] for block in _BLOCKS}
    seen = set()
    header = footer = None
    block = None
    opener = gzip.open if source.name.lower().endswith('.gz') else open
    with opener(source, 'rt', encoding='ascii', errors='replace') as stream:
        for line in stream:
            if line.startswith('%=SNX'):
                if header is not None:
                    raise ValueError('Multiple SINEX headers')
                header = line
            elif line.strip() == '%ENDSNX':
                footer = line
            elif line.startswith('+'):
                if block is not None:
                    raise ValueError(f'Nested SINEX block inside {block}')
                block = line[1:].split()[0]
                if block in retained:
                    if block in seen:
                        raise ValueError(f'Repeated SINEX block {block}')
                    seen.add(block)
                    retained[block].append(line)
            elif line.startswith('-'):
                name = line[1:].split()[0]
                if name != block:
                    raise ValueError(f'Unexpected block trailer {name!r}')
                if block in retained:
                    retained[block].append(line)
                block = None
            elif block in retained:
                if line.startswith('*') or block == 'FILE/REFERENCE':
                    retained[block].append(line)
                elif (block in ('SITE/ID', 'SOLUTION/EPOCHS')
                      and line[:5].strip().upper() in sites):
                    retained[block].append(line)
                elif (block == 'SOLUTION/ESTIMATE'
                      and line[14:18].strip().upper() in sites
                      and line[7:13].strip() in _PARAMETERS):
                    retained[block].append(line)

    if header is None or footer is None or block is not None:
        raise ValueError('Incomplete SINEX: expected header, closed blocks and footer')
    mandatory = {'SITE/ID', 'SOLUTION/EPOCHS', 'SOLUTION/ESTIMATE'}
    missing = mandatory - seen
    if missing:
        raise ValueError(f'SINEX lacks required blocks: {sorted(missing)}')
    counts = {site: 0 for site in sites}
    for row in retained['SOLUTION/ESTIMATE']:
        if row.startswith(' '):
            counts[row[14:18].strip().upper()] += 1
    for site, count in counts.items():
        if count < 6:
            raise LookupError(f'SINEX has fewer than six XYZ/velocity records for {site}')

    new_estimates = []
    n = 0
    for row in retained['SOLUTION/ESTIMATE']:
        if row.startswith(' '):
            n += 1
            new_estimates.append(f'{n:6d}' + row[6:])
        else:
            new_estimates.append(row)
    if n > 99999:
        raise ValueError('SINEX header cannot represent more than 99999 estimates')
    retained['SOLUTION/ESTIMATE'] = new_estimates
    header = header[:60] + f'{n:5d}' + header[65:]

    # All validation above is read-only; only now create the deliverable.
    with destination.open('x', encoding='ascii') as output:
        output.write(header)
        for name in _BLOCKS:
            output.writelines(retained[name])
        output.write(footer)
    return n


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', help='Original .SNX or .SNX.gz')
    parser.add_argument('destination', help='New, uncompressed .snx file')
    parser.add_argument('stations', nargs='+', help='Four-character station IDs')
    args = parser.parse_args()
    n = slim_sinex(args.source, args.destination, args.stations)
    print(f'Wrote {args.destination}: {n} position/velocity estimates')


if __name__ == '__main__':
    main()
