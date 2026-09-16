"""Extract the 24 original files from the CPC ADGJ card-image archive."""

import argparse
from pathlib import Path
import re


def split(source, output):
    marker = re.compile(
        rb"^-- (This is|End of) .*?\b(PSM_(?:EE|EO|OE|OO)|"
        rb"(?:EE|EO|OE|OO)_(?:DATA|PARA|OUT[123]))\s*,"
    )
    files = {}
    name = None
    for line in source.read_bytes().splitlines(keepends=True):
        match = marker.match(line)
        if match:
            action, label = match.groups()
            label = label.decode('ascii')
            if action == b'This is':
                if name is not None or label in files:
                    raise ValueError(f'Duplicate or nested section: {label}')
                name = label
                files[name] = bytearray()
            else:
                if name != label:
                    raise ValueError(f'Mismatched end: {label}, expected {name}')
                name = None
        elif name is not None:
            files[name].extend(line)
    expected = {f'PSM_{kind}' for kind in ('EE', 'EO', 'OE', 'OO')}
    expected |= {f'{kind}_{suffix}' for kind in ('EE', 'EO', 'OE', 'OO')
                 for suffix in ('DATA', 'PARA', 'OUT1', 'OUT2', 'OUT3')}
    if name is not None or set(files) != expected:
        raise ValueError('Unterminated section or unexpected file inventory')
    for label, content in files.items():
        kind = label[-2:] if label.startswith('PSM_') else label[:2]
        filename = label + '.f' if label.startswith('PSM_') else label
        target = output / kind / filename
        if target.exists() and target.read_bytes() != content:
            raise FileExistsError(f'Refusing to overwrite changed file: {target}')
    for label, content in files.items():
        kind = label[-2:] if label.startswith('PSM_') else label[:2]
        filename = label + '.f' if label.startswith('PSM_') else label
        target = output / kind / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
    print(f'Extracted {len(files)} files to {output}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', nargs='?', type=Path,
                        default=Path(__file__).with_name('adgj'))
    parser.add_argument('output', nargs='?', type=Path,
                        default=Path(__file__).with_name('split'))
    args = parser.parse_args()
    split(args.source, args.output)
