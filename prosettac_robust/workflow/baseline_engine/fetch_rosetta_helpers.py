#!/usr/bin/env python3
"""Obtain pinned Rosetta helpers directly from their official source.

Rosetta's own license applies. This script neither grants a license nor
redistributes PyRosetta. Existing files with different content are never replaced.
"""
import argparse
import hashlib
import json
import os
import tempfile
import urllib.request
from pathlib import Path


def fetch(root, source_directory=None):
    vendor = root / 'vendor'
    metadata = json.loads((vendor / 'SOURCE.json').read_text())
    for relative, expected in metadata['files'].items():
        if not (vendor / relative).resolve().is_relative_to(vendor.resolve()):
            raise ValueError('Unsafe helper path')
        targets = [vendor / relative, root / 'runtime/adapter/vendor' / relative]
        for target in targets:
            if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest() != expected:
                raise ValueError('Existing helper differs from pinned source: ' + str(target))
        if all(target.exists() for target in targets):
            continue
        if source_directory:
            data = (source_directory / relative).read_bytes()
        else:
            url = ('https://raw.githubusercontent.com/RosettaCommons/rosetta/'
                   + metadata['commit'] + '/' + metadata['path'] + relative)
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read()
        if hashlib.sha256(data).hexdigest() != expected:
            raise ValueError('Downloaded helper hash mismatch: ' + relative)
        for target in targets:
            if target.exists():
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            fd, temporary = tempfile.mkstemp(dir=target.parent, prefix='.helper-')
            try:
                with os.fdopen(fd, 'wb') as out:
                    out.write(data)
                os.replace(temporary, target)
            finally:
                if os.path.exists(temporary):
                    os.unlink(temporary)
    print('All pinned Rosetta helpers verified; PyRosetta must be installed separately.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-directory', type=Path,
                        help='Local licensed checkout of source/scripts/python/public')
    args = parser.parse_args()
    fetch(Path(__file__).resolve().parent, args.source_directory)
