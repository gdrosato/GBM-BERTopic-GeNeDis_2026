#!/usr/bin/env python3
from pathlib import Path
import hashlib

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

failures = 0
for line in Path('SHA256SUMS').read_text(encoding='utf-8').splitlines():
    if not line.strip():
        continue
    expected, rel = line.split('  ', 1)
    p = Path(rel)
    if not p.exists():
        print('MISSING ', rel)
        failures += 1
        continue
    actual = digest(p)
    if actual != expected:
        print('FAILED  ', rel)
        failures += 1
    else:
        print('OK      ', rel)
raise SystemExit(1 if failures else 0)
