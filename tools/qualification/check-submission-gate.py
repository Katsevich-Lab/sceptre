#!/usr/bin/env python3
"""Black-box gate policy probes, not package/scientific qualification.
Args: real checked archive, real bound diagnostics-summary.json, NEW output.
The fixture summaries are deliberately synthetic; never use them for admission.
"""
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

archive, summary, out = map(Path, sys.argv[1:])
archive, summary, out = archive.resolve(), summary.resolve(), out.resolve()
out.mkdir(parents=True, exist_ok=False)
doc = json.loads(summary.read_text())
invocation = json.loads(summary.with_name('invocation.json').read_text())
assert invocation['archive_sha256'] == hashlib.sha256(archive.read_bytes()).hexdigest()
assert invocation['diagnostics_summary_sha256'] == hashlib.sha256(summary.read_bytes()).hexdigest()
if isinstance(doc['warning'], dict) and 'checkBiocDepsDESC' in doc['warning']:
    severity, key = 'warning', 'checkBiocDepsDESC'
else:
    severity, key = 'note', 'checkCodingPractice'
assert isinstance(doc[severity], dict) and key in doc[severity]
results = []
cases = [
    ('unchanged', 0, ''),
    ('warning_removed', 0, ''),
    ('new_message', 1, 'New diagnostic within inherited category'),
    ('duplicate_message', 1, 'New diagnostic within inherited category'),
    ('wrong_archive', 1, 'Summary/archive identity mismatch'),
    ('stale_summary', 1, 'Diagnostic summary bytes do not match invocation'),
    ('weakened_options', 1, 'Wrong checker options'),
]
for name, expected, reason in cases:
    directory = out / name
    directory.mkdir()
    fixture, binding = copy.deepcopy(doc), copy.deepcopy(invocation)
    if name == 'warning_removed':
        fixture['warning'] = []
        fixture['counts']['warnings'] = 0
    elif name == 'new_message':
        fixture[severity][key].append({key: 'A synthetic newly introduced diagnostic'})
    elif name == 'duplicate_message':
        fixture[severity][key].append(copy.deepcopy(fixture[severity][key][0]))
    elif name == 'wrong_archive':
        binding['archive_sha256'] = '0' * 64
    elif name == 'stale_summary':
        fixture['counts']['notes'] += 1
    elif name == 'weakened_options':
        binding['flags'].append('no-check-coding-practice')
    text = json.dumps(fixture)
    (directory / 'diagnostics-summary.json').write_text(text)
    if name != 'stale_summary':
        binding['diagnostics_summary_sha256'] = hashlib.sha256(text.encode()).hexdigest()
    (directory / 'invocation.json').write_text(json.dumps(binding))
    command = [sys.executable, str(Path(__file__).with_name('submission-gate.py')),
               str(archive), str(directory / 'diagnostics-summary.json'), str(directory / 'gate')]
    process = subprocess.run(command, capture_output=True, text=True, env=os.environ)
    (directory / 'terminal.log').write_text(process.stdout + process.stderr)
    passed = (process.returncode == 0) == (expected == 0)
    if expected:
        passed = passed and reason in process.stderr
    results.append(dict(probe=name, exit_code=process.returncode,
        expected_exit=expected, passed=passed,
        fixture_kind='Synthetic policy fixture; not a package diagnostic result',
        note='No warning remained to remove' if name == 'warning_removed' and not doc['warning'] else ''))
(out / 'probes.json').write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps(results, indent=2))
raise SystemExit(0 if all(row['passed'] for row in results) else 1)
