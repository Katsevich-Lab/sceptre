#!/usr/bin/env python3
"""Post-format only the new long Rcpp wrappers; R AST/argument order is unchanged."""
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
targets = {
    'run_association_test', 'sample_calibration_groups',
    'sample_calibration_pairs', 'evaluate_pair_eligibility_cpp',
    'sparse_pair_counts', 'sparse_group_counts',
}
lines = []
wrapped = 0
for line in path.read_text().splitlines():
    signature = re.fullmatch(r'(\w+) <- function\((.*)\) \{', line)
    call = re.fullmatch(r'(\s*)\.Call\((.*)\)', line)
    if len(line) > 80 and signature and signature[1] in targets:
        arguments = [part.strip() for part in signature[2].split(',')]
        assert all(re.fullmatch(r'[A-Za-z][A-Za-z0-9._]*', a) for a in arguments)
        lines.append(signature[1] + ' <- function(')
        lines.extend('    ' + a + (',' if i < len(arguments)-1 else '')
                     for i, a in enumerate(arguments))
        lines.append(') {')
        wrapped += 1
    elif len(line) > 80 and call:
        arguments = [part.strip() for part in call[2].split(',')]
        symbol = re.fullmatch(r'`_sceptre_(\w+)`', arguments[0])
        if symbol and symbol[1] in targets:
            assert all(re.fullmatch(r'[A-Za-z][A-Za-z0-9._]*', a)
                       for a in arguments[1:])
            lines.append(call[1] + '.Call(')
            lines.extend(call[1] + '    ' + a + (',' if i < len(arguments)-1 else '')
                         for i, a in enumerate(arguments))
            lines.append(call[1] + ')')
            wrapped += 1
        else:
            lines.append(line)
    else:
        lines.append(line)
path.write_text('\n'.join(lines) + '\n')
print(f'Wrapped {wrapped} new generated-binding lines without changing arguments')
