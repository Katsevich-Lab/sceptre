#!/usr/bin/env python3
"""Reviewed-candidate site inventory, not semantic architectural certification."""
import argparse
import hashlib
import json
from pathlib import Path
import re

parser = argparse.ArgumentParser()
parser.add_argument("source", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
patterns = {
    "private_ondisc_calls": r"ondisc:::\s*[A-Za-z_][A-Za-z_0-9]*",
    "raw_storage_access_candidates": r"@(h5_file|ptr|i|j|p|x)\b",
    "view_remapping_candidates": r"cells_in_use\[|all_nt_idxs\[|\[\s*cells_in_use\s*[,\]]|\[\s*all_nt_idxs\s*[,\]]|combined_idxs\s*<-|trt_idxs\s*<-\s*seq\(",
    "response_cache_reset_events": r"@response_precomputations\s*<-\s*list\(\)",
    "response_cache_writebacks": r"@response_precomputations\s*<-",
}
sites = {name: [] for name in patterns}
source_hashes = {}
for path in sorted((args.source / "R").glob("*.R")):
    if path.name == "RcppExports.R":
        continue
    relative = path.relative_to(args.source).as_posix()
    source_hashes[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    for line_number, text in enumerate(path.read_text().splitlines(), 1):
        if text.lstrip().startswith("#"):
            continue
        for name, pattern in patterns.items():
            if name == "response_cache_writebacks" and re.search(patterns["response_cache_reset_events"], text):
                continue
            if re.search(pattern, text):
                sites[name].append(dict(path=relative, line=line_number, text=text.strip()))
report = dict(schema_version=1, rules_version=1,
              scope="Non-generated production R; matching lines, not expression counts. Candidate sites require responsibility review.",
              source_hashes=source_hashes, patterns=patterns,
              counts={name: dict(lines=len(rows), files=len({row['path'] for row in rows})) for name, rows in sites.items()},
              sites=sites,
              unmeasured=["semantic coordinate-conversion count", "complete dynamic call graph", "token/block clone coverage"])
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report["counts"], indent=2))
