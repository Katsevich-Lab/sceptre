#!/usr/bin/env python3
"""Compare frozen-input resource runs; threshold failures require disposition."""
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("reference", type=Path)
parser.add_argument("candidate", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
reference = json.loads(args.reference.read_text())
candidate = json.loads(args.candidate.read_text())
key = lambda row: (row["case"], row["cache"], row["workers"])
original = {key(row): row for row in reference["summary"]}
updated = {key(row): row for row in candidate["summary"]}
if original.keys() != updated.keys():
    raise SystemExit("reference/candidate workloads differ")
rows = []
for setting in original:
    before, after = original[setting], updated[setting]
    row = dict(case=setting[0], cache=setting[1], workers=setting[2])
    for field in ("input_sha256", "B", "passing_pairs"):
        if before[field] != after[field]:
            raise SystemExit(f"{setting}: changed {field}")
    row["scientific_hash_exact"] = before["scientific_sha256"] == after["scientific_sha256"]
    row["within_run_scientific_determinism"] = len(before["scientific_sha256"]) == len(after["scientific_sha256"]) == 1
    row["reference_median_seconds"] = before["median_elapsed_seconds"]
    row["candidate_median_seconds"] = after["median_elapsed_seconds"]
    row["elapsed_ratio"] = after["median_elapsed_seconds"] / before["median_elapsed_seconds"] if before["median_elapsed_seconds"] else None
    row["elapsed_gate"] = ("diagnostic only: reference median under1second" if before["median_elapsed_seconds"] < 1
                           else "pass" if row["elapsed_ratio"] <= 1.1 else "review required")
    row["memory"] = []
    for metric in ("sampled_main_peak_rss_bytes", "sampled_worker_peak_rss_bytes", "sampled_aggregate_tree_peak_rss_bytes"):
        baseline, result = before[metric], after[metric]
        permitted_increase = max(0.1 * baseline, 64 * 1024**2)
        gate = ("not applicable: serial" if "worker_peak" in metric and setting[2] == 1
                else "missing samples" if not baseline or not result
                else "pass" if result - baseline <= permitted_increase else "review required")
        row["memory"].append(dict(metric=metric, reference_bytes=baseline,
                                   candidate_bytes=result, permitted_increase_bytes=permitted_increase,
                                   gate=gate))
    row["fork_observation"] = ("serial" if setting[2] == 1 else "observed2workers"
                              if min(before["max_concurrent_fork_workers"], after["max_concurrent_fork_workers"]) >= 2
                              else "incomplete fork sampling")
    rows.append(row)
report = dict(schema_version=1, reference_source_commit=reference["source_commit"],
              candidate_source_commit=candidate["source_commit"], comparison=rows,
              interpretation="Initial review thresholds, not proof of speedup. Sampled peaks are lower bounds; aggregate RSS counts shared pages repeatedly. Exact output hash differences require numerical/identity investigation, never automatic snapshot replacement.")
args.output.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(rows, indent=2))
