#!/usr/bin/env python3
"""Compare frozen-input resource runs; threshold failures require disposition."""
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("reference", type=Path)
parser.add_argument("candidate", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--reference-environment", type=Path)
parser.add_argument("--candidate-environment", type=Path)
args = parser.parse_args()
if bool(args.reference_environment) != bool(args.candidate_environment):
    parser.error("supply both dependency environment sidecars")
dependency_provenance = None
if args.reference_environment:
    original_environment = json.loads(args.reference_environment.read_text())
    updated_environment = json.loads(args.candidate_environment.read_text())
    original_dependencies = {name: value for name, value in original_environment["packages"].items() if name != "sceptre"}
    updated_dependencies = {name: value for name, value in updated_environment["packages"].items() if name != "sceptre"}
    dependency_provenance = dict(exact=original_dependencies == updated_dependencies,
        packages=sorted(original_dependencies),
        reference=str(args.reference_environment), candidate=str(args.candidate_environment))
reference = json.loads(args.reference.read_text())
candidate = json.loads(args.candidate.read_text())
if not reference.get("build_provenance") or not candidate.get("build_provenance"):
    raise SystemExit("Both reports require verified release build provenance")
build_configuration_exact = all(reference["build_provenance"][field] == candidate["build_provenance"][field]
    for field in ("compiler", "sdk", "native_compile_prefix"))
if not build_configuration_exact:
    raise SystemExit("Compiler, SDK, or native compilation flags differ")
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
    before_samples = [sample for sample in reference["samples"] if key(sample) == setting and not sample["warmup"]]
    after_samples = [sample for sample in candidate["samples"] if key(sample) == setting and not sample["warmup"]]
    row["warning_lists_exact"] = all(
        {json.dumps(sample[field], sort_keys=True) for sample in before_samples} ==
        {json.dumps(sample[field], sort_keys=True) for sample in after_samples}
        for field in ("preparation_warnings", "warm_warnings", "analysis_warnings"))
    row["within_run_configuration_consistency"] = all(
        all(sample[field] == samples[0][field] for sample in samples)
        for samples in (before_samples, after_samples)
        for field in ("input_sha256", "B", "passing_pairs", "threads", "R", "ondisc_version"))
    row["environment_exact"] = all(before_samples[0][field] == after_samples[0][field]
        for field in ("threads", "R", "ondisc_version", "external", "platform"))
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
              candidate_source_commit=candidate["source_commit"], dependency_provenance=dependency_provenance, build_configuration_exact=build_configuration_exact, comparison=rows,
              interpretation="Initial review thresholds, not proof of speedup. Sampled peaks are lower bounds; aggregate RSS counts shared pages repeatedly. Exact output hash differences require numerical/identity investigation, never automatic snapshot replacement.")
args.output.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(rows, indent=2))
