#!/usr/bin/env python3
"""Explicit small-screen timing harness; never run as a package test."""
import argparse
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("library", type=Path)
parser.add_argument("fixtures", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--repetitions", type=int, default=3)
args = parser.parse_args()
if args.repetitions < 3:
    parser.error("At least three measured repetitions are required")
args.output.mkdir(parents=True, exist_ok=True)
worker = Path(__file__).with_name("benchmark-worker.R").resolve()
environment = os.environ.copy()
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    environment[variable] = "1"
samples = []
for case in ("low_complement_perm", "high_complement_crt", "backed_perm"):
    for cache in ("cold", "warm"):
        for repetition in range(args.repetitions + 1):
            stem = args.output / f"{case}-{cache}-{repetition}"
            result_path = stem.with_suffix(".json")
            command = ["Rscript", "--vanilla", str(worker), str(args.library.resolve()),
                       str(args.fixtures.resolve()), case, cache, str(result_path.resolve())]
            if platform.system() == "Darwin":
                command = ["/usr/bin/time", "-l", *command]
            elif platform.system() == "Linux":
                command = ["/usr/bin/time", "-v", *command]
            else:
                raise SystemExit("RSS collection currently supports macOS and Linux")
            with stem.with_suffix(".log").open("w") as log:
                run = subprocess.run(command, env=environment, stdout=log, stderr=log, check=False)
            if run.returncode:
                raise SystemExit(f"Benchmark failed: {stem.with_suffix('.log')}")
            result = json.loads(result_path.read_text())
            log_text = stem.with_suffix(".log").read_text()
            if platform.system() == "Darwin":
                rss = re.search(r"(\d+)\s+maximum resident set size", log_text)
                result["process_peak_rss_bytes"] = int(rss.group(1)) if rss else None
            else:
                rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", log_text)
                result["process_peak_rss_bytes"] = int(rss.group(1)) * 1024 if rss else None
            if result["process_peak_rss_bytes"] is None:
                raise SystemExit(f"RSS missing: {stem.with_suffix('.log')}")
            result.update(repetition=repetition, warmup=repetition == 0)
            samples.append(result)
summary = []
for case in sorted({row["case"] for row in samples}):
    for cache in ("cold", "warm"):
        selected = [row for row in samples if row["case"] == case and row["cache"] == cache and not row["warmup"]]
        summary.append(dict(case=case, cache=cache,
                            median_elapsed_seconds=statistics.median(row["elapsed_seconds"] for row in selected),
                            max_process_peak_rss_bytes=max(row["process_peak_rss_bytes"] for row in selected)))
report = dict(schema_version=1, scope="small synthetic characterization; diagnostic timing, not realistic performance acceptance",
              outstanding=["difficult-fit real sparse screen", "backed workload exceeding declared dense-memory budget", "fork-worker/process-tree RSS"],
              samples=samples, summary=summary)
(args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(summary, indent=2))
