#!/usr/bin/env python3
"""Fresh-process acceptance workload timings and sampled fork-tree RSS.

Run outside the filesystem/process sandbox on macOS. ps RSS is sampled, so peaks
are lower bounds; the aggregate double-counts pages shared by forked processes.
This is not proportional set size, physical allocation, or an OS memory limit.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time


def process_rows():
    result = subprocess.run(["ps", "-axo", "pid=,ppid=,rss="], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    rows = {}
    for line in result.stdout.splitlines():
        pid, parent, rss = (int(value) for value in line.split())
        rows[pid] = (parent, rss * 1024)  # macOS/Linux ps report KiB
    return rows


def descendants(rows, root):
    selected = {root}
    while True:
        found = {pid for pid, (parent, _) in rows.items() if parent in selected}
        combined = selected | found
        if combined == selected:
            return selected & rows.keys()
        selected = combined


def execute(command, environment, log_path, phase_path, interval, capture_timeline=False):
    peaks = dict(sampled_main_peak_rss_bytes=0, sampled_worker_peak_rss_bytes=0,
                 sampled_aggregate_tree_peak_rss_bytes=0, max_concurrent_fork_workers=0,
                 analysis_samples=0, process_peak_samples={})
    timeline = []
    first_sample = None
    with log_path.open("w") as log:
        process = subprocess.Popen(command, env=environment, stdout=log, stderr=log)
        while process.poll() is None:
            if phase_path.exists():
                state = json.loads(phase_path.read_text())
                if state["phase"] == "analysis":
                    rows = process_rows()
                    root = state["pid"]
                    tree = descendants(rows, root)
                    workers = tree - {root}
                    if capture_timeline:
                        if first_sample is None:
                            first_sample = time.monotonic()
                        timeline.append(dict(unix_seconds=time.time(), seconds_from_first_sample=time.monotonic()-first_sample,
                            main_pid=root, rss_bytes_by_pid={str(pid): rows[pid][1] for pid in tree}, fork_workers=len(workers)))
                    peaks["analysis_samples"] += 1
                    main_rss = rows.get(root, (0, 0))[1]
                    peaks["sampled_main_peak_rss_bytes"] = max(peaks["sampled_main_peak_rss_bytes"], main_rss)
                    peaks["sampled_aggregate_tree_peak_rss_bytes"] = max(
                        peaks["sampled_aggregate_tree_peak_rss_bytes"], sum(rows[pid][1] for pid in tree))
                    peaks["max_concurrent_fork_workers"] = max(peaks["max_concurrent_fork_workers"], len(workers))
                    for pid in tree:
                        key = str(pid)
                        peaks["process_peak_samples"][key] = max(peaks["process_peak_samples"].get(key, 0), rows[pid][1])
                    if workers:
                        peaks["sampled_worker_peak_rss_bytes"] = max(
                            peaks["sampled_worker_peak_rss_bytes"], max(rows[pid][1] for pid in workers))
            time.sleep(interval)
        if process.returncode:
            raise RuntimeError(f"R workload failed ({process.returncode}): {log_path}")
    if not peaks["analysis_samples"]:
        # Tiny compact analyses can finish between samples. Their resource/timing
        # values must remain diagnostic, not silently count as a measured peak.
        peaks["memory_observation"] = "analysis finished between RSS samples"
    else:
        peaks["memory_observation"] = "sampled lower bound; shared pages overcounted in aggregate"
    if capture_timeline:
        peaks["analysis_rss_timeline"] = timeline
    return peaks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library", type=Path)
    parser.add_argument("data", type=Path)
    parser.add_argument("fixtures", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--build-provenance", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--cases", nargs="+")
    parser.add_argument("--cache", nargs="+", choices=["cold", "warm"], default=["cold", "warm"])
    parser.add_argument("--poll-interval", type=float, default=0.1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--rss-timeline", action="store_true")
    args = parser.parse_args()
    if args.repetitions < 3:
        parser.error("at least3 measured repetitions plus warm-up are required")
    if platform.system() not in ("Darwin", "Linux"):
        parser.error("process sampler currently supports macOS and Linux")
    if args.poll_interval < 0.02 or args.poll_interval > 1:
        parser.error("poll interval must be between0.02 and1second")
    build = json.loads(args.build_provenance.read_text())
    if not build.get("verified_clean_release_build") or build["library"] != str(args.library.resolve()) or build["source_commit"] != args.source_commit:
        parser.error("verified release build provenance does not match library/source commit")
    installed_dll = args.library / "sceptre" / "libs" / "sceptre.so"
    if hashlib.sha256(installed_dll.read_bytes()).hexdigest() != build["native_dll_sha256"]:
        parser.error("installed native library differs from verified build")
    args.output.mkdir(parents=True, exist_ok=True)
    worker = Path(__file__).with_name("resource-worker.R").resolve()
    environment = os.environ.copy()
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        environment[variable] = "1"
    manifest = json.loads((args.data / "manifest.json").read_text())
    cases = [case["name"] for case in manifest["cases"]]
    if args.cases:
        if set(args.cases) - set(cases):
            parser.error("unknown case")
        cases = args.cases
    plan = dict(label=args.label, source_commit=args.source_commit,
                library=str(args.library.resolve()), R_LIBS=environment.get("R_LIBS"),
                cases=cases, cache=args.cache, workers=args.workers,
                repetitions=args.repetitions, poll_interval=args.poll_interval, rss_timeline=args.rss_timeline,
                worker_sha256=hashlib.sha256(worker.read_bytes()).hexdigest(),
                harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                manifest=manifest, build_provenance=build)
    plan_path = args.output / "plan.json"
    if args.resume and plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        raise RuntimeError("resume plan differs from recorded source, input, or execution configuration")
    plan_path.write_text(json.dumps(plan, indent=2) + "\n")
    samples = []
    for case in cases:
        for cache in args.cache:
            for workers in args.workers:
                for repetition in range(args.repetitions + 1):
                    stem = args.output / f"{case}-{cache}-w{workers}-r{repetition}"
                    result_path = stem.with_suffix(".json")
                    memory_path = stem.with_suffix(".memory.json")
                    phase_path = stem.with_suffix(".phase.json")
                    if not (args.resume and result_path.exists() and memory_path.exists()):
                        phase_path.unlink(missing_ok=True)
                        command = ["Rscript", "--vanilla", str(worker), str(args.library.resolve()),
                                   str((args.data / f"{case}.rds").resolve()), cache, str(workers),
                                   str(result_path.resolve()), str(phase_path.resolve()), str(args.fixtures.resolve())]
                        started = time.monotonic()
                        memory = execute(command, environment, stem.with_suffix(".log"), phase_path, args.poll_interval, args.rss_timeline)
                        memory["fresh_process_wall_seconds"] = time.monotonic() - started
                        memory_path.write_text(json.dumps(memory, indent=2) + "\n")
                    result = json.loads(result_path.read_text())
                    result.update(json.loads(memory_path.read_text()))
                    result.update(repetition=repetition, warmup=repetition == 0)
                    samples.append(result)
                    print(f"{args.label} {case} {cache} workers={workers} rep={repetition}: "
                          f"{result['elapsed_seconds']:.3f}s ({result['analysis_samples']} RSS samples)", flush=True)
    summary = []
    for case in cases:
        for cache in args.cache:
            for workers in args.workers:
                selected = [row for row in samples if row["case"] == case and row["cache"] == cache
                            and row["workers"] == workers and not row["warmup"]]
                row = dict(case=case, cache=cache, workers=workers,
                           median_elapsed_seconds=statistics.median(s["elapsed_seconds"] for s in selected),
                           input_sha256=selected[0]["input_sha256"],
                           scientific_sha256=sorted({s["scientific_sha256"] for s in selected}),
                           B=selected[0]["B"], passing_pairs=selected[0]["passing_pairs"],
                           min_analysis_samples=min(s["analysis_samples"] for s in selected),
                           max_concurrent_fork_workers=max(s["max_concurrent_fork_workers"] for s in selected))
                for metric in ("sampled_main_peak_rss_bytes", "sampled_worker_peak_rss_bytes", "sampled_aggregate_tree_peak_rss_bytes"):
                    row[metric] = max(s[metric] for s in selected)
                row["timing_use"] = "diagnostic only" if row["median_elapsed_seconds"] < 1 else "review against10percent median elapsed budget"
                summary.append(row)
    if hashlib.sha256(installed_dll.read_bytes()).hexdigest() != build["native_dll_sha256"]:
        raise RuntimeError("native library changed during resource trials")
    report = dict(schema_version=1, label=args.label, source_commit=args.source_commit,
                  platform=platform.platform(), poll_interval_seconds=args.poll_interval,
                  memory_scope="analysis-only sampled main/fork descendant RSS; aggregate overcounts shared pages; sampled peaks are lower bounds",
                  dense_budget_scope="declared256MiB response-assay materialization budget; not an enforced OS memory limit and not a claim that512MiB cannot fit the Mac",
                  data_manifest=manifest, build_provenance=build, samples=samples, summary=summary)
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
