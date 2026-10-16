#!/usr/bin/env python3
"""Bounded, paired B6/T11 aggregation experiment. macOS; no fork workers."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import signal
import statistics
import subprocess
import time


def sha(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def headroom():
    raw = subprocess.check_output(["/usr/bin/vm_stat"], text=True)
    page = int(re.search(r"page size of (\d+) bytes", raw).group(1))
    counts = {name: int(re.search(rf"{re.escape(name)}:\s+(\d+)\.", raw).group(1))
              for name in ("Pages free", "Pages inactive", "Pages speculative")}
    return dict(vm_stat=raw, page_size=page, counts=counts,
                available_proxy_bytes=sum(counts.values()) * page,
                meaning="free+inactive+speculative; not an OS guarantee or hard native-memory limit")


def trial(worker, library, fixture, stem, mode, iterations, environment):
    ready = headroom()
    if ready["available_proxy_bytes"] < 8 * 1024**3:
        raise RuntimeError("Less than declared 8GiB headroom proxy; do not launch child.")
    result = stem.with_suffix(".json")
    phase = stem.with_suffix(".phase.json")
    clock = stem.with_suffix(".time.txt")
    command = ["/usr/bin/time", "-l", "-o", str(clock), "Rscript", "--vanilla", str(worker),
               str(library), str(fixture), str(result), str(phase), mode, str(iterations)]
    samples = []
    safety_stop = False
    started = time.monotonic()
    with stem.with_suffix(".log").open("w") as log:
        process = subprocess.Popen(command, env=environment, stdout=log, stderr=log)
        while process.poll() is None:
            if phase.exists():
                state = json.loads(phase.read_text())
                output = subprocess.run(["ps", "-o", "rss=", "-p", str(state["pid"])],
                                        text=True, capture_output=True)
                if output.returncode == 0 and output.stdout.strip():
                    rss = int(output.stdout.strip()) * 1024
                    samples.append(dict(unix_seconds=time.time(), rss_bytes=rss,
                                        phase=state["phase"], iteration=state["iteration"]))
                    if rss >= 2 * 1024**3:
                        safety_stop = True
                        os.kill(state["pid"], signal.SIGTERM)
                    elif state["phase"] in ("pre_call", "aggregation_end"):
                        Path(str(phase)+f".{state['phase']}-{state['iteration']}.ack").touch()
            time.sleep(.02)
        returncode = process.wait()
    memory = dict(headroom_before=ready, samples=samples, safety_stop=safety_stop,
                  returncode=returncode, fresh_process_wall_seconds=time.monotonic()-started)
    memory["os_time_output"] = clock.read_text() if clock.exists() else ""
    found = re.search(r"(\d+)\s+maximum resident set size", memory["os_time_output"])
    memory["os_max_rss_bytes"] = int(found.group(1)) if found else None
    active = [row["rss_bytes"] for row in samples if row["phase"] == "aggregation_start"]
    memory["sampled_active_peak_rss_bytes"] = max(active, default=None)
    before = [row["rss_bytes"] for row in samples if row["phase"] == "pre_call"]
    memory["sampled_pre_call_rss_bytes"] = before[0] if before else None
    after = [row["rss_bytes"] for row in samples if row["phase"] == "aggregation_end"]
    memory["sampled_post_call_rss_bytes"] = after[0] if after else None
    stem.with_suffix(".memory.json").write_text(json.dumps(memory, indent=2) + "\n")
    if returncode or safety_stop or not found:
        raise RuntimeError(f"Child failed/stopped or OS maximum RSS missing: {stem}; preserve artifacts.")
    record = json.loads(result.read_text())
    record["memory"] = memory
    return record


def summarize(rows):
    measured = [row for row in rows if not row["warmup"]]
    return dict(median_elapsed_seconds=statistics.median(row["elapsed_seconds"] for row in measured),
                max_os_peak_rss_bytes=max(row["memory"]["os_max_rss_bytes"] for row in measured),
                median_os_peak_rss_bytes=statistics.median(row["memory"]["os_max_rss_bytes"] for row in measured),
                max_sampled_active_peak_rss_bytes=max((row["memory"]["sampled_active_peak_rss_bytes"] or 0) for row in measured),
                measured_trials=len(measured))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("reference_library", "candidate_library", "reference_build", "candidate_build", "data", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    if platform.system() != "Darwin":
        parser.error("This predeclared protocol uses macOS time -l RSS byte units.")
    args.output.mkdir(parents=True, exist_ok=False)
    libraries = {"reference": args.reference_library.resolve(), "candidate": args.candidate_library.resolve()}
    builds = {"reference": json.loads(args.reference_build.read_text()), "candidate": json.loads(args.candidate_build.read_text())}
    for label in libraries:
        build = builds[label]
        if not build.get("verified_clean_release_build") or build["library"] != str(libraries[label]):
            raise RuntimeError(f"Unverified/mismatched build: {label}")
        if sha(libraries[label]/"sceptre/libs/sceptre.so") != build["native_dll_sha256"]:
            raise RuntimeError(f"Installed DLL differs from recorded immutable build: {label}")
    manifest = json.loads((args.data/"manifest.json").read_text())
    worker = Path(__file__).with_name("aggregation-worker.R").resolve()
    environment = os.environ.copy()
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        environment[variable] = "1"
    environment["R_MAX_VSIZE"] = "1536M"
    plan = dict(protocol="B6/T11; 1 warmup+3 measured fresh-process pairs; optional predeclared20-call timing batches",
                builds=builds, manifest=manifest, threads={key:environment[key] for key in
                ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
                R_LIBS=environment.get("R_LIBS"), R_MAX_VSIZE=environment["R_MAX_VSIZE"],
                harness_sha256=sha(Path(__file__)), worker_sha256=sha(worker),
                RSS_scope="macOS time-l whole child maximum bytes; 20ms sampled R-PID active peak may miss brief peaks; no forks",
                gates=dict(elapsed_increase=.10, RSS_increase_fraction=.10, RSS_increase_min_bytes=64*1024**2))
    (args.output/"plan.json").write_text(json.dumps(plan, indent=2)+"\n")
    results = []
    conclusions = []
    for case in manifest["cases"]:
        name = case["name"]
        fixture = (args.data/case["file"]).resolve()
        if sha(fixture) != case["sha256"]:
            raise RuntimeError(f"Fixture differs from predeclared manifest: {name}")
        correctness = {}
        for label, library in libraries.items():
            stem = args.output/f"{name}-{label}-correctness"
            correctness[label] = trial(worker, library, fixture, stem, "correctness", 1, environment)
        paths = [args.output/f"{name}-{label}-correctness.json.output.rds" for label in libraries]
        if sha(paths[0]) != sha(paths[1]):
            raise RuntimeError(f"Exact complete-output serialization mismatch: {name}; retain both outputs and investigate.")
        if correctness["reference"]["environment"] != correctness["candidate"]["environment"]:
            raise RuntimeError(f"Dependency/numerical/thread environment differs: {name}")
        rows = {label: [] for label in libraries}
        for repetition in range(4):
            order = list(libraries) if repetition % 2 == 0 else list(reversed(libraries))
            for label in order:
                stem = args.output/f"{name}-{label}-single-r{repetition}"
                row = trial(worker, libraries[label], fixture, stem, "timing", 1, environment)
                row.update(label=label, case=name, repetition=repetition, warmup=repetition==0)
                rows[label].append(row); results.append(row)
                print(f"{name} {label} single r{repetition}: {row['elapsed_seconds']:.3f}s", flush=True)
        single = {label:summarize(values) for label,values in rows.items()}
        timing = single
        timing_source = "single call"
        if min(row["median_elapsed_seconds"] for row in single.values()) < 1:
            batches = {label: [] for label in libraries}
            for repetition in range(4):
                order = list(libraries) if repetition % 2 == 0 else list(reversed(libraries))
                for label in order:
                    stem = args.output/f"{name}-{label}-batch20-r{repetition}"
                    row = trial(worker, libraries[label], fixture, stem, "batch", 20, environment)
                    row.update(label=label, case=name, repetition=repetition, warmup=repetition==0)
                    batches[label].append(row); results.append(row)
                    print(f"{name} {label} batch20 r{repetition}: {row['elapsed_seconds']:.3f}s", flush=True)
            timing = {label:summarize(values) for label,values in batches.items()}
            timing_source = "sum of20 independently loaded conversion intervals; single-call RSS remains primary"
        profile = {}
        for label, library in libraries.items():
            stem = args.output/f"{name}-{label}-profile"
            profile[label] = trial(worker, library, fixture, stem, "profile", 1, environment)
        elapsed_ratio = timing["candidate"]["median_elapsed_seconds"] / timing["reference"]["median_elapsed_seconds"]
        meaningful = min(row["median_elapsed_seconds"] for row in timing.values()) >= 1
        reference_peak = single["reference"]["max_os_peak_rss_bytes"]
        candidate_peak = single["candidate"]["max_os_peak_rss_bytes"]
        crossing = candidate_peak-reference_peak > max(.10*reference_peak,64*1024**2)
        active_ref = single["reference"]["max_sampled_active_peak_rss_bytes"]
        active_candidate = single["candidate"]["max_sampled_active_peak_rss_bytes"]
        active_crossing = active_candidate-active_ref > max(.10*active_ref,64*1024**2)
        conclusion = dict(case=name, exact_output_bytes=True, environment_equal=True,
            single=single, timing=timing, timing_source=timing_source,
            elapsed_ratio=elapsed_ratio, timing_gate="pass" if meaningful and elapsed_ratio<=1.10 else "fail" if meaningful else "inconclusive (<1s)",
            memory_gate="fail" if crossing or active_crossing else "pass", profiles=profile)
        conclusions.append(conclusion)
        (args.output/"summary.json").write_text(json.dumps(dict(plan=plan, trials=results, cases=conclusions),indent=2)+"\n")
        if crossing or active_crossing or (meaningful and elapsed_ratio>1.10):
            raise RuntimeError(f"Declared gate crossed: {name}; investigate and seek explicit review, do not hide original trials.")
    for label, library in libraries.items():
        if sha(library/"sceptre/libs/sceptre.so") != builds[label]["native_dll_sha256"]:
            raise RuntimeError(f"Library changed during trials: {label}")


if __name__ == "__main__":
    main()
