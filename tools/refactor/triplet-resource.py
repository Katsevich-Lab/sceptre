#!/usr/bin/env python3
"""Approved paired conversion protocols; one bounded child at a time.

Preparation and measurement require an explicitly confirmed quiet window. OS
maximum RSS is authoritative; 20 ms process samples are phase lower bounds.
"""
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
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def host_memory():
    raw = subprocess.check_output(["vm_stat"], text=True)
    page = int(re.search(r"page size of (\d+) bytes", raw)[1])
    pages = {name: int(count) for name, count in re.findall(r"^([^:\n]+):\s+(\d+)\.", raw, re.M)}
    # Inactive cache and speculative pages can be reclaimed. Purgeable pages
    # overlap these categories and are deliberately not added a second time.
    available = sum(pages.get(name, 0) for name in ["Pages free", "Pages inactive", "Pages speculative"]) * page
    return dict(vm_stat=raw, page_bytes=page, reclaimable_headroom_estimate_bytes=available,
                total_bytes=int(subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True)))


def process_rows():
    text = subprocess.check_output(["ps", "-axo", "pid=,ppid=,rss="], text=True)
    return {int(pid): (int(parent), int(rss) * 1024) for pid, parent, rss in
            (line.split() for line in text.splitlines())}


def execute(command, environment, stem, phase_path=None):
    stem.parent.mkdir(parents=True, exist_ok=True)
    log, oslog = Path(str(stem) + ".log"), Path(str(stem) + ".time.log")
    samples, stop_reason = [], None
    if phase_path is not None:
        phase_path.unlink(missing_ok=True)
    with log.open("w") as output:
        process = subprocess.Popen(["/usr/bin/time", "-l", "-o", str(oslog)] + command,
            env=environment, stdout=output, stderr=output, start_new_session=True)
        next_sample = time.monotonic()
        def abort_monitor(error):
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
            process.wait()
            return "Resource monitoring failed: " + str(error)
        while process.poll() is None:
            try:
                rows = process_rows()
            except Exception as error:
                stop_reason = abort_monitor(error)
                break
            state = None
            if phase_path is not None and phase_path.exists():
                state = json.loads(phase_path.read_text())
            if state:
                pid = state["pid"]
            else:
                children = [pid for pid, (parent, _) in rows.items() if parent == process.pid]
                pid = children[0] if children else process.pid
            rss = rows.get(pid, (0, 0))[1]
            samples.append(dict(unix_seconds=time.time(), monotonic_seconds=time.monotonic(),
                pid=pid, rss_bytes=rss, phase=state["phase"] if state else "startup"))
            if rss >= 2 * 1024**3:
                stop_reason = "2 GiB RSS safety stop; resource qualification incomplete"
                os.killpg(process.pid, signal.SIGTERM)
                process.wait()
                break
            next_sample += 0.02
            time.sleep(max(0, next_sample - time.monotonic()))
        returncode = process.wait()
    timing = oslog.read_text() if oslog.exists() else ""
    match = re.search(r"^\s*(\d+)\s+maximum resident set size", timing, re.M)
    os_peak = int(match[1]) if match else None
    # A short peak can escape the 20 ms sampler. Preserve the completed child,
    # then stop subsequent work if the authoritative whole-process peak crossed
    # the same safety bound; this outcome cannot qualify the experiment.
    if stop_reason is None and os_peak is not None and os_peak >= 2 * 1024**3:
        stop_reason = "OS peak reached 2 GiB RSS safety bound; resource qualification incomplete"
    result = dict(command=command, returncode=returncode, os_peak_rss_bytes=os_peak,
        os_peak_unit="Darwin time -l bytes", nominal_sampling_interval_seconds=0.02,
        sampled_peak_rss_bytes=max((s["rss_bytes"] for s in samples), default=0),
        samples=samples, safety_stop=stop_reason)
    if phase_path is not None:
        events_path = Path(str(stem) + ".json.events.jsonl")
        events = [json.loads(line) for line in events_path.read_text().splitlines()] if events_path.exists() else []
        result["events"] = events
        starts = [e for e in events if e["phase"] == "conversion_start"]
        ends = [e for e in events if e["phase"] == "conversion_end"]
        windows = list(zip(starts, ends))
        selected = [s["rss_bytes"] for s in samples if any(
            start["unix_seconds"] <= s["unix_seconds"] <= end["unix_seconds"] for start, end in windows)]
        result["sampled_conversion_peak_rss_bytes"] = max(selected) if selected else None
        before, after = [], []
        for start, end in windows:
            pre = [s for s in samples if s["unix_seconds"] <= start["unix_seconds"]]
            post = [s for s in samples if s["unix_seconds"] >= end["unix_seconds"]]
            before.append(pre[-1]["rss_bytes"] if pre else None)
            after.append(post[0]["rss_bytes"] if post else None)
        result["pre_conversion_rss_bytes"] = before
        result["post_conversion_rss_bytes"] = after
        result["phase_sampling_limit"] = "Samples can miss short peaks; OS whole-child maximum remains total authority."
    Path(str(stem) + ".memory.json").write_text(json.dumps(result, indent=2) + "\n")
    if returncode or stop_reason or match is None:
        raise RuntimeError(f"Stopped resource child: {stem}; code={returncode}, reason={stop_reason}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["prepare", "correctness", "measure", "summarize"])
    parser.add_argument("output", type=Path)
    parser.add_argument("--recipe", choices=["triplet", "dense"], default="triplet")
    parser.add_argument("--original-library", type=Path, required=True)
    parser.add_argument("--candidate-library", type=Path, required=True)
    parser.add_argument("--original-build", type=Path, required=True)
    parser.add_argument("--candidate-build", type=Path, required=True)
    parser.add_argument("--original-archive", type=Path, required=True)
    parser.add_argument("--candidate-archive", type=Path, required=True)
    parser.add_argument("--quiet-window-confirmed", action="store_true")
    args = parser.parse_args()
    if platform.system() != "Darwin":
        parser.error("This execution protocol currently uses Darwin time -l byte units.")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    tool = Path(__file__).resolve().parent
    protocol = tool / ("DENSE_CONVERSION_COST.md" if args.recipe == "dense" else "cleanup/triplet-resource-plan.md")
    generator = tool / ("dense-conversion-data.R" if args.recipe == "dense" else "triplet-resource-data.R")
    expected_settings = 4 if args.recipe == "dense" else 8
    env = os.environ.copy()
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
        env[name] = "1"
    env["R_MAX_VSIZE"] = "1536M"
    builds = {name: json.loads(path.read_text()) for name, path in
              [("original", args.original_build), ("candidate", args.candidate_build)]}
    libraries = {"original": args.original_library.resolve(), "candidate": args.candidate_library.resolve()}
    archives = {label: dict(path=str(path.resolve()), sha256=sha(path)) for label, path in
        [("original", args.original_archive), ("candidate", args.candidate_archive)]}
    for name in builds:
        if (not builds[name]["verified_clean_release_build"] or
            builds[name]["library"] != str(libraries[name]) or
            sha(libraries[name] / "sceptre/libs/sceptre.so") != builds[name]["native_dll_sha256"]):
            raise RuntimeError("Resource build/library verification failed: " + name)
    if builds["original"]["source_commit"] != "21f9ea098b69c444a884a49bc254de711968e3b5":
        raise RuntimeError("Resource original must be the pinned 21f9ea0 release build.")
    if builds["original"]["native_compile_prefix"] != builds["candidate"]["native_compile_prefix"]:
        raise RuntimeError("Original and candidate native release compiler flags differ.")
    if args.phase != "summarize":
        if not args.quiet_window_confirmed:
            parser.error("Wait for an explicitly confirmed quiet window before executing resource children.")
        memory = host_memory()
        (output / (args.phase + "-host-memory.json")).write_text(json.dumps(memory, indent=2) + "\n")
        if memory["reclaimable_headroom_estimate_bytes"] < 8 * 1024**3:
            raise RuntimeError("Less than 8 GiB headroom: resource qualification incomplete; no child launched.")
    clock = output / "triplet-clock.so"
    if args.phase == "prepare":
        if (output / "plan.json").exists():
            raise RuntimeError("Prepare requires a fresh output directory.")
        # Compile the same benchmark-only monotonic clock once for both builds.
        clock_source = output / "triplet-clock.c"
        clock_source.write_bytes((tool / "triplet-clock.c").read_bytes())
        with (output / "clock-build.log").open("w") as log:
            subprocess.run(["R", "CMD", "SHLIB", "-o", str(clock), str(clock_source)],
                           env=env, cwd=output, stdout=log, stderr=log, check=True)
        execute(["Rscript", "--vanilla", str(generator), str(output / "fixtures")],
                env, output / "fixture-preparation")
        plan = dict(recipe=args.recipe, protocol_sha256=sha(protocol),
            tool_sha256={p.name: sha(p) for p in [Path(__file__), generator,
                tool / "triplet-resource-worker.R", tool / "triplet-resource-compare.R", tool / "triplet-clock.c"]},
            clock_sha256=sha(clock), builds=builds, source_archives=archives, environment={name: env.get(name) for name in
                ["R_LIBS", "R_MAX_VSIZE", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]},
            fixture_manifest=json.loads((output / "fixtures/manifest.json").read_text()),
            cpu=subprocess.check_output(["sysctl", "-n", "machdep.cpu.brand_string"], text=True).strip(),
            available_memory=memory, platform=platform.platform(), pair_orders=["original-first", "candidate-first", "candidate-first", "original-first"])
        copies = output / "harness"
        copies.mkdir()
        for name in plan["tool_sha256"]:
            (copies / name).write_bytes((tool / name).read_bytes())
        (copies / protocol.name).write_bytes(protocol.read_bytes())
        plan["installed_package_sha256"] = {label: {
            str(path.relative_to(libraries[label])): sha(path)
            for path in sorted((libraries[label] / "sceptre").rglob("*")) if path.is_file()}
            for label in libraries}
        (output / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
        return
    plan = json.loads((output / "plan.json").read_text())
    if (plan["recipe"] != args.recipe or plan["protocol_sha256"] != sha(protocol) or
        plan["builds"] != builds or plan["source_archives"] != archives or sha(clock) != plan["clock_sha256"]):
        raise RuntimeError("Prepared source/library/clock provenance changed.")
    for name, expected in plan["tool_sha256"].items():
        if sha(tool / name) != expected:
            raise RuntimeError("Prepared resource tool changed: " + name)
    fixtures = plan["fixture_manifest"]["fixtures"]
    cases = [(name, orientation) for name in fixtures for orientation in ["row", "column"]]
    if len(cases) != expected_settings:
        raise RuntimeError("Fixture count differs from the declared recipe.")
    def child(label, case, orientation, mode, repetition, batch=1):
        if fixtures[case]["nnz"] == 8000000:
            scaling = []
            for memory_path in output.glob("*/1000000-*.memory.json"):
                observed = json.loads(memory_path.read_text())
                pre = next((value for value in observed.get("pre_conversion_rss_bytes", []) if value is not None), None)
                peak = observed.get("os_peak_rss_bytes")
                if pre is not None and peak is not None:
                    scaling.append(dict(trial=str(memory_path), pre_rss_bytes=pre,
                        whole_child_peak_bytes=peak,
                        conservative_8M_extrapolation_bytes=pre + 8 * max(0, peak - pre)))
            if not scaling:
                raise RuntimeError("Missing small-case scaling evidence before 8M conversion.")
            assessment = dict(measurement=False, purpose="conservative safety check before the large case",
                small_case_evidence=scaling, rss_stop_bytes=2 * 1024**3)
            (output / "large-case-scaling-check.json").write_text(json.dumps(assessment, indent=2) + "\n")
            if max(row["conservative_8M_extrapolation_bytes"] for row in scaling) >= 2 * 1024**3:
                raise RuntimeError("Small-case scaling threatens the RSS stop; large qualification remains incomplete.")
        stem = output / mode / f"{case}-{orientation}-{label}-r{repetition}-b{batch}"
        if Path(str(stem) + ".json").exists():
            raise RuntimeError("Refusing to overwrite an existing resource trial: " + str(stem))
        fixture = output / "fixtures" / fixtures[case]["file"]
        if sha(fixture) != fixtures[case]["sha256"]:
            raise RuntimeError("Immutable resource fixture changed.")
        phase = Path(str(stem) + ".phase.json")
        command = ["Rscript", "--vanilla", str(tool / "triplet-resource-worker.R"),
            str(libraries[label]), str(fixture), orientation, mode, str(stem) + ".json", str(phase), str(clock), str(batch)]
        if args.recipe == "dense":
            command.append(str(fixtures[case]["nnz"]))
        execute(command, env, stem, phase)
        print(stem.name, flush=True)
        return json.loads(Path(str(stem) + ".json").read_text())
    if args.phase == "correctness":
        checks = []
        for case, orientation in cases:
            results = {name: child(name, case, orientation, "correctness", 0) for name in builds}
            if not results["candidate"]["intervals"][0]["correctness"]["caller_unchanged"]:
                raise RuntimeError("Candidate mutated triplet input.")
            paths = [output / "correctness" / f"{case}-{orientation}-{name}-r0-b1.json.output.rds" for name in builds]
            comparison = output / "correctness" / f"{case}-{orientation}-comparison.json"
            execute(["Rscript", "--vanilla", str(tool / "triplet-resource-compare.R"),
                *map(str, paths), str(comparison)], env, comparison.with_suffix(""))
            checks.append(dict(case=case, orientation=orientation, results=results,
                               exact_output=json.loads(comparison.read_text())))
        (output / "correctness.json").write_text(json.dumps(checks, indent=2) + "\n")
    if args.phase == "measure":
        for case, orientation in cases:
            for label in builds:
                child(label, case, orientation, "performance", 0)
            for repetition, order in enumerate(plan["pair_orders"], 1):
                for label in (["original", "candidate"] if order == "original-first" else ["candidate", "original"]):
                    child(label, case, orientation, "performance", repetition)
            for label in builds:
                child(label, case, orientation, "profile", 0)
            if args.recipe == "triplet" and fixtures[case]["nnz"] == 8000000:
                medians = [statistics.median(json.loads((output / "performance" /
                    f"{case}-{orientation}-{name}-r{r}-b1.json").read_text())["elapsed_seconds"] for r in range(1, 5)) for name in builds]
                if min(medians) < 1:
                    for repetition, order in enumerate(plan["pair_orders"], 1):
                        for label in (["original", "candidate"] if order == "original-first" else ["candidate", "original"]):
                            child(label, case, orientation, "performance", repetition, batch=5)
    if args.phase in ["measure", "summarize"]:
        trials = []
        summaries = []
        for case, orientation in cases:
            for batch in [1, 5]:
                selected = {}
                for label in builds:
                    rows = []
                    for repetition in range(0 if batch == 1 else 1, 5):
                        stem = output / "performance" / f"{case}-{orientation}-{label}-r{repetition}-b{batch}"
                        path = Path(str(stem) + ".json")
                        if path.exists():
                            row = json.loads(path.read_text()); row["memory"] = json.loads(Path(str(stem) + ".memory.json").read_text())
                            row.update(case=case, label=label, repetition=repetition, warmup=repetition == 0, batch=batch)
                            trials.append(row)
                            if repetition: rows.append(row)
                    selected[label] = rows
                if not selected["original"] and not selected["candidate"]: continue
                if any(len(rows) != 4 for rows in selected.values()):
                    raise RuntimeError("Incomplete paired trials; no passing summary is possible.")
                summary = dict(case=case, orientation=orientation, batch=batch, builds={})
                for label, rows in selected.items():
                    times = [row["elapsed_seconds"] for row in rows]
                    os_peaks = [row["memory"]["os_peak_rss_bytes"] for row in rows]
                    phase_peaks = [row["memory"]["sampled_conversion_peak_rss_bytes"] for row in rows]
                    summary["builds"][label] = dict(median_elapsed_seconds=statistics.median(times), elapsed_range=[min(times), max(times)],
                        max_os_peak_rss_bytes=max(os_peaks), median_os_peak_rss_bytes=statistics.median(os_peaks),
                        max_sampled_conversion_peak_rss_bytes=max(p for p in phase_peaks if p is not None) if any(phase_peaks) else None,
                        median_sampled_conversion_peak_rss_bytes=statistics.median(p for p in phase_peaks if p is not None) if any(phase_peaks) else None)
                original, candidate = [summary["builds"][name] for name in builds]
                meaningful = original["median_elapsed_seconds"] >= 1 and candidate["median_elapsed_seconds"] >= 1
                summary["timing_gate"] = dict(meaningful=meaningful, candidate_fraction_increase=candidate["median_elapsed_seconds"] / original["median_elapsed_seconds"] - 1,
                    passed=(candidate["median_elapsed_seconds"] <= 1.1 * original["median_elapsed_seconds"]) if meaningful else None)
                summary["memory_gates"] = {}
                if batch == 1:
                    for metric in ["max_os_peak_rss_bytes", "max_sampled_conversion_peak_rss_bytes"]:
                        old, new = original[metric], candidate[metric]
                        summary["memory_gates"][metric] = dict(original=old, candidate=new,
                            delta_bytes=new-old if old is not None and new is not None else None,
                            delta_bytes_per_stored_entry=(new-old) / fixtures[case]["nnz"] if old is not None and new is not None else None,
                            delta_bytes_per_dense_input_element=(new-old) / fixtures[case]["elements"] if old is not None and new is not None and args.recipe == "dense" else None,
                            allowed_increase=max(0.1 * old, 64 * 1024**2) if old is not None else None,
                            passed=new <= old + max(0.1 * old, 64 * 1024**2) if old is not None and new is not None else None)
                summaries.append(summary)
        if args.recipe == "triplet":
            for case, orientation in cases:
                single = next(row for row in summaries if row["case"] == case and
                              row["orientation"] == orientation and row["batch"] == 1)
                short_large = fixtures[case]["nnz"] == 8000000 and any(
                    row["median_elapsed_seconds"] < 1 for row in single["builds"].values())
                if short_large and not any(row["case"] == case and row["orientation"] == orientation and
                                           row["batch"] == 5 for row in summaries):
                    raise RuntimeError("Required short-case timing batches are incomplete.")
        profiles = []
        for path in sorted((output / "profile").glob("*-r0-b1.json")):
            record = json.loads(path.read_text())
            stem = str(path)[:-5]
            memory = json.loads(Path(stem + ".memory.json").read_text())
            allocation_path = Path(str(path) + ".rprofmem")
            allocations = [int(match[1]) for line in allocation_path.read_text().splitlines()
                           if (match := re.match(r"^(\d+)\s*:", line))]
            phase_windows = []
            open_phases = {}
            for event in memory["events"]:
                phase = event["phase"]
                if phase.endswith("_start") and phase != "conversion_start":
                    open_phases.setdefault(phase[:-6], []).append(event)
                elif phase.endswith("_end") and phase != "conversion_end":
                    key = phase[:-4]
                    if open_phases.get(key):
                        start = open_phases[key].pop()
                        selected = [sample["rss_bytes"] for sample in memory["samples"] if
                                    start["unix_seconds"] <= sample["unix_seconds"] <= event["unix_seconds"]]
                        before = [sample["rss_bytes"] for sample in memory["samples"] if sample["unix_seconds"] <= start["unix_seconds"]]
                        phase_windows.append(dict(phase=key, start=start, end=event,
                            elapsed_seconds=event["monotonic_seconds"] - start["monotonic_seconds"],
                            sampled_peak_rss_bytes=max(selected) if selected else None,
                            pre_phase_rss_bytes=before[-1] if before else None))
            peak = memory["sampled_conversion_peak_rss_bytes"]
            pre = memory["pre_conversion_rss_bytes"][0]
            nnz = fixtures[Path(record["fixture"]).stem]["nnz"]
            profiles.append(dict(result=record, result_path=str(path), memory=memory,
                stored_entries=nnz,
                allocation_trace=str(allocation_path), allocation_trace_sha256=sha(allocation_path),
                recorded_R_allocation_bytes=sum(allocations), max_recorded_R_allocation_bytes=max(allocations, default=0),
                recorded_R_allocation_bytes_per_entry=sum(allocations) / nnz,
                dense_input_elements=fixtures[Path(record["fixture"]).stem].get("elements"),
                recorded_R_allocation_bytes_per_dense_input_element=sum(allocations) / fixtures[Path(record["fixture"]).stem]["elements"] if args.recipe == "dense" else None,
                net_conversion_rss_increment_bytes=peak-pre if peak is not None and pre is not None else None,
                net_conversion_rss_increment_bytes_per_entry=(peak-pre) / nnz if peak is not None and pre is not None else None,
                sorting_table_and_sort_windows=phase_windows,
                allocation_limit="Rprofmem totals are allocation traffic, not peak live workspace; native allocations may be absent."))
        if len(profiles) != 2 * expected_settings:
            raise RuntimeError("Incomplete allocation profiles; qualification remains incomplete.")
        counts = dict(warmups=sum(row["warmup"] for row in trials),
            measured_single_conversions=sum(not row["warmup"] and row["batch"] == 1 for row in trials),
            measured_batch_children=sum(row["batch"] == 5 for row in trials), profiles=len(profiles))
        correctness = json.loads((output / "correctness.json").read_text())
        if (counts["warmups"] != 2 * expected_settings or
            counts["measured_single_conversions"] != 8 * expected_settings or
            len(correctness) != expected_settings or
            (args.recipe == "dense" and counts["measured_batch_children"] != 0)):
            raise RuntimeError("Incomplete warmups, single-conversion trials or correctness cases; qualification remains incomplete.")
        crossings = [dict(case=row["case"], orientation=row["orientation"], batch=row["batch"],
            timing_gate=row["timing_gate"], memory_gates=row["memory_gates"]) for row in summaries if
            row["timing_gate"]["passed"] is False or any(gate["passed"] is False for gate in row["memory_gates"].values())]
        dense_costs = []
        if args.recipe == "dense":
            for case, orientation in cases:
                paired = {label: next(profile for profile in profiles if
                    Path(profile["result"]["fixture"]).stem == case and
                    profile["result"]["orientation"] == orientation and
                    Path(profile["result"]["package"]).parent.resolve() == libraries[label])
                    for label in builds}
                old, new = [paired[label] for label in builds]
                delta = new["recorded_R_allocation_bytes"] - old["recorded_R_allocation_bytes"]
                dense_costs.append(dict(case=case, orientation=orientation,
                    input_elements=fixtures[case]["elements"], input_payload_bytes=fixtures[case]["input_payload_bytes"],
                    original_recorded_R_allocation_bytes=old["recorded_R_allocation_bytes"],
                    candidate_recorded_R_allocation_bytes=new["recorded_R_allocation_bytes"],
                    delta_recorded_R_allocation_bytes=delta,
                    delta_recorded_R_allocation_bytes_per_input_element=delta / fixtures[case]["elements"],
                    original_max_recorded_R_allocation_bytes=old["max_recorded_R_allocation_bytes"],
                    candidate_max_recorded_R_allocation_bytes=new["max_recorded_R_allocation_bytes"],
                    limit="Profile allocation traffic is not peak live workspace; OS and sampled phase gates remain separately reported."))
        report = dict(protocol=plan, timing_threshold_seconds=1, summaries=summaries, all_timing_trials=trials,
            measurement_status="completed; gate crossings require investigation or an explicit reviewed exception" if crossings else
                "completed; interpret inconclusive short timing gates separately",
            gate_crossings=crossings, accepted_correctness_cost_exception=None,
            counts=counts, correctness=correctness,
            profiles=profiles, dense_cost_characterization=dense_costs,
            limits=["Phase RSS is sampled, not an exact allocator peak.",
                    "Conversion peak minus pre-conversion RSS is a net increment, not isolated sorting allocation.",
                    "Passing this focused gate does not waive integrated workflow/resource/consumer gates."])
        (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
        lines = ["# " + ("Dense general-conversion cost" if args.recipe == "dense" else "Triplet conversion resource measurements"), "",
            "Acceptance remains subject to the approved gate/exception review. Short timing gates are inconclusive.", "",
            "Original: `" + builds["original"]["source_commit"] + "`; candidate: `" + builds["candidate"]["source_commit"] + "`.",
            "Counts: " + ", ".join(f"{key}={value}" for key, value in counts.items()) + ".", "",
            "| Case | View | Batch | Original median s | Candidate median s | Time gate | OS peak delta MiB | OS gate | Sampled phase gate |",
            "| --- | --- | ---: | ---: | ---: | --- | ---: | --- | --- |"]
        for row in summaries:
            old, new = [row["builds"][label] for label in builds]
            timing = row["timing_gate"]["passed"]
            gates = row["memory_gates"]
            os_gate = gates.get("max_os_peak_rss_bytes", {})
            phase_gate = gates.get("max_sampled_conversion_peak_rss_bytes", {})
            delta = os_gate.get("delta_bytes")
            lines.append("| " + " | ".join(map(str, [row["case"], row["orientation"], row["batch"],
                round(old["median_elapsed_seconds"], 6), round(new["median_elapsed_seconds"], 6),
                "inconclusive" if timing is None else timing,
                "n/a" if delta is None else round(delta / 1024**2, 3),
                os_gate.get("passed", "batch excluded"), phase_gate.get("passed", "batch excluded")])) + " |")
        lines += ["", "Gate crossings: " + str(len(crossings)) + ". No correctness-cost exception is automatically granted.", "",
            "RSS is in bytes in JSON; sampled phase peaks can miss short peaks. Rprofmem records allocation traffic, not isolated live sorting workspace.",
            "The full JSON retains every warmup, trial, profile, input mutation/output equality check, build/archive hash and threshold comparison.",
            "Integrated workflow, backed, worker-memory/cache and consumer gates remain separate."]
        (output / "summary.md").write_text("\n".join(lines) + "\n")
    for name in builds:
        if sha(libraries[name] / "sceptre/libs/sceptre.so") != builds[name]["native_dll_sha256"]:
            raise RuntimeError("Installed native library changed during resource protocol.")
        for relative, expected in plan["installed_package_sha256"][name].items():
            if sha(libraries[name] / relative) != expected:
                raise RuntimeError("Installed package file changed during resource protocol: " + relative)
    for name, expected in plan["tool_sha256"].items():
        if sha(tool / name) != expected:
            raise RuntimeError("Resource tool changed during execution: " + name)


if __name__ == "__main__":
    main()
