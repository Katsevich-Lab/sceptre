from pathlib import Path
import datetime
import hashlib
import json
import shutil

root = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
load = lambda path: json.loads((root / path).read_text())
prep = load("preparation.json")
controlled = load("qualification-comparison-arrival-controlled.json")
parallel = load("qualification-comparison-parallel-arrival.json")
serialized = load("qualification-comparison-serialized.json")
arrival = load("arrival-order-diagnostic.json")
assert controlled["scoped_qualification_passed"]
assert arrival["controlled_arrival_plans_exact"]
assert all(value is True for value in arrival["raw_worker_payloads_compared_by_existing_pod_id"].values())
assert not parallel["scoped_qualification_passed"]
assert not serialized["scoped_qualification_passed"]

immutability = {}
for label in ["reference", "candidate", "ondisc"]:
    package = "ondisc" if label == "ondisc" else "sceptre"
    package_dir = Path(prep[label]["library"]) / package
    differences = [name for name, item in prep[label]["installed_artifacts"].items()
                   if sha(package_dir / name) != item["sha256"]]
    assert not differences, (label, differences)
    immutability[label] = {"library": prep[label]["library"], "all_installed_artifacts_unchanged": True,
                           "artifact_count": len(prep[label]["installed_artifacts"])}
helper_readback = load("provenance/helper-consumer-readback.json")
assert sha(root / "provenance/helper-execution-start.R") == prep["helper"]["sha256"]
assert all(helper_readback["consumer_used_function_definitions_exact"].values())
assert sha(prep["helper"]["path"]) == sha(root / "provenance/helper-postrun.R")
assert sha(prep["fixture_manifest"]["path"]) == prep["fixture_manifest"]["sha256"]
for runner in prep["runners"]:
    assert sha(runner["staged_path"]) == runner["sha256"]

summary = {
    "schema_version": 1,
    "observed_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "reference_source": prep["reference"]["source_revision"],
    "candidate_source": prep["candidate"]["source_revision"],
    "build_provenance": "provenance/{reference,candidate}-release-build.json and build-readback.json",
    "build_mode": "Both verified clean preclean -O2 builds; matching compiler/SDK/native optimization flags.",
    "fixed_thread_environment": prep["thread_environment"],
    "installed_artifacts_unchanged": immutability,
    "runner_helper_fixture_hashes": "preparation.json",
    "helper_concurrent_change": {
        "execution_start_sha256": prep["helper"]["sha256"],
        "postrun_sha256": sha(root / "provenance/helper-postrun.R"),
        "both_actual_script_runs_recorded_start_hash": True,
        "consumed_import_prepare_definitions_ast_exact": helper_readback["consumer_used_function_definitions_exact"],
        "explanation": helper_readback["description"],
        "snapshots_diff_readback": "provenance/helper-*"},
    "scientific_float_tolerance": 1e-8,
    "actual_pipeline_scripts": {"reference_passed": 30, "candidate_passed": 30,
                                **controlled["comparisons"]["actual_scripts"]},
    "arrival_controlled_local_nextflow": {
        "status": "qualified",
        "runtime": "Nextflow25.10.4/Java17, array=null, containers/conda disabled, offline cached plugin, R seed271828, queueSize1, maxForks1, fair=true",
        "shared_actual_arrival_plan_exact": arrival["controlled_arrival_plans_exact"],
        "response_pod_maps_exact": controlled["execution"]["exact_response_to_pod_plan"],
        **controlled["comparisons"]["complete_nextflow"]},
    "parallel_two_task_nextflow": {
        "execution": "Both source versions completed standard and nuclear workflows successfully.",
        "ordered_gather_equality": "failed; retained",
        "cause": arrival["cause"],
        "all_five_raw_worker_pod_payloads_exact": True,
        "report": "qualification-comparison-parallel-arrival.json"},
    "serialization_only_nextflow": {
        "execution": "Both source versions completed standard and nuclear workflows successfully.",
        "shared_arrival_plan_exact": arrival["serialized_arrival_plans_exact"],
        "ordered_gather_equality": "failed; serialization alone did not enforce emissions in input order",
        "report": "qualification-comparison-serialized.json"},
    "uncached_book": {
        "reference_passed": 14, "candidate_passed": 14,
        "known_failure_count_each": 1,
        "known_failure": "import-data.qmd:28: installed sceptredata has two 20-byte, zero-line compressed barcode files; original and candidate both reject them.",
        "known_baseline_failure_preserved": controlled["execution"]["known_baseline_book_data_failure_preserved"],
        "whole_book_success_claimed": False,
        **controlled["comparisons"]["completed_uncached_book_cells"]},
    "no_published_output_sorted_or_modified": True,
    "qualified_scopes_complete": True,
    "remaining_gates": ["Default parallel collection order/factor-level sensitivity", "Unmodified default array configuration",
                        "Historical containers", "nf-test/hosted schedulers", "Valid CellRanger book example barcodes",
                        "Whole-book/full repeated-analysis chapters and external rd7 data"],
    "scope_revision_note": "This qualifies immutable production28df390 only. Later source-contract/readability edits require their own relevant equivalence evidence."
}
(root / "qualification-comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
prep["state"] = "completed_with_explicit_external_limitations"
prep["completed_at"] = summary["observed_at"]
prep["installed_artifact_immutability"] = immutability
(root / "preparation.json").write_text(json.dumps(prep, indent=2) + "\n")

durable = root / "durable"
durable.mkdir(exist_ok=False)
def copy(source, target):
    destination = durable / target
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(root / source, destination)

for name in ["preparation.json", "qualification-comparison.json", "qualification-comparison-parallel-arrival.json",
             "qualification-comparison-serialized.json", "qualification-comparison-arrival-controlled.json",
             "arrival-order-diagnostic.json", "serialized-control-plan.json", "arrival-control-plan.json",
             "compare-qualification.R", "compare-qualification-serialized.R", "compare-qualification-arrival-controlled.R",
             "diagnose-arrival.R", "finish-staging.py"]:
    copy(name, name)
for directory in ["runners", "provenance"]:
    shutil.copytree(root / directory, durable / directory)
for side in ["reference", "candidate"]:
    for lane, report_name, payload in [
        ("stages", "report.json", "comparison.rds"),
        ("nextflow", "nextflow-report.json", "nextflow-comparison.rds"),
        ("nextflow-serialized", "nextflow-report.json", "nextflow-comparison.rds"),
        ("nextflow-arrival-controlled", "nextflow-report.json", "nextflow-comparison.rds"),
        ("book", "book-report.json", "book-comparison.rds")]:
        copy(f"{side}/{lane}/{report_name}", f"reports/{side}-{lane}.json")
        copy(f"{side}/{lane}/{payload}", f"payloads/{side}-{lane}.rds")
        if lane.startswith("nextflow"):
            for role in ["standard", "nuclear"]:
                for log in ["stdout.log", "stderr.log", ".nextflow.log"]:
                    copy(f"{side}/{lane}/{role}/{log}", f"logs/{side}/{lane}/{role}/{log}")
                info = root / side / lane / role / "published/pipeline_info"
                for trace in info.glob("execution_trace_*.txt"):
                    copy(str(trace.relative_to(root)), f"logs/{side}/{lane}/{role}/{trace.name}")
    for stage in (root / side / "stages").iterdir():
        if stage.is_dir():
            for name in ["stdout.log", "stderr.log"]:
                if (stage / name).exists():
                    copy(str((stage / name).relative_to(root)), f"logs/{side}/stages/{stage.name}/{name}")
    for name in ["stdout.log", "stderr.log"]:
        copy(f"{side}/book/{name}", f"logs/{side}/book/{name}")
for log in root.glob("*-driver.log"):
    copy(log.name, f"logs/{log.name}")

readme = """# Final consumer qualification: immutable21f9ea0 versus28df390

Read qualification-comparison.json first. Actual30-script payloads and the
explicit fair/serialized local Nextflow control are exactly equal, with exact
identifiers, decisions, schema, ordering and actual pod/arrival plans. Floating
comparisons also pass1e-8. Source, installed-library and native-build provenance
is pinned and verified; seven thread environment variables were fixed1.

Both original two-task Nextflow workflows execute successfully, but their
ordered gathered output comparison fails because completion order changes
factor levels and tied-row ordering. Serialization alone does not fix emission
order. Both failed comparisons are retained; all five raw worker-pod payloads
are exactly equal. The qualified control adds fair=true before the original
collect, without editing pipeline sources or sorting published outputs.

Book runs each pass14/15 unchanged uncached cells. Both preserve the identical
known CellRanger example-data failure; no whole-book success is claimed.

Reports retain absolute paths of original executions. Portable comparison
payloads, source runners/inventory, build manifests/logs and task trace/logs are
copied here for durable review. Large raw work directories/plots/source archives
remain in the parent staging directory. artifact-manifest.json covers every
file here except itself; its external SHA256 is reported to the parent agent.

This evidence is specific to production28df390, not later source edits.
"""
(durable / "README.md").write_text(readme)
files = sorted(path for path in durable.rglob("*") if path.is_file())
manifest = {"schema_version": 1, "observed_at": summary["observed_at"],
            "excludes_self": True, "files": [{"path": str(path.relative_to(durable)),
              "bytes": path.stat().st_size, "sha256": sha(path)} for path in files]}
(durable / "artifact-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
shutil.copyfile(durable / "artifact-manifest.json", root / "artifact-manifest.json")
print(json.dumps({"durable": str(durable), "file_count": len(files),
                  "bytes": sum(path.stat().st_size for path in files),
                  "manifest_sha256": sha(durable / "artifact-manifest.json"),
                  "qualification": summary["qualified_scopes_complete"]}))
