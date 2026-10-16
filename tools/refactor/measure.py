#!/usr/bin/env python3
"""Reproducible static engineering metrics; no package code is executed.

Provide cyclocomp 1.1.2 via R_LIBS and lizard 1.24.1 via PYTHONPATH.
Usage: python3 measure.py SOURCE_ROOT OUTPUT_DIR [--source-commit LABEL]
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys

import lizard


def digest(data):
    return hashlib.sha256(data).hexdigest()


def quantile(values, q):
    """R type=7 / linearly interpolated sample quantile."""
    values = sorted(values)
    if not values:
        return None
    position = (len(values) - 1) * q
    left = int(position)
    fraction = position - left
    return values[left] + fraction * (values[min(left + 1, len(values) - 1)] - values[left])


def distribution(rows, field):
    values = [row[field] for row in rows]
    return dict(count=len(values), maximum=max(values) if values else None,
                p90=quantile(values, .90), p95=quantile(values, .95))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--source-commit", default=None,
                        help="Provenance label supplied by caller; hashes, not this label, verify content")
    args = parser.parse_args()
    root = args.source_root.resolve(strict=True)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    version = importlib.metadata.version("lizard")
    if version != "1.24.1":
        raise SystemExit("Expected lizard 1.24.1, found " + version)
    rscript = Path(__file__).with_name("measure-r.R")
    subprocess.run(["Rscript", "--vanilla", str(rscript), str(root), str(output)], check=True)
    rmetrics = json.loads((output / "r-metrics.json").read_text())
    # 'nd' is Lizard's bundled nesting-depth extension, not a custom brace count.
    analyzer = lizard.FileAnalyzer(lizard.get_extensions(["nd"]))
    probe = analyzer.analyze_source_code("analyzer-probe.cpp", "int f(int x){if(x){while(x){x--;}}return x;}").function_list[0]
    if probe.cyclomatic_complexity != 3 or probe.max_nesting_depth != 2:
        raise SystemExit("Lizard complexity/nesting analyzer contract changed")
    cpp_files, cpp_functions = [], []
    for path in sorted((root / "src").rglob("*")):
        if path.suffix.lower() not in {".c", ".cc", ".cpp", ".cxx", ".h", ".hpp", ".hxx"} or path.name == "RcppExports.cpp":
            continue
        source = path.read_text()
        relative = path.relative_to(root).as_posix()
        result = analyzer.analyze_source_code(relative, source)
        cpp_files.append(dict(path=relative, sha256=digest(path.read_bytes()),
                              physical_lines=len(source.splitlines()), lizard_nloc=result.nloc))
        lines = source.splitlines(keepends=True)
        for function in result.function_list:
            cpp_functions.append(dict(id=f"{relative}::{function.name}@{function.start_line}",
                name=function.name, signature=function.long_name, file=relative,
                start_line=function.start_line, end_line=function.end_line,
                source_line_span_sha256=digest("".join(lines[function.start_line-1:function.end_line]).encode()),
                cyclomatic_complexity=function.cyclomatic_complexity,
                max_nesting_depth=function.max_nesting_depth, nloc=function.nloc,
                token_count=function.token_count, parameter_count=function.parameter_count))
    cppmetrics = dict(schema_version=1, language="C++", analyzer=dict(name="lizard", version=version,
        extensions=["nd"], python=platform.python_version()), rules=dict(
        scope="Handwritten C/C++ source and headers under src; generated RcppExports.cpp excluded; tests and build artifacts excluded.",
        complexity="Unmodified Lizard C-like cyclomatic complexity with bundled default preprocessing/comment/line/token/condition processors.",
        nesting="Bundled Lizard nd extension (same pinned 1.24.1 distribution), per-function max_nesting_depth. Its token rules include control keywords and first logical operator in a condition, with special else-if handling; this differs from the custom R AST depth. No default zero placeholder used.",
        LOC="Lizard nloc and physical lines are supporting diagnostics only.",
        limitations="Static heuristic parsing, not a compiler or cognitive complexity. Semantic duplication/coupling and human comprehension are unmeasured."),
        files=cpp_files, functions=cpp_functions)
    (output / "cpp-metrics.json").write_text(json.dumps(cppmetrics, indent=2) + "\n")
    top = [row for row in rmetrics["functions"] if row["top_level"]]
    named_top = [row for row in top if row["kind"] != "anonymous"]
    summaries = dict(R_exclusive_all_function_units=distribution(rmetrics["functions"], "cyclomatic_exclusive_body"),
        R_exclusive_named_top_level=distribution(named_top, "cyclomatic_exclusive_body"),
        R_inclusive_named_top_level=distribution(named_top, "cyclomatic_inclusive_body"),
        R_exclusive_control_nesting=distribution(rmetrics["functions"], "control_nesting_exclusive_body"),
        CPP_cyclomatic=distribution(cpp_functions, "cyclomatic_complexity"),
        CPP_nesting=distribution(cpp_functions, "max_nesting_depth"))
    file_hashes = {row["path"]: row["sha256"] for row in rmetrics["files"] + cpp_files}
    scripts = {path.name: digest(path.read_bytes()) for path in [Path(__file__), rscript]}
    selected_names = {"run_perm_test_in_memory", "run_crt_in_memory_v2", "run_sceptre_analysis_high_level",
        "discovery_ntcells_perm_test", "discovery_ntcells_crt", "run_low_level_test_full_v4",
        "compute_null_full_statistics", "fit_and_evaluate_skew_normal"}
    selected = [row for row in rmetrics["functions"] + cpp_functions if row["name"] in selected_names]
    source_hash = digest(json.dumps(file_hashes, sort_keys=True, separators=(",", ":")).encode())
    summary = dict(schema_version=1, source_commit_label=args.source_commit,
        source_label_authority="Caller-supplied; content identity established by per-file SHA256 and aggregate handwritten_source_sha256.",
        handwritten_source_sha256=source_hash, source_files_sha256=file_hashes,
        measurement_scripts_sha256=scripts, summary=summaries,
        quantile_rule="R type=7: linear interpolation at (n-1)*p, no complexity sum reported.",
        trouble_spots=selected,
        highest_R_inclusive=sorted(named_top, key=lambda row: (-row["cyclomatic_inclusive_body"], row["id"]))[:10],
        highest_R_exclusive=sorted(rmetrics["functions"], key=lambda row: (-row["cyclomatic_exclusive_body"], row["id"]))[:10],
        highest_CPP=sorted(cpp_functions, key=lambda row: (-row["cyclomatic_complexity"], row["id"]))[:10],
        LOC_diagnostic=dict(R_physical_lines=sum(row["physical_lines"] for row in rmetrics["files"]),
            R_noncomment_token_lines=sum(row["noncomment_token_lines"] for row in rmetrics["files"]),
            CPP_physical_lines=sum(row["physical_lines"] for row in cpp_files),
            CPP_lizard_nloc=sum(row["lizard_nloc"] for row in cpp_files)),
        unmeasured=["Automated cross-file call-graph coupling", "Token/block duplication", "Semantic rule ownership/locality",
            "Cognitive complexity", "Human comprehension"],
        interpretation="Engineering baseline; R and C++ scores use different analyzers and must be compared within language/rule/version. No claim of human comprehension or scientific validity.")
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(dict(source_sha256=source_hash, summary=summaries), indent=2))


if __name__ == "__main__":
    main()
