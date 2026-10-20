#!/usr/bin/env python3
"""Run current-source native probes under supported Clang sanitizer runtimes.

Use a fresh evidence directory. The host R/Rcpp runtime is not rebuilt, and
ASan leak detection is disabled for that host. A deliberately faulty test-only
bridge must fail in each lane, proving the sanitizer is active. No OOM is used.
The source directory and installed libraries are never modified.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("root", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--lanes", nargs="+", choices=("undefined", "address"),
                    default=["undefined", "address"])
args = parser.parse_args()
root = args.root.resolve(strict=True)
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
if any(output.iterdir()):
    parser.error("evidence directory must be fresh")
rscript = shutil.which("Rscript")
if not rscript:
    parser.error("Rscript is required")
source = sorted((root / "src").glob("*")) + sorted((root / "tools/checks").glob("*"))
def hashes():
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in source if p.is_file()}
before = hashes()
summary = {"source_sha256": before, "platform": platform.platform(), "lanes": {},
           "scope": "Current source probes only; host R uninstrumented; ASan leaks disabled; no OS OOM"}
def run(command, env, log):
    with log.open("w") as stream:
        process = subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT)
    return process.returncode
for lane in args.lanes:
    directory = output / lane
    directory.mkdir()
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
               MAKEFLAGS="-j1", PKG_CXXFLAGS="-O1 -g -fno-omit-frame-pointer -fsanitize=" + lane,
               PKG_LIBS="-fsanitize=" + lane)
    if lane == "undefined":
        env["PKG_CXXFLAGS"] += " -fno-sanitize-recover=undefined"
        env["UBSAN_OPTIONS"] = "halt_on_error=1:print_stacktrace=1"
    else:
        env["ASAN_OPTIONS"] = "detect_leaks=0:halt_on_error=1"
        if sys.platform == "darwin":
            resource = subprocess.check_output(["clang++", "--print-resource-dir"], text=True).strip()
            runtime = Path(resource) / "lib/darwin/libclang_rt.asan_osx_dynamic.dylib"
            if not runtime.exists():
                summary["lanes"][lane] = {"status": "unsupported", "reason": str(runtime) + " missing"}
                continue
            env["DYLD_INSERT_LIBRARIES"] = str(runtime)
    def r_call(script, *arguments):
        if lane == "address" and sys.platform == "darwin":
            # Rscript launches R through a shell on this host; macOS strips
            # DYLD variables at that boundary. Inject into the Mach-O R binary
            # itself, before libR/Rcpp are loaded, rather than loading ASan late.
            home = subprocess.check_output(["R", "RHOME"], text=True).strip()
            env["R_HOME"] = home
            return [str(Path(home) / "bin/exec/R"), "--no-echo", "--no-restore",
                    "--file=" + str(script), "--args", *arguments]
        return [rscript, str(script), *arguments]
    # This negative control is separate from all production probes.
    bridge = directory / "negative.cpp"
    if lane == "undefined":
        body = "volatile int n = 2147483647; return n + 1;"
    else:
        body = "volatile int* p = new int[1]; delete[] p; return p[0];"
    bridge.write_text("#include <Rcpp.h>\n// [[Rcpp::export]]\nint sanitizer_negative() { " + body + " }\n")
    runner = directory / "negative.R"
    runner.write_text("Rcpp::sourceCpp(commandArgs(TRUE)[[1L]], rebuild=TRUE)\nsanitizer_negative()\n")
    negative = run(r_call(runner, str(bridge)), env, directory / "negative.log")
    text = (directory / "negative.log").read_text(errors="replace")
    detected = negative != 0 and (("runtime error: signed integer overflow" in text) if lane == "undefined"
                                 else ("ERROR: AddressSanitizer: heap-use-after-free" in text))
    record = {"negative_exit": negative, "negative_detected": detected, "checks": {}}
    summary["lanes"][lane] = record
    if not detected:
        record["status"] = "unsupported_or_failed_activation"
        continue
    for name in ("native-components", "sampling-ownership", "count-ownership",
                 "dispersion-exceptions", "result-ownership"):
        command = r_call(root / ("tools/checks/test-" + name + ".R"), str(root))
        if name != "native-components":
            command.append(str(directory / name))
        code = run(command, env, directory / (name + ".log"))
        record["checks"][name] = code
        if code:
            break
    record["status"] = "passed" if len(record["checks"]) == 5 and not any(record["checks"].values()) else "failed"
summary["source_unchanged"] = before == hashes()
(output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({"source_unchanged": summary["source_unchanged"],
                  "lanes": summary["lanes"]}, indent=2))
sys.exit(0 if summary["source_unchanged"] and all(v["status"] == "passed"
         for v in summary["lanes"].values()) else 1)
