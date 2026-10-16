#!/usr/bin/env python3
"""Verify a fresh archive/install log before accepting a release benchmark."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
for name in ("repository", "source", "library", "install_log", "output"):
    parser.add_argument(name, type=Path)
parser.add_argument("--commit", required=True)
args = parser.parse_args()
commit = subprocess.check_output(["git", "rev-parse", args.commit], cwd=args.repository, text=True).strip()
tree = subprocess.check_output(["git", "ls-tree", "-r", commit, "--", "R", "src", "DESCRIPTION", "NAMESPACE"], cwd=args.repository, text=True)
hashes = {}
cpp = set()
for line in tree.splitlines():
    identity, name = line.split("\t", 1)
    expected = identity.split()[2]
    data = (args.source / name).read_bytes()
    found = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    if found != expected:
        raise SystemExit(f"source differs from{commit}:{name}")
    hashes[name] = found
    if name.startswith("src/") and name.endswith(".cpp"):
        cpp.add(Path(name).name)
log = args.install_log.read_text()
commands = [shlex.split(line) for line in log.splitlines() if " -c " in line]
compiled = {command[command.index("-c") + 1] for command in commands}
if compiled != cpp or "Nothing to be done" in log:
    raise SystemExit("install log does not prove fresh compilation of every native source")
prefixes = {tuple(command[:command.index("-c")]) for command in commands}
if len(prefixes) != 1:
    raise SystemExit("native compilation flags vary within build")
flags = list(next(iter(prefixes)))
if "-O2" not in flags or "-O0" in flags or "-UNDEBUG" in flags:
    raise SystemExit("benchmark requires the selected release-O2 compilation flags")
if "* DONE (sceptre)" not in log or str(args.library.resolve()) not in log:
    raise SystemExit("install did not complete into the declared frozen library")
dll = args.library / "sceptre" / "libs" / "sceptre.so"
report = dict(schema_version=1, source_commit=commit, library=str(args.library.resolve()),
    source=str(args.source.resolve()), source_file_git_blobs=hashes,
    compiler=next(line for line in log.splitlines() if line.startswith("using C++ compiler")),
    sdk=next((line for line in log.splitlines() if line.startswith("using SDK")), None),
    native_compile_prefix=flags, native_compiled_sources=sorted(compiled),
    verified_clean_release_build=True,
    native_dll_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(), native_dll_bytes=dll.stat().st_size,
    install_log_sha256=hashlib.sha256(args.install_log.read_bytes()).hexdigest(),
    note="All tracked production source files match the commit; every native source was compiled in this successful install with the recorded release flags. Source archives and libraries remain immutable during trials.")
args.output.write_text(json.dumps(report, indent=2) + "\n")
