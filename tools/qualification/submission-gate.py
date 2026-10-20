#!/usr/bin/env python3
"""Conservative original-submission diagnostic category AND instance gate.
Run from repo root with original BiocCheck1.49.30 first in R_LIBS.
Args: candidate archive, its matched diagnostics-summary.json, NEW output dir.
"""
import collections, hashlib, json, pathlib, re, subprocess, sys, tarfile
archive, summary, out = map(pathlib.Path, sys.argv[1:])
archive, summary, out = archive.resolve(), summary.resolve(), out.resolve()
out.mkdir(parents=True, exist_ok=False)
root = pathlib.Path.cwd()
base = root / "tools/hardening/submission-baseline"
accepted = root / "tools/hardening/submission-cleanup/qualification/bioccheck-1.49.30"
for path, expected in [
    (base / "official/sceptre_0.99.0.tar.gz", "5300b3e3adfbb73a9a2173384ea3ec7e57b8cad8612f3fdc8116d5db5156eaf4"),
    (accepted / "sceptre_0.99.0.tar.gz", "606363fce32238354fc0d126d79f515973da7c8287ca43b14ce53ce2d00df4b7")]:
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, "Reference archive changed"
old = json.loads((accepted / "diagnostics-summary.json").read_text())
new = json.loads(summary.read_text())
assert new["metadata"]["BiocCheckVersion"] == "1.49.30", "Wrong checker"
assert str(new["metadata"]["BiocVersion"]) == "3.24", "Wrong Bioconductor target"
# No new category or raised severity. The accepted cleanup already maps to the
# actual submitted source; this is a transitive bound, never a count-only pass.
for severity in ["error", "warning", "note"]:
    assert set(new[severity]) <= set(old[severity]), (severity, "new diagnostic category")

def messages(obj, key):
    if isinstance(obj, dict):
        return [obj[key]] if key in obj else sum((messages(v, key) for v in obj.values()), [])
    if isinstance(obj, list):
        return sum((messages(v, key) for v in obj), [])
    return []
for severity, key in [("note", "checkCodingPractice"), ("warning", "checkBiocDepsDESC")]:
    assert messages(new[severity].get(key, []), key) == messages(old[severity].get(key, []), key), \
        "New diagnostic within inherited category"
for role, package in [("candidate", archive), ("submitted", base / "official/sceptre_0.99.0.tar.gz"),
                      ("accepted", accepted / "sceptre_0.99.0.tar.gz")]:
    with tarfile.open(package) as t:
        # Inputs are local pinned package artifacts. Reject traversal/links even
        # on Python3.9, where tarfile's newer extraction filter is unavailable.
        for member in t.getmembers():
            assert not pathlib.PurePosixPath(member.name).is_absolute()
            assert ".." not in pathlib.PurePosixPath(member.name).parts
            assert member.isfile() or member.isdir(), "Unsupported archive entry"
        t.extractall(out / role)
candidate = out / "candidate/sceptre"
subprocess.run(["Rscript", "--vanilla", "-e",
    'stopifnot(as.character(packageVersion("BiocCheck")) == "1.49.30")'], check=True)
subprocess.run(["Rscript", "--vanilla", str(base / "inventories/inventory-candidate.R"),
    str(candidate), str(out / "functions.csv")], check=True)
subprocess.run([sys.executable, str(base / "inventories/gate-candidate.py"),
    str(out / "submitted/sceptre"), str(candidate), str(out / "functions.csv"),
    str(base / "inventories/submitted-functions-1.49.30.csv"), str(out / "instances")], check=True)
# Display-call positions come from the actual checker. Preserve exact flagged
# source lines within each file, with multiplicity, against the accepted cleanup.
# A modified/moved display expression needs explicit reviewed lineage, not a
# silently enlarged allowlist. This gate is deliberately conservative.
def strings(obj):
    if isinstance(obj, str):
        return [obj]
    if isinstance(obj, dict):
        return sum((strings(v) for v in obj.values()), [])
    if isinstance(obj, list):
        return sum((strings(v) for v in obj), [])
    return []
def display(doc, source):
    rows = []
    for value in strings(doc["note"].get("checkCodingPractice", [])):
        m = re.fullmatch(r"(cat|print)\(\) in (R/.*?) \(line (\d+), column (\d+)\)", value)
        if m:
            call, file, line, col = m.groups()
            rows.append((file, call, (source / file).read_text().splitlines()[int(line)-1].strip()))
    return collections.Counter(rows)
prior = display(old, out / "accepted/sceptre")
current = display(new, candidate)
extra = current - prior
result = dict(passed=not extra, display_calls=sum(current.values()),
    new_display_lines=[list(k) + [v] for k, v in extra.items()],
    matched_checker="1.49.30", flags=["no-check-R-ver", "no-check-bioc-help"],
    authority="21f9ea0 original submission via reviewed 01bb39c cleanup lineage",
    limitation="Exact flagged line/file multiplicity; semantic relocation needs explicit review. Not default-environment or scientific qualification.")
(out / "submission-gate.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
assert result["passed"], "New or modified display-call instance"
