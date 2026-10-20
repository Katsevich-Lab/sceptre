#!/usr/bin/env python3
"""Conservative original-submission diagnostic category AND instance gate.
Run from repo root with original BiocCheck1.49.30 first in R_LIBS.
Args: candidate archive, its matched diagnostics-summary.json, NEW output dir.
"""
import collections, csv, hashlib, json, pathlib, re, subprocess, sys, tarfile
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
invocation_path = summary.with_name("invocation.json")
invocation = json.loads(invocation_path.read_text())
archive_hash = hashlib.sha256(archive.read_bytes()).hexdigest()
assert invocation["archive_sha256"] == archive_hash, "Summary/archive identity mismatch"
assert invocation["diagnostics_summary_sha256"] == hashlib.sha256(summary.read_bytes()).hexdigest(), \
    "Diagnostic summary bytes do not match invocation"
assert invocation["checker"] == "1.49.30", "Wrong invocation checker"
assert str(invocation["Bioconductor"]) == "3.24", "Wrong invocation Bioconductor"
assert sorted(invocation["flags"]) == ["no-check-R-ver", "no-check-bioc-help"], \
    "Wrong checker options"
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
            rows.append(dict(file=file, call=call, line=int(line), column=int(col),
                text=(source / file).read_text().splitlines()[int(line)-1].strip()))
    return rows
old_display = display(old, out / "accepted/sceptre")
new_display = display(new, candidate)
key = lambda row: (row["file"], row["call"], row["text"])
prior = collections.Counter(map(key, old_display))
current = collections.Counter(map(key, new_display))
extra = current - prior
available = collections.defaultdict(collections.deque)
for row in old_display:
    available[key(row)].append(row)
with (accepted.parent / "instance-gate/display-call-lineage.csv").open() as f:
    reviewed = {(row["candidate_file"], int(row["candidate_line"])): row
                for row in csv.DictReader(f)}
lineage = []
for row in new_display:
    inherited = available[key(row)].popleft() if available[key(row)] else None
    original = reviewed[(inherited["file"], inherited["line"])] if inherited else {}
    lineage.append(dict(candidate_file=row["file"], candidate_line=row["line"],
        candidate_column=row["column"], call=row["call"], exact_source_line=row["text"],
        accepted_file=inherited["file"] if inherited else "",
        accepted_line=inherited["line"] if inherited else "",
        submitted_file=original.get("submitted_file", ""),
        submitted_line=original.get("submitted_line", ""),
        reviewed_lineage=original.get("lineage", ""),
        status="inherited_exact_line_with_reviewed_submission_lineage" if inherited else "NEW_DISPLAY_INSTANCE"))
if lineage:
    with (out / "display-call-lineage.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(lineage[0]))
        writer.writeheader()
        writer.writerows(lineage)
result = dict(passed=not extra,
    candidate_archive_sha256=archive_hash,
    diagnostics_summary_sha256=hashlib.sha256(summary.read_bytes()).hexdigest(),
    invocation_sha256=hashlib.sha256(invocation_path.read_bytes()).hexdigest(),
    display_calls=sum(current.values()),
    new_display_lines=[list(k) + [v] for k, v in extra.items()],
    matched_checker="1.49.30", flags=["no-check-R-ver", "no-check-bioc-help"],
    authority="21f9ea0 original submission via reviewed 01bb39c cleanup lineage",
    limitation="Exact flagged line/file multiplicity; semantic relocation needs explicit review. Not default-environment or scientific qualification.")
(out / "submission-gate.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
assert result["passed"], "New or modified display-call instance"
