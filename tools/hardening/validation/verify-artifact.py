"""Verify exact packaged implementation/test files against the checkout."""
import hashlib
import json
import pathlib
import subprocess
import sys
import tarfile

archive = pathlib.Path(sys.argv[1]).resolve()
output = pathlib.Path(sys.argv[2])
paths = subprocess.check_output(
    ["git", "ls-files", "-z", "R", "src", "man", "tests"]
).decode().rstrip("\0").split("\0")
rows = []
with tarfile.open(archive) as package:
    members = {member.name: member for member in package.getmembers()}
    for path in paths:
        name = "sceptre/" + path
        if name not in members:
            raise RuntimeError("Packaged source missing: " + path)
        packaged = package.extractfile(members[name]).read()
        local = pathlib.Path(path).read_bytes()
        if packaged != local:
            raise RuntimeError("Packaged source differs: " + path)
        rows.append({"path": path, "sha256": hashlib.sha256(local).hexdigest()})
result = {
    "archive": str(archive),
    "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
    "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip(),
    "all_equal": True,
    "file_count": len(rows),
    "files": rows,
}
output.write_text(json.dumps(result, indent=2) + "\n")
print(f"Exact archive/check-out equality: {len(rows)} R/native/help/test files")
