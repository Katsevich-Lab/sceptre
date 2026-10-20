#!/usr/bin/env bash
set -euo pipefail
proof_dir=$(cd "$(dirname "$0")" && pwd)
source_dir=${1:-$(cd "$proof_dir/../../../.." && pwd)}
proof_tmp=$(mktemp -d "${TMPDIR:-/tmp}/sceptre-format-proof.XXXXXX")
trap 'rm -rf "$proof_tmp"' EXIT
verify_source() {
  python3 - "$proof_dir/source-manifest.json" "$source_dir" <<'PY'
from pathlib import Path
import hashlib, json, sys
manifest=json.loads(Path(sys.argv[1]).read_text())
root=Path(sys.argv[2])
changed=[name for name,sha in manifest.items()
         if not (root/name).is_file() or hashlib.sha256((root/name).read_bytes()).hexdigest()!=sha]
assert not changed, 'Source differs from captured qualification: '+str(changed)
print(f'PASS all {len(manifest)} captured R/man/vignette/DESCRIPTION/NAMESPACE source hashes')
PY
}
verify_source
tar -xzf "$proof_dir/snapshots.tar.gz" -C "$proof_tmp"
Rscript "$proof_dir/prove.R" "$proof_tmp/before" "$source_dir"
Rscript "$proof_dir/prove-docs.R" "$proof_tmp/before" "$source_dir"
Rscript "$proof_dir/prove-signals.R" "$proof_dir/string-wrap-manifest.json"
cmp "$proof_tmp/ast-equivalence.csv" "$proof_dir/ast-equivalence.csv"
cmp "$proof_tmp/roxygen-equivalence.csv" "$proof_dir/roxygen-equivalence.csv"
cp "$proof_tmp/before/R/RcppExports.R" "$proof_tmp/RcppExports.R"
python3 "$proof_dir/wrap-generated-bindings.py" "$proof_tmp/RcppExports.R"
cmp "$proof_tmp/RcppExports.R" "$source_dir/R/RcppExports.R"
verify_source
