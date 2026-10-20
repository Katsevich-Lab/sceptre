#!/usr/bin/env bash
# Usage: container.sh IMAGE PLATFORM SOURCE_ARCHIVE NEW_OUTPUT [present|absent|either]
# Requires a locally admitted image; never pulls or installs dependencies.
set -euo pipefail
image=$1
platform=$2
archive=$(cd "$(dirname "$3")" && pwd)/$(basename "$3")
output=$4
backend=${5:-either}
tools_root=$(cd "$(dirname "$0")" && pwd)
mkdir "$output"
output=$(cd "$output" && pwd)
docker image inspect "$image" > "$output/image.json"
# Container output goes into a new child directory. Keep image receipt outside.
docker run --rm --pull never --platform "$platform" --network none --cpus 2 \
  --memory 6g --pids-limit 512 --cap-drop ALL --security-opt no-new-privileges \
  --mount "type=bind,src=$tools_root,dst=/qualification,readonly" \
  --mount "type=bind,src=$archive,dst=/input/package.tar.gz,readonly" \
  --mount "type=bind,src=$output,dst=/evidence" \
  --env _R_CHECK_FORCE_SUGGESTS_=false \
  --env SCEPTRE_BACKEND_POLICY="$backend" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 --env MAKEFLAGS=-j2 \
  --entrypoint bash "$image" -euo pipefail -c '
    mkdir -p /evidence/run/library
    trap '\''code=$?; printf "{\"status\":\"%s\",\"exit_code\":%d}\\n" \
      "$([ "$code" = 0 ] && echo passed || echo failed)" "$code" > /evidence/run/receipt.json'\'' EXIT
    cp /input/package.tar.gz /evidence/run/package.tar.gz
    sha256sum /input/package.tar.gz > /evidence/run/archive-sha256.txt
    R CMD INSTALL --preclean --install-tests --library=/evidence/run/library \
      /input/package.tar.gz > /evidence/run/install.log 2>&1
    Rscript --vanilla /qualification/installed-suite.R /evidence/run/library \
      /evidence/run "$SCEPTRE_BACKEND_POLICY" > /evidence/run/tests.log 2>&1
  ' > "$output/process.log" 2>&1
