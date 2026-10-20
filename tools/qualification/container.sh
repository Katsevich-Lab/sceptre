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
  --entrypoint python3 "$image" /qualification/qualify.py \
  --archive /input/package.tar.gz --output /evidence/run --ondisc "$backend" \
  > "$output/process.log" 2>&1
