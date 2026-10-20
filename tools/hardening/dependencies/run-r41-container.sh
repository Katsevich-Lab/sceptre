#!/usr/bin/env bash
# Run from the repository root after building Dockerfile.r41's small context.
# Package sources are read-only; all compilation/library writes stay inside the
# disposable container. Only validation receipts are mounted writable.
set -euo pipefail
source_root=$(pwd)
evidence_root="$source_root/tools/hardening/dependencies/evidence/r41"
mkdir -p "$evidence_root"
docker run --rm --platform linux/arm64 --network none \
  -v "$source_root:/source:ro" -v "$evidence_root:/evidence" \
  -e SCEPTRE_FLOOR_LIBRARY=/opt/sceptre-floor-library \
  -e SCEPTRE_EXPECTED_R_MINOR=4.1 \
  -e SCEPTRE_DEPENDENCY_EVIDENCE=/evidence \
  "${SCEPTRE_R41_IMAGE:-sceptre-r41-dependency-floor:testthat-3.1.8}" bash -euo pipefail -c '
    mkdir -p /tmp/sceptre-floor-source
    tar -C /source --exclude=.git --exclude=.github --exclude=tools \
      --exclude="*.o" --exclude="*.so" --exclude="*.dll" \
      --exclude="*.tar.gz" --exclude="*.Rcheck" -cf - . |
      tar -C /tmp/sceptre-floor-source -xf -
    cp /opt/sceptre-floor/r-floor-dependencies.txt /evidence/base-image-dependency-session.txt
    cp /opt/sceptre-floor/r-floor-installed.tsv /evidence/base-image-installed.tsv
    Rscript --vanilla -e '\''stopifnot(all(vapply(setdiff(names(tools:::.split_dependencies(read.dcf("/source/DESCRIPTION")[1,"Imports"])),c("methods","stats")), requireNamespace, logical(1), quietly=TRUE)))'\''
    cd /tmp/sceptre-floor-source
    R CMD INSTALL --preclean --library=/opt/sceptre-floor-library .
    cd /source
    R_DEFAULT_PACKAGES=NULL Rscript --vanilla tools/hardening/dependencies/probe-declared-namespaces.R
    Rscript --vanilla tools/hardening/dependencies/probe-testthat-floor.R
    cd /tmp/sceptre-floor-source
    Rscript --vanilla /source/tools/hardening/dependencies/run-floor-suite.R /evidence
  '
