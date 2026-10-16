#!/bin/sh
set -eu
characterization_root=/private/tmp/sceptre-import-aggregation-characterization-6a43fa5
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
Rscript --vanilla "$characterization_root/rounding_sweep.R"
Rscript --vanilla "$characterization_root/characterize.R"
