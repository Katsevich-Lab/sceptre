utils::globalVariables(c("n_nonzero_trt", "n_nonzero_cntrl", "pair_str", "assignment", "g",
                         "multiple_grnas", "x", "grna_id", "grna_expressions_bin", "bin_counts",
                         "significant", "lab", "p_values", "pass_qc", "fraction_cells_removed",
                         "grna_id", "pass_qc"))

#' sceptre
#'
#' `sceptre` is a statistically principled, fast, memory-light, and user-friendly software for single-cell CRISPR screen analysis. sceptre achieves state-of-the-art calibration and power on single-cell CRISPR screen data by leveraging several methodological and algorithmic advances in assumption-lean and computationally efficient differential expression analysis. `sceptre` includes modules for the analysis of both low MOI and high MOI data.
#'
#' @useDynLib sceptre, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import BH
#' @docType package
#' @name sceptre
NULL
