utils::globalVariables(c("log_theta", "mean_gene_exp", "regularized", "pod_id", "id", "gRNA_id", "gene_id", "offset_file", "size_unreg_file", "precomp_file", "result_file", "geom_mean_file", "curve", "fitted", "gaussian", "lower", "upper", "..density..", "y", "%dopar%", "foreach", "registerDoParallel", "xi", "omega", "alpha", "nu", "gene_size_loc", "gene_matrix", "pvalue", "r", "expected", "clower", "cupper", "gRNA_group", "grna_group", "p_value", "response_id", "B2", "B3", "full_test_stat", "p_adj", "reject", "log_2_fold_change", "density"))

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
