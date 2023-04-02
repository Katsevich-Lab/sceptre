utils::globalVariables(c("log_theta", "mean_gene_exp", "regularized", "pod_id", "id", "gRNA_id", "gene_id", "offset_file", "size_unreg_file", "precomp_file", "result_file", "geom_mean_file", "curve", "fitted", "gaussian", "lower", "upper", "..density..", "y", "%dopar%", "foreach", "registerDoParallel", "xi", "omega", "alpha", "nu", "gene_size_loc", "gene_matrix", "pvalue", "r", "expected", "clower", "cupper", "gRNA_group", c("grna_group", "p_value", "response_id", "B2", "B3", "full_test_stat", "p_adj", "reject", "log_2_fold_change", "density")))

#' sceptre
#'
#' `sceptre` is a resource-light, fast, and statistically principled software for single-cell CRISPR screen analysis. `sceptre` is powered by several methodological and algorithmic advances in robust differential expression analysis and fast permutation testing. `sceptre` has modules for low multiplicity of infection (MOI) and high MOI single-cell CRISPR screen analysis.
#'
#' @useDynLib sceptre, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import BH
#' @docType package
#' @name sceptre
NULL
