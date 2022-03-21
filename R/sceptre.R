utils::globalVariables(c("log_theta", "mean_gene_exp", "regularized", "pod_id", "id", "gRNA_id", "gene_id", "offset_file", "size_unreg_file", "precomp_file", "result_file", "geom_mean_file", "curve", "fitted", "gaussian", "lower", "upper", "..density..", "y", "%dopar%", "foreach", "registerDoParallel", "xi", "omega", "alpha", "nu", "gene_size_loc", "gene_matrix", "pvalue", "r", "expected", "clower", "cupper", "gRNA_group"))

#' sceptre
#'
#' SCEPTRE (analysis of single cell perturbation screens via conditional resampling) is a method for single-cell CRISPR screen analysis. SCEPTRE infers gene-perturbation associations by modeling the random assortment of CRISPR guide RNAs among cells instead of modeling gene expression, thereby remaining valid despite confounder presence and arbitrary misspecification of the gene expression model.
#'
#' @docType package
#' @importFrom magrittr %>%
#' @import methods
#' @name sceptre
NULL
