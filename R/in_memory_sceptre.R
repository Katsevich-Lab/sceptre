#' Run SCEPTRE on high multiplicity-of-infection single-cell CRISPR screen data
#'
#' This function is the core function of the `sceptre` package. The function applies SCEPTRE to test for association between a set of gRNA groups and genes while controlling for technical confounders. The function returns a p-value for each pairwise test of association.
#'
#' @param gene_matrix a gene-by-cell expression matrix; the rows (i.e., gene IDs) and columns (i.e., cell barcodes) should be named
#' @param combined_perturbation_matrix a binary matrix of perturbations (i.e., gRNA group-to-cell assignments); the rows (i.e., gRNA groups) and columns (i.e., cell barcodes) should be named.
#' @param covariate_matrix the cell-specific matrix of technical factors, ideally containing the following covariates: log-transformed gene library size (numeric), log-transformed gRNA library size (numeric), percent mitochondrial reads (numeric), and batch (factor). The rows (i.e., cell barcodes) should be named
#' @param gene_gRNA_group_pairs a data frame specifying the gene-gRNA group pairs to test for association; the data frame should contain columns named `gene_id` and `gRNA_group`.
#' @param side sidedness of the test; one of "both," "left," and "right"
#' @param B number of resamples to draw for the conditional randomization test
#' @param full_output return the full output (TRUE) or a simplified, reduced output (FALSE; default)?
#' @param regularization_amount non-negative number specifying the amount of regularization to apply to the negative binomial dispersion parameter estimates
#' @param storage_dir directory in which to store the intermediate computations
#' @param parallel parallelize execution?
#' @param seed seed to the random number generator
#'
#' @return the `gene_gRNA_group_pairs` data frame with new columns `p_value` and `z_value` appended. See "details" for a description of the output when `full_output` is set to TRUE.
#'
#' @details
#' Details are arranged from most to least important.
#'
#' - `gene_matrix` should be a **raw** (i.e., un-normalized) matrix of UMI (unique molecular identifier) counts.
#' - `combined_perturbation_matrix` should be a "combined perturbation matrix", which can be obtained by applying the functions `threshold_gRNA_matrix` and `combine_perturbations` (in that order) to a raw gRNA count matrix. `combined_perturbation_matrix` optionally can be a raw gRNA expression matrix or an uncombined perturbation matrix, in which case each gRNA is treated as its own group of one. See the tutorial for more details.
#' - The gene IDs (respectively, gRNA groups) within `gene_gRNA_group_pairs` must be a subset of the row names of `gene_matrix` (respectively, `combined_perturbation_matrix`).
#' - The `side` parameter controls the sidedness of the test. The arguments "left" and "right" are appropriate when testing for a decrease and increase in gene expression, respectively. The default argument -- "both" -- is appropriate when testing for an increase *or* decrease in gene expression.
#' - The default value of `regularization_amount` is 0.1, meaning that a small amount of regularization is applied to the estimated negative binomial size parameters, which helps protect against overfitting. When the number of genes is < 50, however, the default value of `regularization_amount` is set to 0 (i.e., no regularization), as regularization is known to be ineffective when there are few genes.
#' - When `full_output` is set to TRUE (as opposed to FALSE, the default), the output is a data frame with the following columns: `gene_id`, `gRNA_id`, `p_value`, `skew_t_fit_success` (if TRUE, *p*-value based on tail probability of fitted skew-t distribution returned; if FALSE, empirical *p*-value returned), `xi`, `omega`, `alpha`, `nu` (fitted parameters of the skew-t distribution; NA if fit failed), `z_value` (z-value obtained on "ground truth" data), and `z_null_1`, ..., `z_null_B` (z-values obtained from resampled datasets).
#' @export
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(magrittr)
#' # 1. load the data
#' data(gene_matrix) # i. gene expression matrix
#' data(gRNA_matrix) # ii. gRNA expression matrix
#' data(covariate_matrix) # iii. covariate matrix
#' data(gRNA_groups_table) # iv. gRNAs grouped by target site
#' data(gene_gRNA_group_pairs) # v. gene-gRNA group pairs to analyze
#'
#' # 2. threshold and combine gRNA matrix
#' combined_perturbation_matrix <- threshold_gRNA_matrix(gRNA_matrix) %>%
#' combine_perturbations(gRNA_groups_table)
#'
#' # 3. select the gene-gRNA group pairs to analyze
#' set.seed(4)
#' gene_gRNA_group_pairs <- gene_gRNA_group_pairs %>% sample_n(25)
#
#' # 3. run method (takes ~40s on an 8-core Macbook Pro)
#' result <- run_sceptre_high_moi(gene_matrix = gene_matrix,
#' combined_perturbation_matrix = combined_perturbation_matrix,
#' covariate_matrix = covariate_matrix,
#' gene_gRNA_group_pairs = gene_gRNA_group_pairs,
#' side = "left")
#' }
run_sceptre_high_moi <- function(gene_matrix, combined_perturbation_matrix, covariate_matrix, gene_gRNA_group_pairs, side = "both", storage_dir = tempdir(), regularization_amount = 0.1, B = 1000, full_output = FALSE, parallel = TRUE, seed = 4) {
  ##################
  # DEFINE CONSTANTS
  ##################
  cat(paste0("Check ", storage_dir, "/logs for more detailed status updates.\n"))
  THRESHOLD <- 3
  MIN_GENE_EXP <- 250
  MIN_GRNA_EXP <- 30
  N_GENES_REGULARIZATION_THRESH <- 50

  #############################
  # BASIC PROCESSING AND CHECKS
  #############################
  cat("Running checks and setting up directory structure. ")

  # 0. Set up parallel, fst, offsite directory structure, pod sizes
  dirs <- initialize_directories(storage_location = storage_dir)
  if (parallel) {
    n_cores <- parallel::detectCores()
    doParallel::registerDoParallel(cores = n_cores)
    foreach_funct <- foreach::`%dopar%`
  } else {
    n_cores <- 1L
    foreach_funct <- foreach::`%do%`
  }
  # 1. threshold gRNA_matrix (using threshold = 3, for now) if necessary
  gRNA_matrix <- combined_perturbation_matrix
  if (max(gRNA_matrix) >= 2) {
    gRNA_matrix <- gRNA_matrix >= THRESHOLD
  }
  # 2. Ensure that cell barcodes coincide across gene and perturbation matrices
  if  (!identical(colnames(gRNA_matrix), colnames(gene_matrix)) || !identical(row.names(covariate_matrix), colnames(gene_matrix))) {
    stop("The column names of `gene_matrix` and `gRNA_matrix` and the row names of `covariate matrix` should be the cell barcodes. Ensure that the cell barcodes are present and identical across these matrices/data frames.")
  }
  if (is.null(row.names(gRNA_matrix)) || is.null(row.names(gene_matrix))) {
    stop("The row names of `gene_matrix` and `gRNA_matrix` should contain the gene IDs and gRNA IDs, respectively.")
  }

  # 3. Remove all genes and gRNAs that have low counts; arrange pairs by gRNA then gene; remove duplicates
  gene_lib_sizes <- Matrix::rowSums(gene_matrix)
  gRNA_lib_sizes <- Matrix::rowSums(gRNA_matrix)
  bad_genes <- names(gene_lib_sizes[gene_lib_sizes < MIN_GENE_EXP])
  bad_gRNAs <- names(gRNA_lib_sizes[gRNA_lib_sizes < MIN_GRNA_EXP])
  if (length(bad_genes) >= 1) {
    warning(paste0("Removing genes with low expressions (UMI count <", MIN_GENE_EXP, ")."))
    gene_gRNA_group_pairs <- gene_gRNA_group_pairs %>% dplyr::filter(!(gene_id %in% bad_genes))
  }
  if (length(bad_gRNAs) >= 1) {
    warning(paste0("Removing perturbations low counts (counts <", MIN_GRNA_EXP, "). Consider grouping together perturbations that target the same gene."))
    gene_gRNA_group_pairs <- gene_gRNA_group_pairs %>% dplyr::filter(!(gRNA_id %in% bad_gRNAs))
  }
  # 4. Make sure genes/gRNAs in the data frame are actually a part of the expression matrices; ensure all rows distinct
  if ("gRNA_group" %in% colnames(gene_gRNA_group_pairs)) {
    gene_gRNA_group_pairs <- dplyr::rename(gene_gRNA_group_pairs, gRNA_id = gRNA_group)
  }
  if (!all(c("gene_id", "gRNA_id") %in% colnames(gene_gRNA_group_pairs))) stop("The columns `gene_id` and `gRNA_id` must be present in the `gene_gRNA_group_pairs` data frame.")
  abs_genes <- gene_gRNA_group_pairs$gene_id[!(gene_gRNA_group_pairs$gene_id %in% row.names(gene_matrix))]
  abs_gRNAs <- gene_gRNA_group_pairs$gRNA_id[!(gene_gRNA_group_pairs$gRNA_id %in% row.names(gRNA_matrix))]
  if (length(abs_genes) >= 1) {
    msg <- paste0("The genes `", paste0(abs_genes, collapse = ", "), "' are present in the `gene_gRNA_group_pairs` data frame but not in the `gene_matrix` matrix. Either remove these genes from `gene_gRNA_group_pairs,` or add these genes to `gene_matrix`.")
    stop(msg)
  }
  if (length(abs_gRNAs) >= 1) {
    msg <- paste0("The perturbations `", paste0(abs_gRNAs, collapse = ", "), "' are present in the `gene_gRNA_group_pairs` data frame but not in the `gRNA_matrix` matrix. Either remove these perturbations from `gene_gRNA_group_pairs,` or add these perturbations to `gRNA_matrix`.")
    stop(msg)
  }
  gene_gRNA_group_pairs <- gene_gRNA_group_pairs %>% dplyr::distinct()
  # 5. Set the pods
  n_genes <- length(unique(gene_gRNA_group_pairs$gene_id))
  n_gRNAs <- length(unique(gene_gRNA_group_pairs$gRNA_id))
  n_pairs <- nrow(gene_gRNA_group_pairs)
  pod_sizes <- c(gene = ceiling(n_genes/(2 * n_cores)),
                 gRNA = ceiling(n_gRNAs/(2 * n_cores)),
                 pair = ceiling(n_pairs/(2 * n_cores)))
  cat(crayon::green(' \u2713\n'))
  # 6. If number of genes is small, set regularization_amount to 0
  if (n_genes < N_GENES_REGULARIZATION_THRESH) regularization_amount <- 0

  ##############
  # METHOD START
  ##############
  # create file dictionaries
  dicts <- suppressMessages(create_and_store_dictionaries(gene_gRNA_group_pairs,
                                         dirs[["gene_precomp_dir"]],
                                         dirs[["gRNA_precomp_dir"]],
                                         dirs[["results_dir"]],
                                         pod_sizes))
  cat("Running gene precomputations. ")
  # run first round of gene precomputations
  foreach_funct(foreach::foreach(pod_id = seq(1, dicts[["gene"]])),
                     run_gene_precomputation_at_scale_round_1(pod_id = pod_id,
                                                              gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                                              gene_matrix = gene_matrix,
                                                              covariate_matrix = covariate_matrix,
                                                              regularization_amount = regularization_amount,
                                                              log_dir = dirs[["log_dir"]]))
  cat(crayon::green(' \u2713\n'))

  # regularize thetas
  regularize_gene_sizes_at_scale(gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                 regularization_amount = regularization_amount,
                                 log_dir = dirs[["log_dir"]])

  # run second round of gene precomputations
  foreach_funct(foreach::foreach(pod_id = seq(1, dicts[["gene"]])),
                     run_gene_precomputation_at_scale_round_2(pod_id = pod_id,
                                                              gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                                              gene_matrix = gene_matrix,
                                                              covariate_matrix = covariate_matrix,
                                                              regularization_amount = regularization_amount,
                                                              log_dir = dirs[["log_dir"]]))

  # run gRNA precomputations
  cat("Running perturbation precomputations. ")
  foreach_funct(foreach::foreach(pod_id = seq(1, dicts[["gRNA"]])),
                     run_gRNA_precomputation_at_scale(pod_id = pod_id,
                                                      gRNA_precomp_dir = dirs[["gRNA_precomp_dir"]],
                                                      combined_perturbation_matrix = combined_perturbation_matrix,
                                                      covariate_matrix = covariate_matrix,
                                                      log_dir = dirs[["log_dir"]],
                                                      B = B,
                                                      seed = seed))
  cat(crayon::green(' \u2713\n'))

  # run at-scale analysis
  cat("Running perturbation-to-gene association tests.")
  foreach_funct(foreach::foreach(pod_id = seq(1, dicts[["pairs"]])),
                     run_gRNA_gene_pair_analysis_at_scale(pod_id = pod_id,
                                                          gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                                          gRNA_precomp_dir = dirs[["gRNA_precomp_dir"]],
                                                          results_dir = dirs[["results_dir"]],
                                                          log_dir = dirs[["log_dir"]],
                                                          gene_matrix = gene_matrix,
                                                          combined_perturbation_matrix = combined_perturbation_matrix,
                                                          covariate_matrix = covariate_matrix,
                                                          regularization_amount = regularization_amount,
                                                          side = side,
                                                          B = B,
                                                          full_output))
  cat(crayon::green(' \u2713\n'))
  # collect results
  cat("Collecting and returning results. ")
  results <- collect_results(results_dir = dirs[["results_dir"]], gene_gRNA_group_pairs)
  cat(crayon::green(' \u2713\n'))
  return(results)
}
