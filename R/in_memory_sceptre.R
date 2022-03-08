#' Run `sceptre` in memory
#'
#' The core function of the `sceptre` package. This function tests for association between a set of gRNAs and a set of genes, returning a p-value for each pairwise test of association.
#'
#' @param gene_matrix a gene-by cell expression matrix; the rows (i.e., gene IDs) should be named
#' @param gRNA_matrix a gRNA-by cell expression matrix; the rows (i.e., gRNA IDs) should be named
#' @param covariate_matrix the cell-specific matrix of technical factors, ideally containing the following covariates: log-transformed gene library size, log-transformed gRNA library size, percent mitochondrial reads, and batch
#' @param gene_gRNA_pairs a data frame specifying the gene-gRNA pairs to test for association; the data frame should contain columns named `gene_id` and `gRNA_id`
#' @param side sidedness of the test; one of "both," "left," and "right"
#' @param B number of resamples to draw for the conditional randomization test
#' @param full_output return the full output (TRUE) or a simplified, reduced output (FALSE)?
#' @param regularization_amount non-negative number specifying the amount of regularization to apply to the negative binomial dispersion parameter estimates
#' @param storage_dir directory in which to store the intermediate computations
#' @param seed seed to the random number generator
#'
#' @return a data frame containing columns `gene_id`, `gRNA_id`, `p_value`, and `z_value`. See
#' @export
#'
#' @examples
#' # 1. load the data
#' data(gene_matrix) # i. gene expression matrix
#' data(gRNA_matrix) # ii. gRNA expression matrix
#' data(covariate_matrix) # iii. covariate matrix
#' data(gRNA_grps) # iv. gRNAs grouped by target site
#' data(gene_gRNA_pairs) # v. gene-gRNA pairs to analyze
#' # 2. (Optional) group together gRNAs that target the same site
#' gRNA_matrix_grouped <- combine_gRNAs(gRNA_matrix, gRNA_grps)
#' # 3. run method (~40s on 8-core Macbook Pro)
#' result <- run_sceptre_in_memory(gene_matrix, gRNA_matrix_grouped, covariate_matrix, gene_gRNA_pairs)
run_sceptre_in_memory <- function(gene_matrix, gRNA_matrix, covariate_matrix, gene_gRNA_pairs, side = "both", storage_dir = tempdir(), regularization_amount = 0.1, B = 1000, full_output = FALSE, seed = 4) {
  ##################
  # DEFINE CONSTANTS
  ##################
  cat(paste0("Check ", storage_dir, "/logs for more detailed status updates.\n"))
  THRESHOLD <- 3
  MIN_GENE_EXP <- 250
  MIN_GRNA_EXP <- 30

  #############################
  # BASIC PROCESSING AND CHECKS
  #############################
  cat("Running checks and setting up directory structure.")
  # 0. Set up parallel, offsite directory structure, pod sizes
  n_cores <- parallel::detectCores()
  doParallel::registerDoParallel(cores = n_cores)
  # print(paste0("Note: check the `log` subdirectory of the storage directory ", as.character(storage_dir), " for updates!"))
  dirs <- initialize_directories(storage_location = storage_dir)

  # 1. threshold gRNA_matrix (using threshold = 3, for now) if necessary
  if (max(gRNA_matrix) >= 2) {
    gRNA_matrix <- gRNA_matrix >= THRESHOLD
  }
  # 2. Remove all genes and gRNAs that have 0 counts; arrange pairs by gRNA then gene; remove duplicates
  gene_lib_sizes <- Matrix::rowSums(gene_matrix)
  gRNA_lib_sizes <- Matrix::rowSums(gRNA_matrix)
  bad_genes <- names(gene_lib_sizes[gene_lib_sizes < MIN_GENE_EXP])
  bad_gRNAs <- names(gRNA_lib_sizes[gRNA_lib_sizes < MIN_GRNA_EXP])
  if (length(bad_genes) >= 1) {
    warning(paste0("Removing genes with low expressions (UMI count <", MIN_GENE_EXP, ")."))
    gene_gRNA_pairs <- gene_gRNA_pairs %>% dplyr::filter(!(gene_id %in% bad_genes))
  }
  if (length(bad_gRNAs) >= 1) {
    warning(paste0("Removing perturbations low counts (counts <", MIN_GRNA_EXP, "). Consider grouping together perturbations that target the same gene."))
    gene_gRNA_pairs <- gene_gRNA_pairs %>% dplyr::filter(!(gRNA_id %in% bad_gRNAs))
  }
  # 3. Make sure genes/gRNAs in the data frame are actually a part of the expression matrices; ensure all rows distinct
  if (!all(c("gene_id", "gRNA_id") %in% colnames(gene_gRNA_pairs))) stop("The columns `gene_id` and `gRNA_id` must be present in the `gene_gRNA_pairs` data frame.")
  abs_genes <- gene_gRNA_pairs$gene_id[!(gene_gRNA_pairs$gene_id %in% row.names(gene_matrix))]
  abs_gRNAs <- gene_gRNA_pairs$gRNA_id[!(gene_gRNA_pairs$gRNA_id %in% row.names(gRNA_matrix))]
  if (length(abs_genes) >= 1) {
    msg <- paste0("The genes `", paste0(abs_genes, collapse = ", "), "' are present in the `gene_gRNA_pairs` data frame but not in the `gene_matrix` matrix. Either remove these genes from `gene_gRNA_pairs,` or add these genes to `gene_matrix`.")
    stop(msg)
  }
  if (length(abs_gRNAs) >= 1) {
    msg <- paste0("The perturbations `", paste0(abs_gRNAs, collapse = ", "), "' are present in the `gene_gRNA_pairs` data frame but not in the `gRNA_matrix` matrix. Either remove these perturbations from `gene_gRNA_pairs,` or add these perturbations to `gRNA_matrix`.")
    stop(msg)
  }
  gene_gRNA_pairs <- gene_gRNA_pairs %>% dplyr::distinct()

  # 4. Ensure that cell barcodes coincide across gene and perturbation matrices
  if  (!identical(colnames(gRNA_matrix), colnames(gene_matrix))) stop("The cell barcodes in the perturbation and expression matrices do not coincide. Ensure that these matrices have identical cell barcodes in the same order.")
  # 5. Set the pods
  n_genes <- length(unique(gene_gRNA_pairs$gene_id))
  n_gRNAs <- length(unique(gene_gRNA_pairs$gRNA_id))
  n_pairs <- nrow(gene_gRNA_pairs)
  pod_sizes <- c(gene = ceiling(n_genes/(2 * n_cores)),
                 gRNA = ceiling(n_gRNAs/(2 * n_cores)),
                 pair = ceiling(n_pairs/(3 * n_cores)))
  cat(crayon::green(' \u2713\n'))

  ##############
  # METHOD START
  ##############
  # create file dictionaries
  dicts <- create_and_store_dictionaries(gene_gRNA_pairs,
                                         dirs[["gene_precomp_dir"]],
                                         dirs[["gRNA_precomp_dir"]],
                                         dirs[["results_dir"]],
                                         pod_sizes)
  cat("Running gene precomputations. ")
  # run first round of gene precomputations
  foreach::`%dopar%`(foreach::foreach(pod_id = seq(1, dicts[["gene"]])),
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
  foreach::`%dopar%`(foreach::foreach(pod_id = seq(1, dicts[["gene"]])),
                     run_gene_precomputation_at_scale_round_2(pod_id = pod_id,
                                                              gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                                              gene_matrix = gene_matrix,
                                                              covariate_matrix = covariate_matrix,
                                                              regularization_amount = regularization_amount,
                                                              log_dir = dirs[["log_dir"]]))

  # run gRNA precomputations
  cat("Running perturbation precomputations. ")
  foreach::`%dopar%`(foreach::foreach(pod_id = seq(1, dicts[["gRNA"]])),
                     run_gRNA_precomputation_at_scale(pod_id = pod_id,
                                                      gRNA_precomp_dir = dirs[["gRNA_precomp_dir"]],
                                                      gRNA_matrix = gRNA_matrix,
                                                      covariate_matrix = covariate_matrix,
                                                      log_dir = dirs[["log_dir"]],
                                                      B = B,
                                                      seed = seed))
  cat(crayon::green(' \u2713\n'))

  # run at-scale analysis
  cat("Running perturbation-to-gene association tests.")
  foreach::`%dopar%`(foreach::foreach(pod_id = seq(1, dicts[["pairs"]])),
                     run_gRNA_gene_pair_analysis_at_scale(pod_id = pod_id,
                                                          gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                                          gRNA_precomp_dir = dirs[["gRNA_precomp_dir"]],
                                                          results_dir = dirs[["results_dir"]],
                                                          log_dir = dirs[["log_dir"]],
                                                          gene_matrix = gene_matrix,
                                                          gRNA_matrix = gRNA_matrix,
                                                          covariate_matrix = covariate_matrix,
                                                          regularization_amount = regularization_amount,
                                                          side = side,
                                                          B = B,
                                                          full_output))
  cat(crayon::green(' \u2713\n'))
  # collect results
  cat("Collecting and returning results. ")
  results <- collect_results(results_dir = dirs[["results_dir"]])
  cat(crayon::green(' \u2713\n'))
  return(results)
}
