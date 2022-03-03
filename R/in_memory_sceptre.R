#' Run sceptre in memory
#'
#' Runs sceptre in memory. This function is appropriate for data of intermediate size (i.e., the data are small enough to fit into memory, and the analysis does not need to be run across nodes on a computer cluster).
#'
#' @param storage_dir name of a directory in which to store logs, intermediate computations, and final results
#' @param expression_matrix a matrix of gene expressions (in UMI counts); rows correspond to genes and columns to cells; the rows should be named.
#' @param perturbation_matrix a binary matrix of gRNA perturbations; rows correspond to gRNAs and columns to cells; the rows likewise should be names.
#' @param covariate_matrix a data frame of covariates; rows correspond to cells.
#' @param gene_gRNA_pairs a data frame with named columns "gene_id" and "gRNA_id" giving the names of the gRNAs and genes to analyze.
#' @param side sidedness of the test; one of "left", "right," and "both".  "left" is most appropriate for experiments in which cis-regulatory relationships are tested by perturbing putative enhancers with CRISPRi.
#' @param pod_sizes an integer vector giving the size of the "gene", "gRNA", and "pair" pods. `sceptre` groups the genes, gRNAs, and gene-gRNA pairs into distinct "pods" and runs computations on these pods in parallel. Smaller pod sizes give rise to greater parallelization. `pod_sizes` should be vector of length three with names "gene", "gRNA", and "pair" and integer values giving the pod sizes.
#' @param regularization_amount (optional; default 0.1) the amount of regularization to apply to the estimated negative binomial size parameters, where 0 corresponds to no regularization at all.
#' @param seed (optional; default 4) seed to pass to the random number generator.
#' @param B (optional; default 500) number of random samples to draw in the conditional randomization test.
#'
#' @return a data frame containing the results. Includes columns gene_id, gRNA_id, p_value, skew_t_fit_success (indicating whether the skew-t fit succeeded), xi, omega, alpha, nu (parameters of the skew-t distribution), z_value (the ground truth negative binomial test statistic), and n_successful_resamples (indicates how many of the B resamples were successful).
#' @export
#'
#' @examples
#' \dontrun{
#' # generate random perturbation and expression data in memory
#' library(dplyr)
#' set.seed(4)
#' n_gRNAs <- 50
#' n_genes <- 40
#' n_cells <- 5000
#' # create perturbation, expression, and covariate matrices
#' perturbation_matrix <- replicate(n = n_gRNAs, rbinom(n_cells, 1, 0.05)) %>% t()
#' expression_matrix <- replicate(n = n_genes, rpois(n_cells, 1)) %>% t()
#' covariate_matrix <- data.frame(p_mito = runif(n = n_cells, min = 0, max = 10),
#'                                lg_umi_count = log(rpois(n = n_cells, lambda =  500)))
#' # assign column names to the perturbation and expression matrices
#' row.names(perturbation_matrix) <- paste0("gRNA", seq(1, n_gRNAs))
#' row.names(expression_matrix) <- paste0("gene", seq(1, n_genes))
#' # select 5000 random gene-gRNA pairs to analyze; create the gene_gRNA_pairs data frame
#' gene_gRNA_pairs <- expand.grid(gene_id = row.names(expression_matrix),
#'                                gRNA_id = row.names(perturbation_matrix)) %>% slice_sample(n = 90)
#' # let the temporary directory be the storage dir
#' storage_dir <- tempdir()
#' # set the remaining parameters: test sidedness, pod_sizes, regularization_amount, seed, and B.
#' side <- "left"
#' pod_sizes = c(gene = 10, gRNA = 10, pair = 15)
#' regularization_amount <- 0.1
#' seed <- 4
#' B <- 500
#' result <- run_sceptre_in_memory(storage_dir,
#'                      expression_matrix,
#'                      perturbation_matrix,
#'                      covariate_matrix,
#'                      gene_gRNA_pairs,
#'                      side,
#'                      pod_sizes,
#'                      regularization_amount,
#'                      seed,
#'                      B)
#' }
run_sceptre_in_memory <- function(storage_dir, expression_matrix, perturbation_matrix, covariate_matrix, gene_gRNA_pairs, side, pod_sizes = c(gene = 10, gRNA = 10, pair = 10), regularization_amount = 0.1, seed = 4, B = 500) {
  print(paste0("Note: check the `log` subdirectory of the storage directory ", as.character(storage_dir), " for updates!"))
  # set up parallel
  doParallel::registerDoParallel()
  # create the offsite directory structure
  dirs <- initialize_directories(storage_location = storage_dir)

  # BASIC PROCESSING AND CHECKS
  # 1. threshold perturbation_matrix (using threshold = 3, for now) if necessary
  if (max(perturbation_matrix) >= 2) {
    perturbation_matrix <- perturbation_matrix >= 3
  }
  # 2. Remove all genes and gRNAs that have 0 counts
  gene_lib_sizes <- Matrix::rowSums(expression_matrix)
  gRNA_lib_sizes <- Matrix::rowSums(perturbation_matrix)
  bad_genes <- names(gene_lib_sizes[gene_lib_sizes == 0])
  bad_gRNAs <- names(gRNA_lib_sizes[gRNA_lib_sizes == 0])
  gene_gRNA_pairs <- gene_gRNA_pairs %>% dplyr::filter(!(gene_id %in% bad_genes),
                                                       !(gRNA_id %in% bad_gRNAs))
  # 3. Make sure genes/gRNAs in the data frame are actually a part of the expression matrices
  abs_genes <- gene_gRNA_pairs$gene_id[!(gene_gRNA_pairs$gene_id %in% row.names(expression_matrix))]
  abs_gRNAs <- gene_gRNA_pairs$gRNA_id[!(gene_gRNA_pairs$gRNA_id %in% row.names(perturbation_matrix))]
  if (length(abs_genes) >= 1) {
    msg <- paste0("The genes `", paste0(abs_genes, collapse = ", "), "' are present in the `gene_gRNA_pairs` data frame but not in the `expression_matrix` matrix. Either remove these genes from `gene_gRNA_pairs,` or add these genes to `expression_matrix`.")
    stop(msg)
  }
  if (length(abs_gRNAs) >= 1) {
    msg <- paste0("The perturbations `", paste0(abs_gRNAs, collapse = ", "), "' are present in the `gene_gRNA_pairs` data frame but not in the `perturbation_matrix` matrix. Either remove these perturbations from `gene_gRNA_pairs,` or add these perturbations to `perturbation_matrix`.")
    stop(msg)
  }

  # create file dictionaries
  dicts <- create_and_store_dictionaries(gene_gRNA_pairs,
                                         dirs[["gene_precomp_dir"]],
                                         dirs[["gRNA_precomp_dir"]],
                                         dirs[["results_dir"]],
                                         pod_sizes)

  # run first round of gene precomputations
  foreach::`%dopar%`(foreach::foreach(pod_id = seq(1, dicts[["gene"]])),
                     run_gene_precomputation_at_scale_round_1(pod_id = pod_id,
                                                              gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                                              expression_matrix = expression_matrix,
                                                              covariate_matrix = covariate_matrix,
                                                              regularization_amount = regularization_amount,
                                                              log_dir = dirs[["log_dir"]]))

  # regularize thetas
  regularize_gene_sizes_at_scale(gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                 regularization_amount = regularization_amount,
                                 log_dir = dirs[["log_dir"]])

  # run second round of gene precomputations
  foreach::`%dopar%`(foreach::foreach(pod_id = seq(1, dicts[["gene"]])),
                     run_gene_precomputation_at_scale_round_2(pod_id = pod_id,
                                                              gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                                              expression_matrix = expression_matrix,
                                                              covariate_matrix = covariate_matrix,
                                                              regularization_amount = regularization_amount,
                                                              log_dir = dirs[["log_dir"]]))

  # run gRNA precomputations
  foreach::`%dopar%`(foreach::foreach(pod_id = seq(1, dicts[["gRNA"]])),
                     run_gRNA_precomputation_at_scale(pod_id = pod_id,
                                                      gRNA_precomp_dir = dirs[["gRNA_precomp_dir"]],
                                                      perturbation_matrix = perturbation_matrix,
                                                      covariate_matrix = covariate_matrix,
                                                      log_dir = dirs[["log_dir"]]))

  # run at-scale analysis
  foreach::`%dopar%`(foreach::foreach(pod_id = seq(1, dicts[["pairs"]])),
                     run_gRNA_gene_pair_analysis_at_scale(pod_id = pod_id,
                                                          gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                                          gRNA_precomp_dir = dirs[["gRNA_precomp_dir"]],
                                                          results_dir = dirs[["results_dir"]],
                                                          log_dir = dirs[["log_dir"]],
                                                          expression_matrix = expression_matrix,
                                                          perturbation_matrix = perturbation_matrix,
                                                          covariate_matrix = covariate_matrix,
                                                          regularization_amount = regularization_amount,
                                                          side = side,
                                                          seed = seed,
                                                          B = B))

  # collect results
  results <- collect_results(results_dir = dirs[["results_dir"]])
  return(results)
}
