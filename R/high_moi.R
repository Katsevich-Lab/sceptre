################################################################################
# All high MOI functions in the sceptre package are located in this script.
# We are working on integrating the high MOI portion of the package with the
# low MOI portion of the package and improving the performance of the high MOI
# functions by rewriting them in C++. Thus, most or all of the functions in
# this script will be deprecated in the near future.
################################################################################

#' Activate and deactivate sink
#'
#' A convenience function for programs that run at scale. .highmoi_activate_sink activates the sink to the log file, and .highmoi_deactivate_sink deactivates the sink.
#'
#' @param log_file_name the name of the log file in which to sink the output (including the file path)
#' @noRd
.highmoi_activate_sink <- function(log_file_name) {
  if (file.exists(log_file_name)) file.remove(log_file_name)
  sink(log_file_name)
  sink(stdout(), type = "message")
}


#' @rdname .highmoi_activate_sink
#' @noRd
.highmoi_deactivate_sink <- function() {
  sink(NULL, type = "message")
  sink()
}


#' Initialize directories
#'
#' Initializes the offsite directory structure.
#'
#' @param storage_location file path to a directory in which to store the intermediate and final results
#'
#' @return named vector containing the file paths of gene_precomp_dir, gRNA_precomp_dir, results_dir, and log_dir.
#' @noRd
.highmoi_initialize_directories <- function(storage_location) {
  if (!dir.exists(storage_location)) dir.create(storage_location, recursive = TRUE)
  sub_dirs <- c("gene_precomp", "gRNA_precomp", "results", "logs")
  dirs_to_create <- paste0(storage_location, "/", sub_dirs)
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) dir.create(dir)
  }
  return(c(gene_precomp_dir = dirs_to_create[1], gRNA_precomp_dir = dirs_to_create[2], results_dir = dirs_to_create[3], log_dir = dirs_to_create[4]))
}


#' Create dictionary
#'
#' Create a dictionary that maps named elements (gRNAs, genes, or gRNA pairs) to files for a given "pod" size.
#'
#' @param ids the names of elements (gRNAs or genes)
#' @param pod_size size of the pods
#'
#' @noRd
.highmoi_create_dictionary <- function(ids, pod_size) {
  if (length(ids) <= pod_size) {
    out <- data.frame(id = ids, pod_id = 1)
  } else {
    n_pods_minus_1 <- floor(length(ids)/pod_size)
    pod_id <- rep(seq(1, n_pods_minus_1), each = pod_size)
    pod_id_append <- rep(n_pods_minus_1 + 1, times = length(ids) - length(pod_id))
    pod_id <- c(pod_id, pod_id_append)
    out <- data.frame(id = ids, pod_id = pod_id)
  }
  return(out)
}


#' Create and store dictionaries
#'
#' We require a few bookkeeping files (that we call "dictionaries") for the gRNA and gene precomputation steps.
#' This file creates those dictionaries and stores them in the appropriate locations on disk.
#'
#' @inheritParams run_sceptre_high_moi
#' @param gene_precomp_dir the directory in which to store the gene precomputations
#' @param gRNA_precomp_dir the directory in which to store the gRNA precomputations
#' @param results_dir the directory in which to store the results
#' @param pod_sizes an integer vector with three named elements: gRNA, gene, and pair. These elements give the sizes of the respective "pods."
#'
#' @return an integer vector containing the number of pods in the gene, gRNA, and pairs dictionaries.
#' @noRd
.highmoi_create_and_store_dictionaries <- function(gene_gRNA_group_pairs, gene_precomp_dir, gRNA_precomp_dir, results_dir, pod_sizes) {
  # 0. Clear out contents of gRNA precomp, gene precomp, and results directories
  for (direct in c(gRNA_precomp_dir, gene_precomp_dir)) {
    x <- file.remove(list.files(direct, full.names = TRUE))
  }
  file_names <- list.files(results_dir)
  to_delete <- c(grep(pattern = 'result_[0-9]+.fst', x = file_names, value = TRUE), "results_dictionary.fst")
  for (file in to_delete) {
    full_file <- paste0(results_dir, "/", file)
    if (file.exists(full_file)) file.remove(full_file)
  }

  # 1. genes
  genes <- unique(gene_gRNA_group_pairs$gene_id)
  gene_dictionary <- .highmoi_create_dictionary(ids = genes, pod_size = pod_sizes[["gene"]]) |>
    dplyr::mutate(size_unreg_file = paste0(gene_precomp_dir, "/gene_size_unreg_", pod_id, ".rds") |> factor(),
                  geom_mean_file = paste0(gene_precomp_dir, "/gene_geom_mean_", pod_id, ".rds") |> factor(),
                  offset_file = paste0(gene_precomp_dir, "/gene_offsets_", pod_id, ".fst") |> factor(),
                  size_reg_file =  paste0(gene_precomp_dir, "/size_reg_file.rds") |> factor())
  gene_dictionary_fp <- paste0(gene_precomp_dir, "/gene_dictionary.fst")
  fst::write_fst(gene_dictionary, gene_dictionary_fp)

  # 2. gRNA
  gRNAs <- unique(gene_gRNA_group_pairs$gRNA_id)
  gRNA_dictionary <- .highmoi_create_dictionary(ids = gRNAs, pod_size = pod_sizes[["gRNA"]]) |>
    dplyr::mutate(precomp_file = factor(paste0(gRNA_precomp_dir, "/gRNA_precomp_", id, ".rds")))
  gRNA_dictionary_fp <- paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst")
  fst::write_fst(gRNA_dictionary, gRNA_dictionary_fp)

  # 3. gRNA-gene pairs
  pairs_dictionary <- .highmoi_create_dictionary(ids = gene_gRNA_group_pairs$gRNA_id, pod_size = pod_sizes[["pair"]]) |>
    dplyr::rename(gRNA_id = id) |> dplyr::mutate(gene_id = gene_gRNA_group_pairs$gene_id) |> dplyr::select(gRNA_id, gene_id, pod_id) |>
    dplyr::mutate(result_file = paste0(results_dir, "/result_", pod_id, ".fst") |> factor())
  pairs_dictionary_fp <- paste0(results_dir, "/results_dictionary.fst")
  fst::write_fst(pairs_dictionary, pairs_dictionary_fp)

  out <- c(gene = gene_dictionary$pod_id[nrow(gene_dictionary)], gRNA = gRNA_dictionary$pod_id[nrow(gRNA_dictionary)], pairs = pairs_dictionary$pod_id[nrow(pairs_dictionary)])
  return(out)
}


#' Run gRNA precomputation at scale
#'
#' This function runs the gRNA precomputation on a selected "pod" of gRNAs (as identified by the pod_id).
#' It stores the result in the gRNA precomputation directory.
#'
#' @inheritParams run_sceptre_high_moi
#' @param pod_id ID of the pod for which to do the precomputation
#' @param gRNA_precomp_dir file path to the gRNA precomputation directory
#' @param log_dir file path to the log directory
#'
#' @return NULL
#' @noRd
.highmoi_run_gRNA_precomputation_at_scale <- function(pod_id, gRNA_precomp_dir, combined_perturbation_matrix, covariate_matrix, log_dir, B, seed) {
  # Activate the sink for the log file
  if (!is.null(log_dir)) .highmoi_activate_sink(paste0(log_dir, "/gRNA_precomp_", pod_id, ".Rout"))
  # determine the gRNAs on which to run the precomputation
  gRNA_dictionary <- fst::read_fst(paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst")) |> dplyr::filter(pod_id == !!pod_id)
  gRNA_ids <- gRNA_dictionary |> dplyr::pull(id) |> as.character()
  # run the precomputation for each of these gRNAs
  for (i in seq(1, length(gRNA_ids))) {
    gRNA_id <- gRNA_ids[i]
    cat(paste0("Running precomputation for gRNA ", gRNA_id, ".\n"))
    gRNA_indicators <- if (any(methods::is(gene_matrix) %in% c("matrix", "Matrix"))) {
      combined_perturbation_matrix[gRNA_id,]
    } else {
      combined_perturbation_matrix[[gRNA_id,]] |> as.numeric()
    }
    synth_data <- .highmoi_run_gRNA_precomputation(gRNA_indicators, covariate_matrix, B, seed)
    precomp_matrix_fp <- (gRNA_dictionary |> dplyr::pull(precomp_file))[i] |> as.character()
    saveRDS(object = synth_data, precomp_matrix_fp)
  }
  if (!is.null(log_dir)) .highmoi_deactivate_sink()
}


#' .highmoi_run_gene_precomputation_at_scale_round_1
#'
#' This function runs the first round of gene precomputations. In particular, it computes the raw dispersion estimate
#' and log geometric mean of each gene. It saves the results in the gene_precomp directory.
#'
#' @inheritParams run_sceptre_high_moi
#' @param pod_id pod id
#' @param gene_precomp_dir location of the gene precomputation directory
#' @param log_dir directory in which to sink the log file
#' @noRd
.highmoi_run_gene_precomputation_at_scale_round_1 <- function(pod_id, gene_precomp_dir, gene_matrix, covariate_matrix, regularization_amount, log_dir) {
  if (!is.null(log_dir)) .highmoi_activate_sink(paste0(log_dir, "/gene_precomp_round_1_pod_", pod_id, ".Rout"))

  gene_dictionary <- fst::read_fst(paste0(gene_precomp_dir, "/gene_dictionary.fst")) |> dplyr::filter(pod_id == !!pod_id)
  gene_ids <- gene_dictionary |> dplyr::pull(id) |> as.character()
  precomps <- purrr::map(gene_ids, function(gene_id) {
    cat(paste0("Running precomputation round 1 for gene ", gene_id, ".\n"))
    expressions <- if (any(methods::is(gene_matrix) %in% c("matrix", "Matrix"))) {
      gene_matrix[gene_id,]
    } else {
      gene_matrix[[gene_id,]] |> as.numeric()
    }
    unreg_size <- .highmoi_run_gene_precomputation(expressions = expressions, covariate_matrix = covariate_matrix, gene_precomp_size = NULL)[["gene_precomp_size"]]
    out <- list()
    out$unreg_size <- unreg_size
    if (regularization_amount > 0) {
      gene_log_geom_mean <- .highmoi_log_geom_mean(expressions)
      out$gene_log_geom_mean <- gene_log_geom_mean
    }
    return(out)
  })

  names(precomps) <- gene_ids
  out_unreg_sizes <- purrr::map_dbl(precomps, function(l) l$unreg_size)
  if (regularization_amount > 0) out_log_geom_means <- purrr::map_dbl(precomps, function(l) l$gene_log_geom_mean)

  unreg_sizes_save_fp <- (gene_dictionary |> dplyr::pull(size_unreg_file))[1] |> as.character()
  saveRDS(object = out_unreg_sizes, file = unreg_sizes_save_fp)
  if (regularization_amount > 0) {
    geom_means_save_fp <- (gene_dictionary |> dplyr::pull(geom_mean_file))[1] |> as.character()
    saveRDS(object = out_log_geom_means, file = geom_means_save_fp)
  }

  if (!is.null(log_dir)) .highmoi_deactivate_sink()
}


#' Regularize genes at scale
#'
#' Regularizes the estimated gene sizes.
#'
#' @param gene_precomp_dir location of the gene precomputation directory
#' @param log_dir location of the directory in which to sink the log file
#' @param regularization_amount amount of regularization to apply to thetas (>= 0)
#' @noRd
.highmoi_regularize_gene_sizes_at_scale <- function(gene_precomp_dir, regularization_amount, log_dir) {
  if (regularization_amount > 0) {
    if (!is.null(log_dir)) .highmoi_activate_sink(paste0(log_dir, "/gene_precomp_regularize_sizes", ".Rout"))
    file_names <- list.files(gene_precomp_dir)
    gene_size_unreg_files <- paste0(gene_precomp_dir, "/", grep(pattern = 'gene_size_unreg_[0-9]+.rds', x = file_names, value = TRUE))
    gene_geom_mean_files <- paste0(gene_precomp_dir, "/", grep(pattern = 'gene_geom_mean_[0-9]+.rds', x = file_names, value = TRUE))
    load_distributed_vector <- function(f_names) f_names |> purrr::map(readRDS) |> purrr::reduce(c)
    gene_sizes_unreg <- load_distributed_vector(gene_size_unreg_files)
    gene_geom_means <- load_distributed_vector(gene_geom_mean_files)
    if (!all(names(gene_sizes_unreg) == names(gene_geom_means))) {
      gene_sizes_unreg <- gene_sizes_unreg[order(names(gene_sizes_unreg))]
      gene_geom_means <- gene_geom_means[order(names(gene_geom_means))]
    }
    cat("Regularizing gene sizes.\n")
    sizes_reg <- .highmoi_regularize_thetas(genes_log_gmean = gene_geom_means, theta = gene_sizes_unreg, plot_me = FALSE)
    names(sizes_reg) <- names(gene_sizes_unreg)
    sizes_reg_save_fp <- (fst::read_fst(paste0(gene_precomp_dir, "/gene_dictionary.fst"), "size_reg_file") |> dplyr::pull())[1] |> as.character()
    saveRDS(object = sizes_reg, file = sizes_reg_save_fp)
    if (!is.null(log_dir)) .highmoi_deactivate_sink()
  }
}


#' Run gene precomputation at scale round 2
#'
#' Runs the second round of gene precomputations.
#'
#' @inheritParams run_sceptre_high_moi
#' @param pod_id pod id
#' @param gene_precomp_dir gene precomp dir
#' @param gene_matrix the gene expression matrix, stored as an ondisc_matrix
#' @param log_dir directory in which to sink the logs
#' @noRd
.highmoi_run_gene_precomputation_at_scale_round_2 <- function(pod_id, gene_precomp_dir, gene_matrix, covariate_matrix, regularization_amount, log_dir) {
  if (!is.null(log_dir)) .highmoi_activate_sink(paste0(log_dir, "/gene_precomp_round_2_pod_", pod_id, ".Rout"))

  gene_dictionary <- fst::read_fst(paste0(gene_precomp_dir, "/gene_dictionary.fst")) |> dplyr::filter(pod_id == !!pod_id)
  gene_ids <- gene_dictionary |> dplyr::pull(id) |> as.character()
  if (regularization_amount > 0) {
    gene_sizes <- readRDS(as.character(gene_dictionary$size_reg_file[1]))[gene_ids]
  } else {
    gene_sizes <- readRDS(as.character(gene_dictionary$size_unreg_file[1]))
  }

  offsets <- sapply(gene_ids, function(gene_id) {
    cat(paste0("Running precomputation round 2 for gene ", gene_id, ".\n"))
    expressions <- if (any(methods::is(gene_matrix) %in% c("matrix", "Matrix"))) {
      gene_matrix[gene_id,]
    } else {
      gene_matrix[[gene_id,]] |> as.numeric()
    }
    dist_offsets <- .highmoi_run_gene_precomputation(expressions = expressions, covariate_matrix = covariate_matrix, gene_precomp_size = gene_sizes[[gene_id]])[["gene_precomp_offsets"]]
    return(dist_offsets)
  }) |> as.data.frame()

  offsets_save_fp <- (gene_dictionary |> dplyr::pull(offset_file))[1] |> as.character()
  fst::write_fst(x = offsets, path = offsets_save_fp)
  if (!is.null(log_dir)) .highmoi_deactivate_sink()
}


#' Run gRNA-gene pair analysis at scale
#'
#' Runs the gene-gRRNA pair ananysis across an entire pod of gene-gRNA pairs.
#'
#' @inheritParams run_sceptre_high_moi
#' @param pod_id id of the pod to compute
#' @param gene_precomp_dir the directory containing the results of the gene precomputation
#' @param gRNA_precomp_dir the directory containing the results of the gRNA precomputation
#' @param results_dir directory in which to store the results
#' @param log_dir (optional) directory in which to sink the log file
#' @return NULL
#' @noRd
.highmoi_run_gRNA_gene_pair_analysis_at_scale <- function(pod_id, gene_precomp_dir, gRNA_precomp_dir, results_dir, log_dir, gene_matrix, combined_perturbation_matrix, covariate_matrix, regularization_amount, side, B, full_output) {
  if (!is.null(log_dir)) .highmoi_activate_sink(paste0(log_dir, "/result_", pod_id, ".Rout"))

  results_dict <- fst::read_fst(paste0(results_dir, "/results_dictionary.fst")) |> dplyr::filter(pod_id == !!pod_id)
  gene_dict <- fst::read_fst(paste0(gene_precomp_dir, "/gene_dictionary.fst"))
  gRNA_dict <- fst::read_fst(paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst"))
  if (regularization_amount > 0) {
    regularized_gene_sizes <- readRDS(gene_dict$size_reg_file[1] |> as.character())[as.character(results_dict$gene_id)]
  }
  out_l <- vector(mode = "list", length = nrow(results_dict))

  for (i in seq(1, nrow(results_dict))) {
    curr_gene <- results_dict[[i, "gene_id"]] |> as.character()
    curr_gRNA <- results_dict[[i, "gRNA_id"]] |> as.character()
    cat(paste0("Running distilled CRT on gene ", curr_gene, " and gRNA ", curr_gRNA, ".\n"))

    # Load the appropriate data from disk into memory
    # 1. gene
    if (i == 1 || results_dict[[i, "gene_id"]] != results_dict[[i - 1, "gene_id"]]) {
      gene_precomp_locs <- dplyr::filter(gene_dict, id == curr_gene)
      gene_offset_loc <- gene_precomp_locs |> dplyr::pull(offset_file) |> as.character()
      gene_size_loc <- gene_precomp_locs |> dplyr::pull(size_unreg_file) |> as.character()
      gene_precomp_offsets <- fst::read_fst(path = gene_offset_loc, columns = curr_gene) |> dplyr::pull()
      if (regularization_amount > 0) {
        gene_precomp_size <- regularized_gene_sizes[[curr_gene]]
      } else {
        gene_precomp_size <- readRDS(file = gene_size_loc)[[curr_gene]]
      }
      expressions <- if (any(methods::is(gene_matrix) %in% c("matrix", "Matrix"))) {
        gene_matrix[curr_gene,]
      } else {
        gene_matrix[[curr_gene,]] |> as.numeric()
      }
    }

    # 2. gRNA
    if (i == 1 || results_dict[[i, "gRNA_id"]] != results_dict[[i - 1, "gRNA_id"]]) {
      gRNA_prcomp_loc <- dplyr::filter(gRNA_dict, id == curr_gRNA) |> dplyr::pull(precomp_file) |> as.character()
      gRNA_precomp <- readRDS(file = gRNA_prcomp_loc)
      gRNA_indicators <- if (any(methods::is(combined_perturbation_matrix) %in% c("matrix", "Matrix"))) {
        combined_perturbation_matrix[curr_gRNA,]
      } else {
        combined_perturbation_matrix[[curr_gRNA,]] |> as.numeric()
      }
    }

    # Run the dCRT
    out_l[[i]] <- .highmoi_run_sceptre_using_precomp_fast(expressions, gRNA_indicators, gRNA_precomp, side,
                                                 gene_precomp_size, gene_precomp_offsets, full_output)
  }
  # Create and save the result dataframe
  out_df <- data.table::rbindlist(out_l) |>
    dplyr::mutate(gRNA_id = results_dict$gRNA_id, gene_id = results_dict$gene_id) |>
    dplyr::relocate(gene_id, gRNA_id)
  out_fp <- (results_dict |> dplyr::pull(result_file))[1] |> as.character()
  fst::write_fst(out_df, out_fp)
  if (!is.null(log_dir)) .highmoi_deactivate_sink()
}


#' Collect results
#'
#' Collates the individual results files into a single result file.
#'
#' @param results_dir the directory containing the results
#' @param gene_gRNA_group_pairs the gene-gRNA groups data frame
#'
#' @return the collated results data frame.
#' @noRd
.highmoi_collect_results <- function(results_dir, gene_gRNA_group_pairs) {
  file_names <- list.files(results_dir)
  to_load <- grep(pattern = 'result_[0-9]+.fst', x = file_names, value = TRUE)
  all_results <- results_dir |> paste0("/", to_load) |> purrr::map(fst::read_fst) |>  data.table::rbindlist()
  if (ncol(gene_gRNA_group_pairs) >= 3) {
    all_results <- dplyr::left_join(gene_gRNA_group_pairs, all_results)
  }
  saveRDS(object = all_results, file = paste0(results_dir, "/all_results.rds"))
  return(all_results)
}


#' Run `sceptre` on high MOI data
#'
#' The function applies SCEPTRE to test for association between a set of gRNA groups and genes in a high MOI experiment. The function returns a p-value and log fold change estimate for each pairwise test of association.
#'
#' @param gene_matrix a gene-by-cell expression matrix; the rows (i.e., gene IDs) and columns (i.e., cell barcodes) should be named
#' @param combined_perturbation_matrix a binary matrix of perturbations (i.e., gRNA group-to-cell assignments); the rows (i.e., gRNA groups) and columns (i.e., cell barcodes) should be named.
#' @param covariate_matrix the cell-specific matrix of technical factors, ideally containing the following covariates: log-transformed gene library size (numeric), log-transformed gRNA library size (numeric), percent mitochondrial reads (numeric), and batch (factor). The rows (i.e., cell barcodes) should be named
#' @param gene_gRNA_group_pairs a data frame specifying the gene-gRNA group pairs to test for association; the data frame should contain columns named `gene_id` and `gRNA_group`.
#' @param side sidedness of the test; one of "both," "left," and "right". (default "both")
#' @param B number of resamples to draw for the conditional randomization test. (default 1000)
#' @param full_output return the full output (TRUE) or a simplified, reduced output (FALSE)? (default FALSE)
#' @param regularization_amount non-negative number specifying the amount of regularization to apply to the negative binomial dispersion parameter estimates (default 0)
#' @param storage_dir directory in which to store the intermediate computations (default tempdir)
#' @param parallel parallelize execution? (default TRUE)
#' @param seed seed to the random number generator (default 4)
#'
#' @return the `gene_gRNA_group_pairs` data frame with new columns `p_value`, `z_value`, and `log_fold_change` appended. See "details" for a description of the output when `full_output` is set to TRUE.
#'
#' @details
#' Details are arranged from most to least important.
#'
#' - `gene_matrix` should be a **raw** (i.e., un-normalized) matrix of UMI (unique molecular identifier) counts.
#' - `combined_perturbation_matrix` should be a "combined perturbation matrix", which can be obtained by applying the functions `threshold_gRNA_matrix` and `combine_perturbations` (in that order) to a raw gRNA count matrix. `combined_perturbation_matrix` optionally can be a raw gRNA expression matrix or an uncombined perturbation matrix, in which case each gRNA is treated as its own group of one. See the tutorial for more details.
#' - The gene IDs (respectively, gRNA groups) within `gene_gRNA_group_pairs` must be a subset of the row names of `gene_matrix` (respectively, `combined_perturbation_matrix`).
#' - The `side` parameter controls the sidedness of the test. The arguments "left" and "right" are appropriate when testing for a decrease and increase in gene expression, respectively. The default argument -- "both" -- is appropriate when testing for an increase *or* decrease in gene expression.
#' - The default value of `regularization_amount` is 0.0, meaning that zero regularization is applied to the estimated negative binomial size parameters. One can increase the value of this parameter to protect against overfitting, which can be useful when there are many genes.
#' - When `full_output` is set to TRUE (as opposed to FALSE, the default), the output is a data frame with the following columns: `gene_id`, `gRNA_id`, `p_value`, `skew_t_fit_success` (if TRUE, *p*-value based on tail probability of fitted skew-t distribution returned; if FALSE, empirical *p*-value returned), `xi`, `omega`, `alpha`, `nu` (fitted parameters of the skew-t distribution; NA if fit failed), `z_value` (z-value obtained on "ground truth" data), and `z_null_1`, ..., `z_null_B` (z-values obtained from resampled datasets).
#' @export
#' @examples
#' \dontrun{
#' library(dplyr)
#' # 1. load the data
#' data(gene_matrix) # i. gene expression matrix
#' data(gRNA_matrix) # ii. gRNA expression matrix
#' data(covariate_matrix) # iii. covariate matrix
#' data(gRNA_groups_table) # iv. gRNAs grouped by target site
#' data(gene_gRNA_group_pairs) # v. gene-gRNA group pairs to analyze
#'
#' # 2. threshold and combine gRNA matrix
#' combined_perturbation_matrix <- threshold_gRNA_matrix(gRNA_matrix) |>
#' combine_perturbations(gRNA_groups_table)
#'
#' # 3. select the gene-gRNA group pairs to analyze
#' set.seed(4)
#' gene_gRNA_group_pairs <- gene_gRNA_group_pairs |> sample_n(25)
#
#' # 3. run method (takes ~40s on an 8-core Macbook Pro)
#' result <- run_sceptre_highmoi(gene_matrix = gene_matrix,
#' combined_perturbation_matrix = combined_perturbation_matrix,
#' covariate_matrix = covariate_matrix,
#' gene_gRNA_group_pairs = gene_gRNA_group_pairs,
#' side = "left")
#' }
run_sceptre_highmoi <- function(gene_matrix, combined_perturbation_matrix, covariate_matrix, gene_gRNA_group_pairs, side = "both", storage_dir = tempdir(), regularization_amount = 0.0, B = 1000, full_output = FALSE, parallel = FALSE, seed = 4) {
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
  dirs <- .highmoi_initialize_directories(storage_location = storage_dir)
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
    warning(paste0("Removing the following genes with low expressions (UMI count <", MIN_GENE_EXP, "): ", paste0(bad_gRNAs, collapse = " ")))
    gene_gRNA_group_pairs <- gene_gRNA_group_pairs |> dplyr::filter(!(gene_id %in% bad_genes))
  }
  if (length(bad_gRNAs) >= 1) {
    warning(paste0("Removing the following perturbations with low counts (counts <", MIN_GRNA_EXP, "): ", paste0(bad_gRNAs, collapse = " "), ". Consider grouping together perturbations that target the same site to circumvent this problem."))
    gene_gRNA_group_pairs <- gene_gRNA_group_pairs |> dplyr::filter(!(gRNA_group %in% bad_gRNAs))
  }
  # 4. Make sure genes/gRNAs in the data frame are actually a part of the expression matrices; ensure all rows distinct
  if (!all(c("gene_id", "gRNA_group") %in% colnames(gene_gRNA_group_pairs))) stop("The columns `gene_id` and `gRNA_group` must be present in the `gene_gRNA_group_pairs` data frame.")
  # swap "gRNA_group" to "gRNA_id"
  gene_gRNA_group_pairs <- dplyr::rename(gene_gRNA_group_pairs, gRNA_id = gRNA_group)
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
  gene_gRNA_group_pairs <- gene_gRNA_group_pairs |> dplyr::distinct()
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
  dicts <- suppressMessages(.highmoi_create_and_store_dictionaries(gene_gRNA_group_pairs,
                                                          dirs[["gene_precomp_dir"]],
                                                          dirs[["gRNA_precomp_dir"]],
                                                          dirs[["results_dir"]],
                                                          pod_sizes))
  cat("Running gene precomputations. ")
  # run first round of gene precomputations
  foreach_funct(foreach::foreach(pod_id = seq(1, dicts[["gene"]])),
                .highmoi_run_gene_precomputation_at_scale_round_1(pod_id = pod_id,
                                                         gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                                         gene_matrix = gene_matrix,
                                                         covariate_matrix = covariate_matrix,
                                                         regularization_amount = regularization_amount,
                                                         log_dir = dirs[["log_dir"]]))
  cat(crayon::green(' \u2713\n'))

  # regularize thetas
  .highmoi_regularize_gene_sizes_at_scale(gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                          regularization_amount = regularization_amount,
                                          log_dir = dirs[["log_dir"]])

  # run second round of gene precomputations
  foreach_funct(foreach::foreach(pod_id = seq(1, dicts[["gene"]])),
                .highmoi_run_gene_precomputation_at_scale_round_2(pod_id = pod_id,
                                                         gene_precomp_dir = dirs[["gene_precomp_dir"]],
                                                         gene_matrix = gene_matrix,
                                                         covariate_matrix = covariate_matrix,
                                                         regularization_amount = regularization_amount,
                                                         log_dir = dirs[["log_dir"]]))

  # run gRNA precomputations
  cat("Running perturbation precomputations. ")
  foreach_funct(foreach::foreach(pod_id = seq(1, dicts[["gRNA"]])),
                .highmoi_run_gRNA_precomputation_at_scale(pod_id = pod_id,
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
                .highmoi_run_gRNA_gene_pair_analysis_at_scale(pod_id = pod_id,
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
  results <- .highmoi_collect_results(results_dir = dirs[["results_dir"]], gene_gRNA_group_pairs)
  cat(crayon::green(' \u2713\n'))
  return(results)
}


#' Fit skew-t
#'
#' Fits a skew-t distribution on a set of resampled test statistics.
#'
#' @param t_nulls sampled statistics that form a null distribution
#' @param t_star the observed test statistic
#' @param side the side of the test (left, right, of both)
#'
#' @return
#' A list containing (i) skew_t_fit_success (boolean), (ii) out_p (the p-value), and (iii) skew-t mle (a vector containing the fitted MLEs, NA if fit failed).
#' @noRd
.highmoi_fit_skew_t <- function(t_nulls, t_star, side) {
  p_value_skew_t <- NA
  skew_t_fit <- tryCatch(
    R.utils::withTimeout(expr = sn::selm(t_nulls ~ 1, family = "ST"), timeout = 3),
    error = function(e) return(NA),
    warning = function(w) return(NA))
  if (methods::is(skew_t_fit, "selm")) { # If the fit worked,
    dp <- skew_t_fit@param$dp # then extract the parameters.
    if (!any(is.na(dp))) { # If all the fitted parameters are numbers,
      p_value_skew_t <- switch(side,
                               'left' = pmax(.Machine$double.eps, sn::pst(x = t_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)), # then compute the skew t-based p-value. pst(x = t_star, dp = dp)
                               'right' = pmax(.Machine$double.eps, 1 - sn::pst(x = t_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)),
                               'both' = pmax(.Machine$double.eps, sn::pst(x = -abs(t_star), dp = dp, method = 2, rel.tol = .Machine$double.eps) +
                                               (1 - sn::pst(x = abs(t_star), dp = dp, method = 2, rel.tol = .Machine$double.eps)))
      )
    }
  }
  # check if the skew-t fit worked
  skew_t_fit_success <- !is.na(p_value_skew_t)
  if (skew_t_fit_success) {
    out_p <- p_value_skew_t
    skew_t_mle <- dp
  } else {
    out_p <- switch(side,
                    'left' = mean(c(-Inf, t_nulls) <= t_star),
                    'right' = mean(c(Inf, t_nulls) >= t_star),
                    'both' = mean(c(Inf, abs(t_nulls)) >= abs(t_star)))
    skew_t_mle <- c(xi = NA, omega = NA, alpha = NA, nu = NA)
  }
  return(list(skew_t_fit_success = skew_t_fit_success, out_p = out_p, skew_t_mle = skew_t_mle))
}


#' Run sceptre using precomputations for gRNAs and genes.
#'
#' This function is the workhorse function of the sceptre package. It runs a distilled CRT using a negative binomial test statistic based on an expression vector, a gRNA indicator vector, an offset vector (from the distillation step), gRNA conditional probabilities, an estimate of the negative binomial size parameter, and the number of resampling replicates.
#'
#' @param expressions a vector of gene expressions (in UMI counts)
#' @param gRNA_indicators a vector of gRNA indicators
#' @param gRNA_precomp a n_cells x B matrix of synthetic indicators
#' @param gene_precomp_size the pre-computed size
#' @param gene_precomp_offsets the pre-computed distillation offsets
#' @param side an argument to set test for left-sided, right-sided or two-tailed. Default as 'left' and can take 'left', 'right' or 'both'.#'
#' @noRd
#' @return if reduced_output, a dataframe with the p-value, test statistic, and other helpful values; if !reduced_output, a list containing all of the above, plus the resampled statistics.
.highmoi_run_sceptre_using_precomp_fast <- function(expressions, gRNA_indicators, gRNA_precomp, side, gene_precomp_size, gene_precomp_offsets, full_output = FALSE) {
  exp_gene_offsets <- exp(gene_precomp_offsets)
  # compute test statistic on real data
  y <- expressions[gRNA_indicators == 1]
  exp_o <- exp_gene_offsets[gRNA_indicators == 1]
  log_fold_change <- log(mean(y)) - log(mean(exp_o))
  z_star <- .highmoi_compute_nb_test_stat_fast_score(y, exp_o, gene_precomp_size)
  # compute test statistics on the resampled statistics
  z_null <- apply(X = gRNA_precomp, MARGIN = 2, FUN = function(col) {
    .highmoi_compute_nb_test_stat_fast_score(expressions[col], exp_gene_offsets[col], gene_precomp_size)
  })
  # fit skew-t
  skew_t_fit <- .highmoi_fit_skew_t(z_null, z_star, side)
  # create output
  if (full_output) {
    B <- length(z_null)
    z_df <- as.data.frame(matrix(z_null, nrow = 1, ncol = B))
    colnames(z_df) <- paste0("z_null_", seq(1, B))
    out <- data.frame(p_value = skew_t_fit$out_p, skew_t_fit_success = skew_t_fit$skew_t_fit_success,
                      xi = skew_t_fit$skew_t_mle[["xi"]], omega = skew_t_fit$skew_t_mle[["omega"]],
                      alpha = skew_t_fit$skew_t_mle[["alpha"]], nu = skew_t_fit$skew_t_mle[["nu"]],
                      z_value = z_star, log_fold_change = log_fold_change) |> dplyr::mutate(z_df)
  } else {
    out <- data.frame(p_value = skew_t_fit$out_p, z_value = z_star, log_fold_change = log_fold_change)
  }
  return(out)
}


.highmoi_compute_nb_test_stat_fast_score <- function(y, exp_o, gene_precomp_size) {
  r_exp_o <- gene_precomp_size * exp_o
  y_exp_o <- y * exp_o
  r_plus_exp_o <- gene_precomp_size + exp_o
  sum_y <- sum(y)
  top <- (y_exp_o + r_exp_o)/r_plus_exp_o
  bottom <- r_exp_o/r_plus_exp_o
  z <- (sum_y - sum(top))/sqrt(sum(bottom))
  return(z)
}


#' Run gRNA precomputation
#'
#' This function runs the precomputation for a given gRNA.
#'
#' @param gRNA_indicators a vector of gRNA indicators
#' @param covariate_matrix the cell-specific covariate matrix
#'
#' @noRd
#' @return the fitted probabilities
.highmoi_run_gRNA_precomputation <- function(gRNA_indicators, covariate_matrix, B, seed) {
  set.seed(seed)
  # first, fit a logistic regression model to estimate perturbation probabilities
  fit_model_grna <- stats::glm(gRNA_indicators ~ ., family = stats::binomial(), data = covariate_matrix)
  out <- as.numeric(stats::fitted(fit_model_grna))
  # generate synthetic data
  draws_idxs <- lapply(X = out, FUN = function(prob) {
    draws <- stats::rbinom(n = B, size = 1, prob = prob)
    which(draws == 1)
  })
  n_cells <- length(gRNA_indicators)
  row_idxs <- rep(x = seq(1L, n_cells), times = sapply(draws_idxs, length))
  col_idxs <- unlist(draws_idxs)
  synth_data <- Matrix::sparseMatrix(i = row_idxs, j = col_idxs, dims = c(n_cells, B))
  return(synth_data)
}


#' Run gene precomputation
#'
#' This function runs the precomputation for a given gene. It fits an NB regression of expression on covariates. The estimate of theta (i.e., the NB size parameter) is obtained from the glm.nb function. Offsets are obtained by log-transforming the fitted values.
#'
#' @param expressions the vector of gene expressions
#' @param covariate_matrix the cell-specific covariate matrix
#' @param gene_precomp_size the pre-computed size parameter
#'
#' @return a named list containing two items: offsets and size.
#' @noRd
.highmoi_run_gene_precomputation <- function(expressions, covariate_matrix, gene_precomp_size) {
  # cases on gene_precomp_size
  if (is.null(gene_precomp_size)) {
    # no size supplied; use glm.nb to estimate size and fit model
    # second backup: method of moments on poisson reg
    backup_2 <- function(pois_fit) {
      MASS::theta.mm(y = expressions, mu = pois_fit$fitted.values, dfr = pois_fit$df.residual)
    }

    # first backup: MLE on poisson reg
    backup <- function() {
      pois_fit <- stats::glm(expressions ~ ., data = covariate_matrix, family = stats::poisson())
      gene_precomp_size_out <- tryCatch({
        MASS::theta.ml(expressions, pois_fit$fitted.values, limit = 50)[1]
      }, error = function(e) backup_2(pois_fit), warning = function(w) backup_2(pois_fit))
      fit_nb <- VGAM::vglm(formula = expressions ~ ., family = VGAM::negbinomial.size(gene_precomp_size_out), data = covariate_matrix)
      fitted_vals <- as.numeric(VGAM::fittedvlm(fit_nb))
      list(fitted_vals = fitted_vals, gene_precomp_size_out = gene_precomp_size_out)
    }

    # try to fit a negative binomial GLM with unknown dispersion
    result <- tryCatch({
      fit_nb <- MASS::glm.nb(formula = expressions ~ ., data = covariate_matrix)
      fitted_vals <- as.numeric(fit_nb$fitted.values)
      gene_precomp_size_out <- fit_nb$theta
      list(fitted_vals = fitted_vals, gene_precomp_size_out = gene_precomp_size_out)
    }, error = function(e) backup(), warning = function(w) backup())

    fitted_vals <- result$fitted_vals; gene_precomp_size_out <- result$gene_precomp_size_out

  } else { # size supplied; use vglm to fit model
    gene_precomp_size_out <- gene_precomp_size
    fit_nb <- VGAM::vglm(formula = expressions ~ ., family = VGAM::negbinomial.size(gene_precomp_size_out), data = covariate_matrix)
    fitted_vals <- as.numeric(VGAM::fittedvlm(fit_nb))
  }

  gene_precomp_offsets_out <- log(fitted_vals)
  out <- list(gene_precomp_offsets = gene_precomp_offsets_out, gene_precomp_size = gene_precomp_size_out)
  return(out)
}


#' Run `sceptre` on a gRNA-gene pair
#'
#' Runs `sceptre` on a given gRNA-gene pair.
#'
#' This function is for demonstration purposes only. The function `run_sceptre_in_memory` should be used when testing for association between multiple genes and gRNAs.
#'
#' @param gene_expressions (numeric vector) a vector of gene expressions
#' @param gRNA_expressions (numeric vector) a vector of gRNA expressions (or, optionally, a user-thresholded vector of binary gRNA indicators)
#' @param covariate_matrix (data frame) the cell-specific matrix of covariates
#' @param side sidedness of the test, one of "left," "right," and "both"
#' @param B number of resamples for the conditional randomization test
#' @param full_output return the full output (TRUE) or a streamlined, reduced output (FALSE)?
#' @param seed seed to the random number generator
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' # load the example data
#' data(gene_matrix); data(gRNA_matrix); data(covariate_matrix)
#' gene_expressions <- gene_matrix[1,]
#' gRNA_expressions <- gRNA_matrix[1,]
#' # run method
#' result <- .highmoi_run_sceptre_gRNA_gene_pair(gene_expressions, gRNA_expressions, covariate_matrix, "left")
#' }
.highmoi_run_sceptre_gRNA_gene_pair <- function(gene_expressions, gRNA_expressions, covariate_matrix, side = "both", B = 1500, full_output = FALSE, seed = 4) {
  THRESHOLD <- 3

  cat("Running perturbation precomputation.")
  # threshold gRNA_expressions, if necessary
  if (max(gRNA_expressions) >= 2) gRNA_expressions <- as.integer(gRNA_expressions >= THRESHOLD)
  gene_precomp <- .highmoi_run_gene_precomputation(gene_expressions, covariate_matrix, NULL)
  cat(crayon::green(' \u2713\n'))

  cat("Running gRNA precomputation.")
  gRNA_precomp <- .highmoi_run_gRNA_precomputation(gRNA_expressions, covariate_matrix, B, seed)
  cat(crayon::green(' \u2713\n'))

  cat("Running perturbation-to-gene association analysis.")
  out <- .highmoi_run_sceptre_using_precomp_fast(gene_expressions,
                                        gRNA_expressions,
                                        gRNA_precomp,
                                        side,
                                        gene_precomp$gene_precomp_size,
                                        gene_precomp$gene_precomp_offsets,
                                        full_output)
  cat(crayon::green(' \u2713\n'))
  return(out)
}


#' Plot result
#'
#' For a given gRNA-gene pair analyzed by `sceptre`, plots the resampled test statistics alongside the "ground truth" test statistic derived from the raw data.
#'
#' @param row single row of the data frame outputted by `run_sceptre_high_moi`, when `full_output` is set to TRUE
#'
#' @return a ggplot2 object containing the plot
#' @examples
#' # RUN THE METHOD
#' set.seed(4)
#' library(dplyr)
#' data(gene_matrix)
#' data(gRNA_matrix)
#' data(covariate_matrix)
#' data(gRNA_groups_table)
#' data(gene_gRNA_group_pairs)
#' combined_perturbation_matrix <- threshold_gRNA_matrix(gRNA_matrix) |>
#' combine_perturbations(gRNA_groups_table)
#' gene_gRNA_group_pairs <- gene_gRNA_group_pairs |> sample_n(1)
#' result <- run_sceptre_high_moi(gene_matrix = gene_matrix,
#' combined_perturbation_matrix = combined_perturbation_matrix,
#' covariate_matrix = covariate_matrix,
#' gene_gRNA_group_pairs = gene_gRNA_group_pairs,
#' side = "left", parallel = FALSE, full_output = TRUE)
#'
#' # CREATE THE PLOT
#' plot_result(result[1,])
#' @noRd
plot_result <- function(row) {
  resampled_zvalues <- row |> dplyr::select(dplyr::starts_with("z_null_")) |> as.numeric()
  original_zvalue <- row$z_value
  interval <- range(c(resampled_zvalues, original_zvalue)) + c(-0.5, 0.5)
  dp <- row |> dplyr::select(xi, omega, alpha, nu) |> as.numeric()
  z <- seq(interval[1], interval[2], length.out = 1000)
  df_curves <- data.frame(z = z, fitted = sn::dst(x = z, dp = dp), gaussian = stats::dnorm(z)) |>
    tidyr::gather("curve", "y", fitted, gaussian) |>
    dplyr::mutate(curve = factor(curve, levels = c("fitted","gaussian"), labels = c("Conditional\nrandomization","Negative\nbinomial")))
  df_ribbon <- df_curves |>
    dplyr::filter(z <= original_zvalue, curve == "Conditional\nrandomization") |>
    dplyr::mutate(lower = 0, upper = y) |>
    dplyr::select(z, lower, upper)
  p <- ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = z, y = ..density..),
                            data = data.frame(z = resampled_zvalues),
                            boundary = 0, colour = "black", fill = "lightgray", bins = 15) +
    ggplot2::geom_line(ggplot2::aes(x = z, y = y, group = curve, colour = curve, linetype = curve), data = df_curves) +
    ggplot2::geom_vline(xintercept = original_zvalue, colour = "firebrick3", linetype = "solid") +
    ggplot2::geom_ribbon(ggplot2::aes(x = z, ymin = lower, ymax = upper), fill = "darkorchid2", alpha = 0.5, data = df_ribbon) +
    ggplot2::scale_colour_manual(values = c("darkorchid2", "black"), name = "Null distribution") +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed"), name = "Null distribution") +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::xlab("Negative binomial z-value") + ggplot2::theme_bw() +
    ggplot2::theme(legend.position = c(0.85, 0.8),
                   legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
  return(p)
}


#' Make QQ-plot
#'
#' Makes a QQ-plot for a set of values hypothesized to follow a uniform distribution (e.g., *p*-values).
#'
#' @param p_values a vector values -- hypothesized to follow a uniform distribution -- from which to construct the QQ-plot. This vector typically will be a set of *p*-values.
#' @param ci_level level of the pointwise confidence band (default 0.95)
#' @param point_col color of the plotted points
#' @param alpha transparency of the plotted points
#'
#' @return a ggplot object containing the QQ-plot
#' @export
#' @examples
#' set.seed(4)
#' p_vals <- runif(5000)
#' make_qq_plot(p_vals)
make_qq_plot <- function(p_values, ci_level = 0.95, point_col = "royalblue4", alpha = 0.9) {
  p_thresh <- 1e-8
  to_plot <- data.frame(pvalue = p_values) |>
    dplyr::mutate(r = rank(pvalue), expected = stats::ppoints(dplyr::n())[r],
                  clower = stats::qbeta(p = (1 - ci_level)/2, shape1 = r, shape2 = dplyr::n() + 1 - r),
                  cupper = stats::qbeta(p = (1 + ci_level)/2, shape1 = r, shape2 = dplyr::n()+ 1 - r)) |>
    # dplyr::filter(-log10(expected) > 2) |>
    dplyr::mutate(pvalue = ifelse(pvalue < p_thresh, p_thresh, pvalue))
  p <- ggplot2::ggplot(data = to_plot, mapping = ggplot2::aes(x = expected, y = pvalue, ymin = clower, ymax = cupper)) +
    ggplot2::geom_point(size = 1, alpha = alpha, col = point_col) +
    ggplot2::geom_ribbon(alpha = 0.25) +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::scale_x_continuous(trans = revlog_trans(base = 10)) + ggplot2::scale_y_continuous(trans = revlog_trans(base = 10)) +
    ggplot2::xlab(expression(paste("Expected null p-value"))) +
    ggplot2::ylab(expression(paste("Observed p-value"))) +
    ggplot2::ggtitle("QQ-plot of p-values") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = c(0.25, 0.85),
                                         legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
                                         legend.title = ggplot2::element_blank(),
                                         panel.grid = ggplot2::element_blank(),
                                         strip.background = ggplot2::element_blank(),
                                         panel.border = ggplot2::element_blank(),
                                         axis.line = ggplot2::element_line(),
                                         plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_vline(xintercept = 0.1, col = "darkred", linetype = "dashed")
  return(p)
}


#### NOTE: These functions are taken from the sctransform package. Credit goes to Christoph Hafemeister and Rahul Satija, whom we cited in our high MOI sceptre paper. We had to rework these functions slightly to improve their scalability.
#' Regularize thetas
#'
#' @param genes_log_gmean the log-transformed vector of the geometric mean of the gene expressions
#' @param theta the vector of unregularized gene sizes
#' @param theta_regularization transformation to apply to the thetas before fitting nonparametric regression; options "log_theta" and "od_factor;" defaults to "log_theta"
#' @param bw_adjust amount of regularization (greater value -> more regularization). 0 corresponds to no regularization at all. Default value 3.
#' @param plot_me create informative plots? (default FALSE)
#'
#' @return the regularized vector of thetas
#' @noRd
#'
#' @examples
#' genes_log_gmean <- c(runif(1000, 2, 5))
#' names(genes_log_gmean) <- paste0("gene_", 1:1000)
#' theta <- genes_log_gmean * 0.5 + rnorm(100)
#' theta_regularization <- "log_theta"
#' bw_adjust <- 3
#' theta_reg <- .highmoi_regularize_thetas(genes_log_gmean, theta, theta_regularization, 0, TRUE)
#### NOTE: These functions are taken from the sctransform package. Credit goes to Christoph Hafemeister and Rahul Satija, whom we cited in our high MOI sceptre paper. We had to rework these functions slightly to improve their scalability.

#' Regularize thetas
#'
#' @param genes_log_gmean the log-transformed vector of the geometric mean of the gene expressions
#' @param theta the vector of unregularized gene sizes
#' @param theta_regularization transformation to apply to the thetas before fitting nonparametric regression; options "log_theta" and "od_factor;" defaults to "log_theta"
#' @param bw_adjust amount of regularization (greater value -> more regularization). 0 corresponds to no regularization at all. Default value 3.
#' @param plot_me create informative plots? (default FALSE)
#'
#' @return the regularized vector of thetas
#' @noRd
#'
#' @examples
#' genes_log_gmean <- c(runif(1000, 2, 5))
#' names(genes_log_gmean) <- paste0("gene_", 1:1000)
#' theta <- genes_log_gmean * 0.5 + rnorm(100)
#' theta_regularization <- "log_theta"
#' bw_adjust <- 3
#' theta_reg <- .highmoi_regularize_thetas(genes_log_gmean, theta, theta_regularization, 0, TRUE)
.highmoi_regularize_thetas <- function(genes_log_gmean, theta, theta_regularization = "log_theta", bw_adjust = 3, plot_me = FALSE) {
  if (bw_adjust == 0) {
    theta_out <- theta
  } else {
    theta <- pmax(theta, 1e-7)
    dispersion_par <- switch(theta_regularization, log_theta = log10(theta), od_factor = log10(1 + 10^genes_log_gmean/theta), stop("theta_regularization ", theta_regularization, " unknown - only log_theta and od_factor supported at the moment"))
    outliers <- .highmoi_is_outlier(dispersion_par, genes_log_gmean)

    dispersion_par <- dispersion_par[!outliers]
    genes_log_gmean_step1 <- genes_log_gmean[!outliers]

    bw <- stats::bw.SJ(genes_log_gmean_step1) * bw_adjust
    x_points <- pmax(genes_log_gmean, min(genes_log_gmean_step1))
    x_points <- pmin(x_points, max(genes_log_gmean_step1))
    o <- order(x_points)
    fitted_dispersion_par <- numeric(length = length(x_points))
    fitted_dispersion_par[o] <- stats::ksmooth(x = genes_log_gmean_step1, y = dispersion_par, x.points = x_points, bandwidth = bw, kernel = "normal")$y
    theta_out <- switch(theta_regularization, log_theta = 10^fitted_dispersion_par, od_factor = 10^genes_log_gmean/(10^fitted_dispersion_par - 1))
    problem_idxs <- is.na(theta_out)
    theta_out[problem_idxs] <- theta[problem_idxs]
    attr(theta_out, "outliers") <- outliers

    if (plot_me) {
      df1 <- data.frame(log_theta = c(dispersion_par, fitted_dispersion_par), mean_gene_exp = c(genes_log_gmean_step1, x_points), regularized = c(rep(FALSE, length(genes_log_gmean_step1)), rep(TRUE, length(x_points))))
      p1 <- ggplot2::ggplot(data = df1, mapping = ggplot2::aes(x = mean_gene_exp, y = log_theta, col = regularized)) + ggplot2::geom_point() + ggplot2::theme_bw()
      df2 <- data.frame(theta = theta, theta_out = theta_out)
      p2 <- ggplot2::ggplot(data = df2, mapping = ggplot2::aes(x = theta, y = theta_out)) + ggplot2::geom_point() + ggplot2::theme_bw()
      print(p1); print(p2)
    }
  }
  return(theta_out)
}


#' Internal sctransform function
#' @noRd
.highmoi_robust_scale_binned <- function (y, x, breaks) {
  bins <- cut(x = x, breaks = breaks, ordered_result = TRUE)
  tmp <- stats::aggregate(x = y, by = list(bin = bins), FUN = .highmoi_robust_scale)
  score <- rep(0, length(x))
  o <- order(bins)
  if (inherits(x = tmp$x, what = "list")) {
    score[o] <- unlist(tmp$x)
  }
  else {
    score[o] <- as.numeric(t(tmp$x))
  }
  return(score)
}


#' Internal sctransform function
#' @noRd
.highmoi_robust_scale <- function(x) (x - stats::median(x))/(stats::mad(x) + .Machine$double.eps)


#' Log geometric mean
#'
#' @param x a numeric vector
#' @param eps small number to add to each value (default 1)
#'
#' @return the log-transformed geometric means
#' @noRd
.highmoi_log_geom_mean <- function(x, eps = 1) {log10(exp(mean(log(x + eps))) - eps)}

#' Internal sctransform function
#' @noRd
.highmoi_is_outlier <- function(y, x, th = 10) {
  bin.width <- (max(x) - min(x)) * stats::bw.SJ(x)/2
  eps <- .Machine$double.eps * 10
  breaks1 <- seq(from = min(x) - eps, to = max(x) + bin.width,
                 by = bin.width)
  breaks2 <- seq(from = min(x) - eps - bin.width/2, to = max(x) +
                   bin.width, by = bin.width)
  score1 <- .highmoi_robust_scale_binned(y, x, breaks1)
  score2 <- .highmoi_robust_scale_binned(y, x, breaks2)
  return(pmin(abs(score1), abs(score2)) > th)
}


#' Internal sctransform function
#' @noRd
.highmoi_robust_scale_binned <- function (y, x, breaks) {
  bins <- cut(x = x, breaks = breaks, ordered_result = TRUE)
  tmp <- stats::aggregate(x = y, by = list(bin = bins), FUN = .highmoi_robust_scale)
  score <- rep(0, length(x))
  o <- order(bins)
  if (inherits(x = tmp$x, what = "list")) {
    score[o] <- unlist(tmp$x)
  }
  else {
    score[o] <- as.numeric(t(tmp$x))
  }
  return(score)
}


#' Internal sctransform function
#' @noRd
.highmoi_robust_scale <- function(x) (x - stats::median(x))/(stats::mad(x) + .Machine$double.eps)


#' Log geometric mean
#'
#' @param x a numeric vector
#' @param eps small number to add to each value (default 1)
#'
#' @return the log-transformed geometric means
#' @noRd
.highmoi_log_geom_mean <- function(x, eps = 1) {log10(exp(mean(log(x + eps))) - eps)}

####### END NOTE

#' Run gene precomputation (v2)
#'
#' Runs precomputation on a gene
#'
#' @param expressions the numeric vector of gene expressions
#' @param covariate_matrix the covariate matrix on which to regress (NOTE: should contain an interecept term)
#' @param fam a string indicating the family to use in the regression (either "nb" or "poisson")
#'
#' @return a named vector of fitted coefficients, alongside the fitted size parameter (named "gene_theta") if fam == "nb"
#' @noRd
.highmoi_run_gene_precomputation_v2 <- function(expressions, covariate_matrix, fam) {
  # second backup: method of moments
  backup_2 <- function(pois_fit) {
    MASS::theta.mm(y = expressions, mu = pois_fit$fitted.values, dfr = pois_fit$df.residual)
  }

  # first backup: MLE on poisson reg
  backup <- function() {
    pois_fit <- stats::glm(expressions ~ . + 0, data = covariate_matrix, family = stats::poisson())
    gene_theta <- tryCatch({
      MASS::theta.ml(expressions, pois_fit$fitted.values, limit = 50)[1]
    }, error = function(e) backup_2(pois_fit), warning = function(w) backup_2(pois_fit))
    gene_theta <- max(gene_theta, 0.1)
    fit_nb <- VGAM::vglm(formula = expressions ~ . + 0, family = VGAM::negbinomial.size(gene_theta), data = covariate_matrix)
    fitted_coefs <- stats::coef(fit_nb)
    return(c(fitted_coefs = fitted_coefs, gene_theta = gene_theta))
  }

  # try to fit a negative binomial GLM with unknown dispersion
  if (fam == "nb") {
    result <- tryCatch({
      fit_nb_init <- MASS::glm.nb(formula = expressions ~ . + 0, data = covariate_matrix)
      gene_theta <- max(fit_nb_init$theta, 0.1)
      fit_nb <- VGAM::vglm(formula = expressions ~ . + 0, family = VGAM::negbinomial.size(gene_theta), data = covariate_matrix)
      fitted_coefs <- stats::coef(fit_nb)
      return(c(fitted_coefs, gene_theta = gene_theta))
    }, error = function(e) backup(), warning = function(w) backup())
  }
  if (fam == "poisson") {
    pois_fit <- stats::glm(expressions ~ . + 0, data = covariate_matrix, family = stats::poisson())
    result <- stats::coef(pois_fit)
  }
  return(result)
}


#' Run grna precomputation (v2)
#'
#' Runs precomputation on a vector of grna indicators
#'
#' @param indicators the binary vector of grna indicators
#' @param covariate_matrix the covariate matrix on which to regress (NOTE: should contain an interecept term)
#'
#' @return a named vector of fitted coefficients
#' @noRd
.highmoi_run_grna_precomputation_v2 <- function(indicators, covariate_matrix) {
  fit <- stats::glm(formula = indicators ~ . + 0, family = "binomial", data = covariate_matrix)
  ret <- stats::coef(fit)
  return(ret)
}


#' Generate synthetic data
#'
#' Generate a (sparse) matrix of synthetic grna indicators from a vector of fitted probabilitites.
#'
#' @param fitted_probs a vector of fitted grna presence/absence indicators
#' @param B the number of resampled vectors to create
#'
#' @return a sparse matrix (of dimension n x B) containing the resampled indicators
#' @noRd
.highmoi_generate_synthetic_grna_data <- function(fitted_probs, B) {
  draws_idxs <- lapply(X = fitted_probs, FUN = function(prob) {
    draws <- stats::rbinom(n = B, size = 1, prob = prob)
    which(draws == 1)
  })
  n_cells <- length(fitted_probs)
  row_idxs <- rep(x = seq(1L, n_cells), times = sapply(draws_idxs, length))
  col_idxs <- unlist(draws_idxs)
  synth_data <- Matrix::sparseMatrix(i = row_idxs, j = col_idxs, dims = c(n_cells, B))
  return(synth_data)
}


#' Combine perturbations
#'
#' Combines perturbations by collapsing gRNAs within the same group into a single "combined" gRNA.
#'
#' The function combines binary (i.e., 0/1) perturbation vectors via a "max" operation. In other words, if `x1`, `x2`, ..., `xp` are binary perturbation vectors, then the "combined" perturbation vector `v` is `v = pmax(x1, x2, ..., xp)`, where `pmax` is the element-wise maximum function.
#'
#' @param perturbation_matrix a perturbation matrix (i.e., a binary matrix of gRNA-to-cell assignments) stored as a sparse matrix (as implemented by the Matrix package) or a dense matrix (as implemented by base R); the row names should be the gRNA IDs.
#' @param gRNA_groups_table a data frame with columns `gRNA_id` and `gRNA_group`, mapping each gRNA to its gRNA group
#'
#' @return a "combined" perturbation matrix, where gRNAs within the same gRNA group have been collapsed into a single row
#' @export
#'
#' @examples
#' library(Matrix)
#' data("gRNA_matrix_highmoi")
#' data("gRNA_groups_table_highmoi")
#' gRNA_matrix <- gRNA_matrix_highmoi
#' gRNA_groups_table <- gRNA_groups_table_highmoi
#' perturbation_matrix <- threshold_gRNA_matrix(gRNA_matrix)
#' combined_perturbation_matrix <- combine_perturbations(perturbation_matrix, gRNA_groups_table)
combine_perturbations <- function(perturbation_matrix, gRNA_groups_table) {
  # check that gRNA_id and gRNA_group are present in `gRNA_groups` df
  if (!all(c("gRNA_id", "gRNA_group") %in% colnames(gRNA_groups_table))) {
    stop("`gRNA_id` and `gRNA_group` should be columns of `gRNA_groups`.")
  }

  # convert the data frame to a named list
  gRNA_groups <- as.character(unique(gRNA_groups_table$gRNA_group))
  gRNA_groups <- stats::setNames(gRNA_groups, gRNA_groups)
  gRNA_grp_list <- lapply(X = gRNA_groups, FUN = function(gRNA_group) {
    dplyr::filter(gRNA_groups_table, gRNA_group == !!gRNA_group) |> dplyr::pull(gRNA_id)
  })
  # Determine which gRNAs are not contained within gRNA_grps
  leftover_gRNAs <- setdiff(row.names(perturbation_matrix), unlist(gRNA_grp_list))
  out_leftover <- perturbation_matrix[leftover_gRNAs,]
  out_grped <- sapply(X = names(gRNA_grp_list), FUN = function(grp_name) {
    mat_sub <- perturbation_matrix[gRNA_grp_list[[grp_name]],,drop=FALSE]
    Matrix::colSums(mat_sub)
  }) |> Matrix::t()
  out_grped <- out_grped >= 1
  out <- rbind(out_leftover, out_grped)
}

#' Threshold gRNA count matrix
#'
#' Thresholds a gRNA count matrix, producing a perturbation matrix.
#'
#' @param gRNA_matrix a gRNA-by-cell expression matrix; the matrix can be represented as a sparse matrix (as implemented by the Matrix package) or a dense matrix (as implemented by base R)
#' @param threshold the threshold used to assign perturbation indicators to cells; counts above the threshold are set to 1 (indicating "perturbed"), and counts below the threshold are set to 0 (indicating "unperturbed")
#'
#' @return a binary matrix of perturbation assignments
#' @export
#'
#' @examples
#' data(gRNA_matrix_highmoi)
#' gRNA_matrix <- gRNA_matrix_highmoi
#' perturbation_matrix <- threshold_gRNA_matrix(gRNA_matrix)
threshold_gRNA_matrix <- function(gRNA_matrix, threshold = 3) {
  return(gRNA_matrix >= threshold)
}
