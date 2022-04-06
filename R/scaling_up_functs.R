#' Activate and deactivate sink
#'
#' A convenience function for programs that run at scale. activate_sink activates the sink to the log file, and deactivate_sink deactivates the sink.
#'
#' @param log_file_name the name of the log file in which to sink the output (including the file path)
#' @noRd
activate_sink <- function(log_file_name) {
  if (file.exists(log_file_name)) file.remove(log_file_name)
  sink(log_file_name)
  sink(stdout(), type = "message")
}


#' @rdname activate_sink
#' @noRd
deactivate_sink <- function() {
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
initialize_directories <- function(storage_location) {
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
create_dictionary <- function(ids, pod_size) {
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
create_and_store_dictionaries <- function(gene_gRNA_group_pairs, gene_precomp_dir, gRNA_precomp_dir, results_dir, pod_sizes) {
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
  gene_dictionary <- create_dictionary(ids = genes, pod_size = pod_sizes[["gene"]]) %>%
    dplyr::mutate(size_unreg_file = paste0(gene_precomp_dir, "/gene_size_unreg_", pod_id, ".rds") %>% factor(),
           geom_mean_file = paste0(gene_precomp_dir, "/gene_geom_mean_", pod_id, ".rds") %>% factor(),
           offset_file = paste0(gene_precomp_dir, "/gene_offsets_", pod_id, ".fst") %>% factor(),
           size_reg_file =  paste0(gene_precomp_dir, "/size_reg_file.rds") %>% factor())
  gene_dictionary_fp <- paste0(gene_precomp_dir, "/gene_dictionary.fst")
  fst::write_fst(gene_dictionary, gene_dictionary_fp)

  # 2. gRNA
  gRNAs <- unique(gene_gRNA_group_pairs$gRNA_id)
  gRNA_dictionary <- create_dictionary(ids = gRNAs, pod_size = pod_sizes[["gRNA"]]) %>%
    dplyr::mutate(precomp_file = factor(paste0(gRNA_precomp_dir, "/gRNA_precomp_", id, ".rds")))
  gRNA_dictionary_fp <- paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst")
  fst::write_fst(gRNA_dictionary, gRNA_dictionary_fp)

  # 3. gRNA-gene pairs
  pairs_dictionary <- create_dictionary(ids = gene_gRNA_group_pairs$gRNA_id, pod_size = pod_sizes[["pair"]]) %>%
    dplyr::rename(gRNA_id = id) %>% dplyr::mutate(gene_id = gene_gRNA_group_pairs$gene_id) %>% dplyr::select(gRNA_id, gene_id, pod_id) %>%
    dplyr::mutate(result_file = paste0(results_dir, "/result_", pod_id, ".fst") %>% factor())
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
run_gRNA_precomputation_at_scale <- function(pod_id, gRNA_precomp_dir, combined_perturbation_matrix, covariate_matrix, log_dir, B, seed) {
  # Activate the sink for the log file
  if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/gRNA_precomp_", pod_id, ".Rout"))
  # determine the gRNAs on which to run the precomputation
  gRNA_dictionary <- fst::read_fst(paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst")) %>% dplyr::filter(pod_id == !!pod_id)
  gRNA_ids <- gRNA_dictionary %>% dplyr::pull(id) %>% as.character()
  # run the precomputation for each of these gRNAs
  for (i in seq(1, length(gRNA_ids))) {
    gRNA_id <- gRNA_ids[i]
    cat(paste0("Running precomputation for gRNA ", gRNA_id, ".\n"))
    gRNA_indicators <- if (any(is(gene_matrix) %in% c("matrix", "Matrix"))) {
      combined_perturbation_matrix[gRNA_id,]
    } else {
      combined_perturbation_matrix[[gRNA_id,]] %>% as.numeric()
    }
    synth_data <- run_gRNA_precomputation(gRNA_indicators, covariate_matrix, B, seed)
    precomp_matrix_fp <- (gRNA_dictionary %>% dplyr::pull(precomp_file))[i] %>% as.character
    saveRDS(object = synth_data, precomp_matrix_fp)
  }
  if (!is.null(log_dir)) deactivate_sink()
}


#' Run_gene_precomputation_at_scale_round_1
#'
#' This function runs the first round of gene precomputations. In particular, it computes the raw dispersion estimate
#' and log geometric mean of each gene. It saves the results in the gene_precomp directory.
#'
#' @inheritParams run_sceptre_high_moi
#' @param pod_id pod id
#' @param gene_precomp_dir location of the gene precomputation directory
#' @param log_dir directory in which to sink the log file
run_gene_precomputation_at_scale_round_1 <- function(pod_id, gene_precomp_dir, gene_matrix, covariate_matrix, regularization_amount, log_dir) {
  if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/gene_precomp_round_1_pod_", pod_id, ".Rout"))

  gene_dictionary <- fst::read_fst(paste0(gene_precomp_dir, "/gene_dictionary.fst")) %>% dplyr::filter(pod_id == !!pod_id)
  gene_ids <- gene_dictionary %>% dplyr::pull(id) %>% as.character()
  precomps <- purrr::map(gene_ids, function(gene_id) {
    cat(paste0("Running precomputation round 1 for gene ", gene_id, ".\n"))
    expressions <- if (any(is(gene_matrix) %in% c("matrix", "Matrix"))) {
      gene_matrix[gene_id,]
    } else {
      gene_matrix[[gene_id,]] %>% as.numeric()
    }
    unreg_size <- run_gene_precomputation(expressions = expressions, covariate_matrix = covariate_matrix, gene_precomp_size = NULL)[["gene_precomp_size"]]
    out <- list()
    out$unreg_size <- unreg_size
    if (regularization_amount > 0) {
      gene_log_geom_mean <- log_geom_mean(expressions)
      out$gene_log_geom_mean <- gene_log_geom_mean
    }
    return(out)
  })

  names(precomps) <- gene_ids
  out_unreg_sizes <- purrr::map_dbl(precomps, function(l) l$unreg_size)
  if (regularization_amount > 0) out_log_geom_means <- purrr::map_dbl(precomps, function(l) l$gene_log_geom_mean)

  unreg_sizes_save_fp <- (gene_dictionary %>% dplyr::pull(size_unreg_file))[1] %>% as.character()
  saveRDS(object = out_unreg_sizes, file = unreg_sizes_save_fp)
  if (regularization_amount > 0) {
    geom_means_save_fp <- (gene_dictionary %>% dplyr::pull(geom_mean_file))[1] %>% as.character()
    saveRDS(object = out_log_geom_means, file = geom_means_save_fp)
  }

  if (!is.null(log_dir)) deactivate_sink()
}


#' Regularize genes at scale
#'
#' Regularizes the estimated gene sizes.
#'
#' @param gene_precomp_dir location of the gene precomputation directory
#' @param log_dir location of the directory in which to sink the log file
#' @param regularization_amount amount of regularization to apply to thetas (>= 0)
regularize_gene_sizes_at_scale <- function(gene_precomp_dir, regularization_amount, log_dir) {
  if (regularization_amount > 0) {
    if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/gene_precomp_regularize_sizes", ".Rout"))
    file_names <- list.files(gene_precomp_dir)
    gene_size_unreg_files <- paste0(gene_precomp_dir, "/", grep(pattern = 'gene_size_unreg_[0-9]+.rds', x = file_names, value = TRUE))
    gene_geom_mean_files <- paste0(gene_precomp_dir, "/", grep(pattern = 'gene_geom_mean_[0-9]+.rds', x = file_names, value = TRUE))
    load_distributed_vector <- function(f_names) f_names %>% purrr::map(readRDS) %>% purrr::reduce(c)
    gene_sizes_unreg <- load_distributed_vector(gene_size_unreg_files)
    gene_geom_means <- load_distributed_vector(gene_geom_mean_files)
    if (!all(names(gene_sizes_unreg) == names(gene_geom_means))) {
      gene_sizes_unreg <- gene_sizes_unreg[order(names(gene_sizes_unreg))]
      gene_geom_means <- gene_geom_means[order(names(gene_geom_means))]
    }
    cat("Regularizing gene sizes.\n")
    sizes_reg <- regularize_thetas(genes_log_gmean = gene_geom_means, theta = gene_sizes_unreg, plot_me = FALSE)
    names(sizes_reg) <- names(gene_sizes_unreg)
    sizes_reg_save_fp <- (fst::read_fst(paste0(gene_precomp_dir, "/gene_dictionary.fst"), "size_reg_file") %>% dplyr::pull())[1] %>% as.character()
    saveRDS(object = sizes_reg, file = sizes_reg_save_fp)
    if (!is.null(log_dir)) deactivate_sink()
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
run_gene_precomputation_at_scale_round_2 <- function(pod_id, gene_precomp_dir, gene_matrix, covariate_matrix, regularization_amount, log_dir) {
  if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/gene_precomp_round_2_pod_", pod_id, ".Rout"))

  gene_dictionary <- fst::read_fst(paste0(gene_precomp_dir, "/gene_dictionary.fst")) %>% dplyr::filter(pod_id == !!pod_id)
  gene_ids <- gene_dictionary %>% dplyr::pull(id) %>% as.character()
  if (regularization_amount > 0) {
    gene_sizes <- readRDS(as.character(gene_dictionary$size_reg_file[1]))[gene_ids]
  } else {
    gene_sizes <- readRDS(as.character(gene_dictionary$size_unreg_file[1]))
  }

  offsets <- sapply(gene_ids, function(gene_id) {
    cat(paste0("Running precomputation round 2 for gene ", gene_id, ".\n"))
    expressions <- if (any(is(gene_matrix) %in% c("matrix", "Matrix"))) {
      gene_matrix[gene_id,]
    } else {
      gene_matrix[[gene_id,]] %>% as.numeric()
    }
    dist_offsets <- run_gene_precomputation(expressions = expressions, covariate_matrix = covariate_matrix, gene_precomp_size = gene_sizes[[gene_id]])[["gene_precomp_offsets"]]
    return(dist_offsets)
  }) %>% as.data.frame()

  offsets_save_fp <- (gene_dictionary %>% dplyr::pull(offset_file))[1] %>% as.character()
  fst::write_fst(x = offsets, path = offsets_save_fp)
  if (!is.null(log_dir)) deactivate_sink()
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
run_gRNA_gene_pair_analysis_at_scale <- function(pod_id, gene_precomp_dir, gRNA_precomp_dir, results_dir, log_dir, gene_matrix, combined_perturbation_matrix, covariate_matrix, regularization_amount, side, B, full_output) {
  if (!is.null(log_dir)) activate_sink(paste0(log_dir, "/result_", pod_id, ".Rout"))

  results_dict <- fst::read_fst(paste0(results_dir, "/results_dictionary.fst")) %>% dplyr::filter(pod_id == !!pod_id)
  gene_dict <- fst::read_fst(paste0(gene_precomp_dir, "/gene_dictionary.fst"))
  gRNA_dict <- fst::read_fst(paste0(gRNA_precomp_dir, "/gRNA_dictionary.fst"))
  if (regularization_amount > 0) {
    regularized_gene_sizes <- readRDS(gene_dict$size_reg_file[1] %>% as.character())[as.character(results_dict$gene_id)]
  }
  out_l <- vector(mode = "list", length = nrow(results_dict))

  for (i in seq(1, nrow(results_dict))) {
    curr_gene <- results_dict[[i, "gene_id"]] %>% as.character()
    curr_gRNA <- results_dict[[i, "gRNA_id"]] %>% as.character()
    cat(paste0("Running distilled CRT on gene ", curr_gene, " and gRNA ", curr_gRNA, ".\n"))

    # Load the appropriate data from disk into memory
    # 1. gene
    if (i == 1 || results_dict[[i, "gene_id"]] != results_dict[[i - 1, "gene_id"]]) {
      gene_precomp_locs <- dplyr::filter(gene_dict, id == curr_gene)
      gene_offset_loc <- gene_precomp_locs %>% dplyr::pull(offset_file) %>% as.character
      gene_size_loc <- gene_precomp_locs %>% dplyr::pull(size_unreg_file) %>% as.character
      gene_precomp_offsets <- fst::read_fst(path = gene_offset_loc, columns = curr_gene) %>% dplyr::pull()
      if (regularization_amount > 0) {
        gene_precomp_size <- regularized_gene_sizes[[curr_gene]]
      } else {
        gene_precomp_size <- readRDS(file = gene_size_loc)[[curr_gene]]
      }
      expressions <- if (any(is(gene_matrix) %in% c("matrix", "Matrix"))) {
        gene_matrix[curr_gene,]
      } else {
        gene_matrix[[curr_gene,]] %>% as.numeric()
      }
    }

    # 2. gRNA
    if (i == 1 || results_dict[[i, "gRNA_id"]] != results_dict[[i - 1, "gRNA_id"]]) {
      gRNA_prcomp_loc <- dplyr::filter(gRNA_dict, id == curr_gRNA) %>% dplyr::pull(precomp_file) %>% as.character()
      gRNA_precomp <- readRDS(file = gRNA_prcomp_loc)
      gRNA_indicators <- if (any(is(combined_perturbation_matrix) %in% c("matrix", "Matrix"))) {
        combined_perturbation_matrix[curr_gRNA,]
      } else {
        combined_perturbation_matrix[[curr_gRNA,]] %>% as.numeric()
      }
    }

    # Run the dCRT
    out_l[[i]] <- run_sceptre_using_precomp_fast(expressions, gRNA_indicators, gRNA_precomp, side,
                                                 gene_precomp_size, gene_precomp_offsets, full_output)
  }
  # Create and save the result dataframe
  out_df <- data.table::rbindlist(out_l) %>%
    dplyr::mutate(gRNA_id = results_dict$gRNA_id, gene_id = results_dict$gene_id) %>%
    dplyr::relocate(gene_id, gRNA_id)
  out_fp <- (results_dict %>% dplyr::pull(result_file))[1] %>% as.character()
  fst::write_fst(out_df, out_fp)
  if (!is.null(log_dir)) deactivate_sink()
}


#' Collect results
#'
#' Collates the individual results files into a single result file.
#'
#' @param results_dir the directory containing the results
#' @param gene_gRNA_group_pairs the gene-gRNA groups data frame
#'
#' @return the collated results data frame.
collect_results <- function(results_dir, gene_gRNA_group_pairs) {
  file_names <- list.files(results_dir)
  to_load <- grep(pattern = 'result_[0-9]+.fst', x = file_names, value = TRUE)
  all_results <- results_dir %>% paste0("/", to_load) %>% purrr::map(fst::read_fst) %>%  data.table::rbindlist()
  if (ncol(gene_gRNA_group_pairs) >= 3) {
    all_results <- dplyr::left_join(gene_gRNA_group_pairs, all_results)
  }
  saveRDS(object = all_results, file = paste0(results_dir, "/all_results.rds"))
  return(all_results)
}
