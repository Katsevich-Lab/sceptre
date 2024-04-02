set_matrix_accessibility <- function(matrix_in, make_row_accessible = TRUE) {
  # if a logical matrix, convert to the corresponding numeric matrix; consider more efficient implementation later
  if (methods::is(matrix_in, "lgRMatrix")) {
    attr(matrix_in, "class") <- "dgRMatrix"
    matrix_in@x <- rep(1.0, length(matrix_in@j))
  }
  if (methods::is(matrix_in, "lgCMatrix")) {
    attr(matrix_in, "class") <- "dgCMatrix"
    matrix_in@x <- rep(1.0, length(matrix_in@i))
  }
  if (methods::is(matrix_in, "lgTMatrix")) {
    attr(matrix_in, "class") <- "dgTMatrix"
    matrix_in@x <- rep(1.0, length(matrix_in@j))
  }

  if (methods::is(matrix_in, "dgRMatrix") && make_row_accessible) {
    out <- matrix_in
  } else if (methods::is(matrix_in, "dgCMatrix") && !make_row_accessible) {
    out <- matrix_in
  } else if (methods::is(matrix_in, "odm")) {
    out <- matrix_in
  } else {
    # basic info: get matrix_in, n_responses, response_ids, and n_cells
    n_responses <- nrow(matrix_in)
    response_ids <- rownames(matrix_in)
    n_cells <- ncol(matrix_in)
    if (methods::is(matrix_in, "dgTMatrix")) {
      i <- matrix_in@i
      j <- matrix_in@j
    } else if (methods::is(matrix_in, "matrix")) {
      matrix_in <- methods::as(matrix_in, "TsparseMatrix")
      i <- matrix_in@i
      j <- matrix_in@j
    } else if (methods::is(matrix_in, "dgCMatrix")) {
      i <- matrix_in@i
      j <- convert_pointer_to_index_vector_v2(matrix_in@p)
    } else if (methods::is(matrix_in, "dgRMatrix")) {
      i <- convert_pointer_to_index_vector_v2(matrix_in@p)
      j <- matrix_in@j
    } else {
      stop("Class not recognized.")
    }
    x <- matrix_in@x

    # perform the sorting operation
    sort_in_place <- methods::is(matrix_in, "matrix") || methods::is(matrix_in, "dgTMatrix")
    if (sort_in_place) {
      l <- list(i = i, j = j, x = x)
      dt <- data.table::setDT(l) # assign by reference
    } else {
      dt <- data.table::data.table(i = i, j = j, x = x) # deep copy
    }

    # finally, initialize the output
    if (make_row_accessible) {
      data.table::setorderv(x = dt, cols = c("i", "j"))
      p <- obtain_pointer_vector(i = dt$i, dim = n_responses)
      out <- Matrix::sparseMatrix(j = 1, p = c(0, 1), x = 1, repr = "R")
      out@j <- dt$j
    } else {
      data.table::setorderv(x = dt, cols = c("j", "i"))
      p <- obtain_pointer_vector(i = dt$j, dim = n_cells)
      out <- Matrix::sparseMatrix(i = 1, p = c(0, 1), x = 1, repr = "C")
      out@i <- dt$i
    }
    out@p <- p
    out@x <- dt$x
    out@Dim <- c(n_responses, n_cells)
    rownames(out) <- response_ids
  }
  return(out)
}


convert_pointer_to_index_vector_v2 <- function(p) {
  l <- diff(p)
  rep(x = seq(0, length(l) - 1L), times = l)
}

convert_covariate_df_to_design_matrix <- function(covariate_data_frame, formula_object) {
  global_cell_covariates_new <- stats::model.matrix(object = formula_object, data = covariate_data_frame)
  # verify that the global cell covariates are OK after transformation
  if (nrow(global_cell_covariates_new) != nrow(covariate_data_frame)) {
    stop("Some rows of `covariate_data_frame` were lost after applying the `formula object`. This is likely due to NAs, which should be removed first.")
  }
  for (col_name in colnames(global_cell_covariates_new)) {
    vect <- global_cell_covariates_new[, col_name]
    if (any(is.infinite(vect)) || any(is.na(vect))) {
      stop("The column `", col_name, "` of the `covariate_data_frame` after the `formula object` has been applied contains entries that are -Inf, Inf, or NA. Remove these entries.")
    }
  }
  # verify that matrix is not rank-deficient
  rank_matrix <- Matrix::rankMatrix(global_cell_covariates_new)
  if (rank_matrix != ncol(global_cell_covariates_new)) {
    stop("The `formula_object` contains redundant information. This often occurs when one variable is `nested` inside another. For example, if the data frame contains a column `lane` with levels `lane_1`, `lane_2`, `lane_3`, and `lane_4` and a column `batch` with levels `batch_1` and `batch_2`, and if `lane_1` and `lane_2` are contained entirely within `batch_1` while `lane_3` and `lane_4` are contained entirely within `batch`2`, then `batch` is redundant. In this case `batch` should be removed from the formula object.")
  }
  return(global_cell_covariates_new)
}


compute_cell_covariates <- function(matrix_in, feature_names, compute_p_mito, compute_max_feature) {
  # make response matrix column accessible
  matrix_in <- set_matrix_accessibility(matrix_in, make_row_accessible = FALSE)
  # get MT gene idxs
  if (compute_p_mito) {
    mt_gene_idxs <- grep(pattern = "^MT-", x = feature_names, ignore.case = TRUE)
    compute_p_mito <- length(mt_gene_idxs) >= 1
  } else {
    mt_gene_idxs <- integer()
  }
  # call the low-level function
  out <- compute_cell_covariates_cpp(
    i = matrix_in@i,
    p = matrix_in@p,
    x = matrix_in@x,
    n_genes = nrow(matrix_in),
    n_cells = ncol(matrix_in),
    mt_gene_idxs = mt_gene_idxs,
    compute_p_mito = compute_p_mito,
    compute_max_feature = compute_max_feature
  )
  ret <- data.frame(n_nonzero = out$n_nonzero, n_umis = out$n_umi)
  if (compute_p_mito) ret$p_mito <- out$p_mito
  if (compute_max_feature) {
    ret$frac_umis_max_feature <- out$frac_umis_max_feature
    ret$feature_w_max_expression <- rownames(matrix_in)[out$max_feature + 1L]
  }
  return(ret)
}


get_synthetic_permutation_idxs <- function(grna_assignments, B, calibration_check, control_group_complement, calibration_group_size, n_cells) {
  if (calibration_check && !control_group_complement) { # 1. calibration check, nt cells (low MOI only)
    indiv_nt_sizes <- vapply(grna_assignments$indiv_nt_grna_idxs, length, FUN.VALUE = integer(1)) |> sort(decreasing = TRUE)
    M <- sum(indiv_nt_sizes[seq(1, calibration_group_size)])
    n_control_cells <- length(grna_assignments$all_nt_idxs)
    out <- fisher_yates_samlper(n_tot = n_control_cells, M = M, B = B)
  } else if (calibration_check && control_group_complement) { # 2. calibration check, complement (low and high MOI)
    indiv_nt_sizes <- vapply(grna_assignments$indiv_nt_grna_idxs, length, FUN.VALUE = integer(1)) |> sort(decreasing = TRUE)
    M <- sum(indiv_nt_sizes[seq(1, calibration_group_size)])
    out <- fisher_yates_samlper(n_tot = n_cells, M = M, B = B)
  } else if (!calibration_check && !control_group_complement) { # 3. discovery, nt cells (low MOI only)
    grna_group_sizes <- vapply(grna_assignments$grna_group_idxs, length, FUN.VALUE = integer(1))
    grna_group_sizes <- grna_group_sizes[grna_group_sizes != 0L]
    range_grna_group_sizes <- range(grna_group_sizes)
    n_control_cells <- length(grna_assignments$all_nt_idxs)
    out <- hybrid_fisher_iwor_sampler(
      N = n_control_cells,
      m = range_grna_group_sizes[1],
      M = range_grna_group_sizes[2],
      B = B
    )
  } else if (!calibration_check && control_group_complement) { # 4. discovery, complement (low and high MOI)
    max_cells_per_grna_group <- vapply(grna_assignments$grna_group_idxs, length, FUN.VALUE = integer(1)) |> max()
    out <- fisher_yates_samlper(n_tot = n_cells, M = max_cells_per_grna_group, B = B)
  }

  return(out)
}


update_dfs_based_on_grouping_strategy <- function(sceptre_object) {
  grouping_strategy <- sceptre_object@grna_integration_strategy
  if (grouping_strategy == "union") {
    sceptre_object@grna_target_data_frame$grna_group <- sceptre_object@grna_target_data_frame$grna_target
    sceptre_object@discovery_pairs$grna_group <- sceptre_object@discovery_pairs$grna_target
    sceptre_object@positive_control_pairs$grna_group <- sceptre_object@positive_control_pairs$grna_target
  }
  if (grouping_strategy %in% c("singleton", "bonferroni")) {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame
    discovery_pairs <- sceptre_object@discovery_pairs
    positive_control_pairs <- sceptre_object@positive_control_pairs

    # function to update targeting pairs
    update_targeting_pairs <- function(targeting_pairs, updated_grna_target_data_frame) {
      dplyr::left_join(targeting_pairs,
        updated_grna_target_data_frame |> dplyr::select(grna_target, grna_group),
        by = "grna_target", relationship = "many-to-many"
      )
    }

    # update each data frame, appending a grna group column
    updated_grna_target_data_frame <- grna_target_data_frame |> dplyr::mutate(grna_group = ifelse(grna_target == "non-targeting", "non-targeting", grna_id))
    sceptre_object@grna_target_data_frame <- updated_grna_target_data_frame
    sceptre_object@discovery_pairs <- update_targeting_pairs(discovery_pairs, updated_grna_target_data_frame)
    sceptre_object@positive_control_pairs <- update_targeting_pairs(positive_control_pairs, updated_grna_target_data_frame)
  }
  return(sceptre_object)
}
