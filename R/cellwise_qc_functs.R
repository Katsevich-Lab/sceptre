determine_cells_to_retain <- function(sceptre_object, response_n_umis_range, response_n_nonzero_range, p_mito_threshold, additional_cells_to_remove) {
  # 1. compute the cells to retain based on n_umis and n_nonzero
  exclude_cells_by_clipping_counts <- function(v, percentile_range) {
    cutoffs <- stats::quantile(x = v, probs = percentile_range)
    which(v < cutoffs[1] | v > cutoffs[2])
  }
  cells_to_exclude_n_umis <- exclude_cells_by_clipping_counts(sceptre_object@covariate_data_frame$response_n_umis, response_n_umis_range)
  n_cells_rm_n_umis <- length(cells_to_exclude_n_umis)
  cells_to_exclude_n_responses <- exclude_cells_by_clipping_counts(sceptre_object@covariate_data_frame$response_n_nonzero, response_n_nonzero_range)
  n_cells_rm_n_responses <- length(cells_to_exclude_n_responses)

  # 2. compute cells to retain based on p_mito
  if ("response_p_mito" %in% colnames(sceptre_object@covariate_data_frame)) {
    cells_to_exclude_p_mito <- which(sceptre_object@covariate_data_frame$response_p_mito > p_mito_threshold)
    n_cells_rm_p_mito <- length(cells_to_exclude_p_mito)
  } else {
    cells_to_exclude_p_mito <- integer()
    n_cells_rm_p_mito <- 0L
  }

  # 3. compute cells to retain based on cells containing multiple grnas (if in low MOI)
  cells_to_exclude_multiple_grnas <- sceptre_object@cells_w_multiple_grnas
  n_cells_rm_multiple_grnas <- length(sceptre_object@cells_w_multiple_grnas)

  # 4. remove additional cells specified by the user
  cells_to_exclude_user_specified <- additional_cells_to_remove
  n_cells_rm_user_specified <- length(cells_to_exclude_user_specified)

  # 5. finally, determine the set of cells to retain, and update the sceptre_object
  cells_to_exclude <- c(cells_to_exclude_n_umis, cells_to_exclude_n_responses,
                        cells_to_exclude_p_mito, cells_to_exclude_multiple_grnas,
                        cells_to_exclude_user_specified) |> unique()
  n_cells <- ncol(sceptre_object@response_matrix)
  cells_to_retain <- if (length(cells_to_exclude) == 0L) seq(1L, n_cells) else seq(1L, n_cells)[-cells_to_exclude]
  sceptre_object@cells_in_use <- cells_to_retain
  n_cells_rm_total <- n_cells - length(cells_to_retain)
  cell_removal_metrics <- c(n_cells_rm_n_umis = n_cells_rm_n_umis,
                            n_cells_rm_n_responses = n_cells_rm_n_responses,
                            n_cells_rm_p_mito = n_cells_rm_p_mito,
                            n_cells_rm_multiple_grnas = n_cells_rm_multiple_grnas,
                            n_cells_rm_user_specified = n_cells_rm_user_specified,
                            n_cells_rm_total = n_cells_rm_total)
  sceptre_object@cell_removal_metrics <- cell_removal_metrics
  return(sceptre_object)
}


update_grna_assignments_given_qc <- function(sceptre_object) {
  # 0. define several varaibles
  cells_in_use <- sceptre_object@cells_in_use
  grna_assignments_raw <- sceptre_object@grna_assignments_raw
  n_cells <- ncol(sceptre_object@response_matrix)

  # 1. define update idxs funct
  update_idxs <- function(v, cells_in_use, n_cells) {
    v_log <- logical(length = n_cells)
    v_log[v] <- TRUE
    v_log_sub <- v_log[cells_in_use] # subset step
    v_updated <- which(v_log_sub) # determine the positions of the nonzero entries
    return(v_updated)
  }

  # 2. for each gRNA group and NT gRNA, subset the grna vector and the update the indices
  grna_group_idxs_new <- sapply(grna_assignments_raw$grna_group_idxs, function(v) update_idxs(v, cells_in_use, n_cells), simplify = FALSE)
  nt_idxs_new <- sapply(grna_assignments_raw$indiv_nt_grna_idxs, function(v) update_idxs(v, cells_in_use, n_cells), simplify = FALSE)
  # remove those nt grnas with 0 cells (after QC)
  nt_idxs_new <- nt_idxs_new[sapply(nt_idxs_new, length) != 0L]
  grna_assignments <- list(grna_group_idxs = grna_group_idxs_new, indiv_nt_grna_idxs = nt_idxs_new)

  # 3. if using the NT cells, update indiv gRNA indices so that they are relative to all NTs
  if (!sceptre_object@control_group_complement) {
    l <- update_indiv_grna_assignments_for_nt_cells(grna_assignments$indiv_nt_grna_idxs)
    grna_assignments$indiv_nt_grna_idxs <- l$indiv_nt_grna_idxs
    grna_assignments$all_nt_idxs <- l$all_nt_idxs
  }

  # 4. update the sceptre_object and return
  sceptre_object@grna_assignments <- grna_assignments
  return(sceptre_object)
}


update_indiv_grna_assignments_for_nt_cells <- function(indiv_nt_grna_idxs) {
  out <- list()
  nt_grnas <- names(indiv_nt_grna_idxs)
  all_nt_idxs <- unique(stats::setNames(unlist(indiv_nt_grna_idxs), NULL))
  n_cells_per_nt <- sapply(indiv_nt_grna_idxs, length)
  stop <- cumsum(n_cells_per_nt)
  start <- c(0L, stop[-length(stop)]) + 1L
  indiv_nt_grna_idxs <- lapply(seq(1, length(nt_grnas)), function(i) {
    seq(start[i], stop[i])
  }) |> stats::setNames(nt_grnas)
  out$indiv_nt_grna_idxs <- indiv_nt_grna_idxs
  out$all_nt_idxs <- all_nt_idxs
  return(out)
}
