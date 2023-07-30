construct_negative_control_pairs <- function(n_calibration_pairs, calibration_group_size, grna_assignments, response_matrix, grna_group_data_frame, response_grna_group_pairs, control_group_complement, low_moi) {
  # 1. set a few variables
  N_POSSIBLE_GROUPS_THRESHOLD <- 100L
  nt_grna_names <- names(grna_assignments[["indiv_nt_grna_idxs"]])
  n_nt_grnas <- length(nt_grna_names)

  # 2. compute nt_n_nonzero_matrix and n pairs
  if (is.na(n_calibration_pairs)) {
    grna_group_idxs <- grna_assignments$grna_group_idxs
    dt <- data.table::data.table(response_idx = match(x = response_grna_group_pairs$response_id,
                                                      table = rownames(response_matrix)),
                                 grna_idx = match(x = response_grna_group_pairs$grna_group,
                                                  table = names(grna_group_idxs))) |>
      data.table::setorder(cols = "response_idx")
    to_analyze_response_idxs = dt$response_idx
    to_analyze_grna_idxs = dt$grna_idx
  } else {
    to_analyze_response_idxs <- to_analyze_grna_idxs <- grna_group_idxs <- integer()
  }
  out <- compute_nt_nonzero_matrix_and_n_ok_pairs_v2(j = response_matrix@j,
                                                     p = response_matrix@p,
                                                     n_cells = ncol(response_matrix),
                                                     grna_group_idxs = grna_group_idxs,
                                                     indiv_nt_grna_idxs = grna_assignments$indiv_nt_grna_idxs,
                                                     all_nt_idxs = if (!control_group_complement) grna_assignments$all_nt_idxs else integer(),
                                                     to_analyze_response_idxs = to_analyze_response_idxs,
                                                     to_analyze_grna_idxs = to_analyze_grna_idxs,
                                                     compute_n_ok_pairs = is.na(n_calibration_pairs),
                                                     control_group_complement = control_group_complement)
  n_nonzero_m <- out$n_nonzero_mat
  n_nonzero_tot <- out$n_nonzero_tot

  # 3. compute the number of possible groups; compute the possible groups
  n_possible_groups <- choose(n_nt_grnas, calibration_group_size)
  if (n_possible_groups <= N_POSSIBLE_GROUPS_THRESHOLD) {
    # iterate over the entire set of combinations
    n_possible_groups <- as.integer(n_possible_groups)
    possible_groups_m <- iterate_over_combinations(n_nt_grnas, calibration_group_size, n_possible_groups)
  } else {
    # sample from the set of combinations
    possible_groups_m <- sample_combinations(calibration_group_size, n_calibration_pairs, n_nonzero_trt_thresh,
                                             n_nonzero_cntrl_thresh, n_possible_groups, n_nonzero_m, n_nonzero_tot,
                                             N_POSSIBLE_GROUPS_THRESHOLD)
  }

  # 5. sample WOR from the set of undercover pairs
  samp <- sample_undercover_pairs(n_nonzero_m, n_nonzero_tot, possible_groups_m, n_calibration_pairs,
                                  n_nonzero_trt_thresh, n_nonzero_cntrl_thresh, low_moi,
                                  j = response_matrix@j, p = response_matrix@p, n_cells = ncol(response_matrix),
                                  n_genes = nrow(response_matrix), indiv_nt_grna_idxs = grna_assignments$indiv_nt_grna_idxs)

  # 6. construct the data frame of negative control pairs
  response_ids <- factor(rownames(response_matrix))
  increment_matrix(possible_groups_m)
  undercover_groups <- factor(get_undercover_group_names(possible_groups_m, nt_grna_names))
  df <- data.frame(response_id = response_ids[samp$response_idxs],
                   grna_group = undercover_groups[samp$grna_group_idxs],
                   n_nonzero_trt = samp$n_nonzero_trt_v,
                   n_nonzero_cntrl = samp$n_nonzero_cntrl_v)

  # 7. check the number of rows of df, and issue a warning if less than the required amount
  if (nrow(df) < n_calibration_pairs) {
    warning(paste0("Unable to generate ", n_calibration_pairs, " negative control pairs (i.e., the number of discovery pairs that passes pairwise QC). Consider increasing `calibration_group_size`. If possible, consider rerunning the experiment with more negative control gRNAs."))
  }

  return(df)
}


compute_calibration_group_size <- function(grna_group_data_frame) {
    median_group_size <- grna_group_data_frame |>
      dplyr::filter(grna_group != "non-targeting") |>
      dplyr::pull(grna_group) |> table() |> stats::median() |> round()
    n_ntc <- grna_group_data_frame |> dplyr::filter(grna_group == "non-targeting") |> nrow()
    calibration_group_size <- min(median_group_size, floor(n_ntc/2))
    return(calibration_group_size)
}


get_undercover_group_names <- function(possible_groups_m, nt_grna_names) {
  apply(X = possible_groups_m, MARGIN = 1, FUN = function(r) {
    paste0(nt_grna_names[r], collapse = "&")
  })
}


