construct_negative_control_pairs <- function(n_calibration_pairs, calibration_group_size,
                                             grna_assignments, response_matrix,
                                             grna_group_data_frame, low_moi,
                                             n_nonzero_trt_thresh, n_nonzero_cntrl_thresh,
                                             n_nonzero_m, n_nonzero_tot) {
  # 1. set a few variables
  N_POSSIBLE_GROUPS_THRESHOLD <- 100L
  nt_grna_names <- names(grna_assignments[["indiv_nt_grna_idxs"]])
  n_nt_grnas <- length(nt_grna_names)
  if (is.na(n_calibration_pairs)) n_calibration_pairs <- sceptre_object@n_ok_discovery_pairs

  # 2. compute the number of possible groups; compute the possible groups
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

  # 3. sample WOR from the set of undercover pairs
  samp <- sample_undercover_pairs(n_nonzero_m, n_nonzero_tot, possible_groups_m, n_calibration_pairs,
                                  n_nonzero_trt_thresh, n_nonzero_cntrl_thresh, low_moi,
                                  j = response_matrix@j, p = response_matrix@p, n_cells = ncol(response_matrix),
                                  n_genes = nrow(response_matrix), indiv_nt_grna_idxs = grna_assignments$indiv_nt_grna_idxs)

  # 4. construct the data frame of negative control pairs
  response_ids <- factor(rownames(response_matrix))
  increment_matrix(possible_groups_m)
  undercover_groups <- factor(get_undercover_group_names(possible_groups_m, nt_grna_names))
  df <- data.frame(response_id = response_ids[samp$response_idxs],
                   grna_group = undercover_groups[samp$grna_group_idxs],
                   n_nonzero_trt = samp$n_nonzero_trt_v,
                   n_nonzero_cntrl = samp$n_nonzero_cntrl_v)

  # 5. check the number of rows of df, and issue a warning if less than the required amount
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
