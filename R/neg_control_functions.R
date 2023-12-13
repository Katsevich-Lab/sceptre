construct_negative_control_pairs_v2 <- function(sceptre_object, n_calibration_pairs, calibration_group_size) {
  cat("Constructing negative control pairs.")
  grna_assignments <- sceptre_object@grna_assignments
  response_matrix <- sceptre_object@response_matrix
  grna_target_data_frame <- sceptre_object@grna_target_data_frame
  low_moi <- sceptre_object@low_moi
  n_nonzero_trt_thresh <- sceptre_object@n_nonzero_trt_thresh
  n_nonzero_cntrl_thresh <- sceptre_object@n_nonzero_cntrl_thresh
  n_nonzero_m <- sceptre_object@M_matrix
  n_nonzero_tot <- sceptre_object@n_nonzero_tot_vector
  cells_in_use <- sceptre_object@cells_in_use
  discovery_pairs_with_info <- sceptre_object@discovery_pairs_with_info

  # 1. set a few variables
  N_POSSIBLE_GROUPS_THRESHOLD <- 100L
  nt_grna_names <- names(grna_assignments[["indiv_nt_grna_idxs"]])
  n_nt_grnas <- length(nt_grna_names)
  n_genes <- nrow(response_matrix)

  # 2. compute the number of possible groups; compute the possible groups
  n_possible_groups <- choose(n_nt_grnas, calibration_group_size)
  if (n_possible_groups <= N_POSSIBLE_GROUPS_THRESHOLD) {
    # iterate over the entire set of combinations
    n_possible_groups <- as.integer(n_possible_groups)
    possible_groups_m <- iterate_over_combinations(n_nt_grnas, calibration_group_size, n_possible_groups)
  } else {
    # sample from the set of combinations
    p_hat <- mean(discovery_pairs_with_info$pass_qc)
    possible_groups_m <- sample_combinations_v2(calibration_group_size, n_calibration_pairs, n_possible_groups,
                                                n_nt_grnas, as.numeric(n_genes), N_POSSIBLE_GROUPS_THRESHOLD, p_hat)
  }

  # 3. sample WOR from the set of undercover pairs
  calculate_ess_using_m_matrix <- sceptre_object@low_moi || calibration_group_size == 1L
  if (calculate_ess_using_m_matrix) {
    samp <- sample_undercover_pairs_v2(n_nonzero_m = n_nonzero_m, n_nonzero_tot = n_nonzero_tot,
                                       possible_groups_m = possible_groups_m, n_genes = n_genes,
                                       n_calibration_pairs = n_calibration_pairs,
                                       n_nonzero_trt_thresh = n_nonzero_trt_thresh,
                                       n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh)
  } else {
    nt_grna_group_idxs <- lapply(X = seq(1L, nrow(possible_groups_m)), function(i) {
      nt_grna_idxs <- possible_groups_m[i,] + 1L
      grna_assignments$indiv_nt_grna_idxs[nt_grna_idxs] |> unlist() |> unique()
    })
    if (methods::is(response_matrix, "odm")) {
      n_nonzero_m <- ondisc:::compute_n_trt_cells_matrix_ondisc(file_name_in = response_matrix@h5_file,
                                                                f_row_ptr = response_matrix@ptr,
                                                                n_cells_orig = ncol(response_matrix),
                                                                n_cells_sub = length(cells_in_use),
                                                                n_genes = n_genes,
                                                                nt_grna_group_idxs = nt_grna_group_idxs,
                                                                cells_in_use = cells_in_use)
    } else {
      n_nonzero_m <- compute_n_trt_cells_matrix(j = response_matrix@j, p = response_matrix@p,
                                                n_cells_orig = ncol(response_matrix),
                                                n_cells_sub = length(cells_in_use),
                                                n_genes = n_genes,
                                                nt_grna_group_idxs = nt_grna_group_idxs,
                                                cells_in_use = cells_in_use)
    }
    dummy_possibly_groups_m <- matrix(seq_along(nt_grna_group_idxs) - 1L, ncol = 1L)
    samp <- sample_undercover_pairs_v2(n_nonzero_m = n_nonzero_m, n_nonzero_tot = n_nonzero_tot,
                                       possible_groups_m = dummy_possibly_groups_m, n_genes = n_genes,
                                       n_calibration_pairs = n_calibration_pairs,
                                       n_nonzero_trt_thresh = n_nonzero_trt_thresh,
                                       n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh)
  }

  # 4. construct the data frame of negative control pairs
  if (length(samp$response_idxs) == 0L) {
    stop("No negative control pair passes pairwise QC.")
  }
  response_ids <- rownames(response_matrix)
  increment_matrix(possible_groups_m)
  undercover_groups <- get_undercover_group_names(possible_groups_m, nt_grna_names)
  df <- data.frame(response_id = response_ids[samp$response_idxs],
                   grna_group = undercover_groups[samp$grna_group_idxs],
                   n_nonzero_trt = samp$n_nonzero_trt_v,
                   n_nonzero_cntrl = samp$n_nonzero_cntrl_v,
                   pass_qc = TRUE)
  cat(crayon::green(' \u2713\n'))

  # 5. check the number of rows of df, and issue a warning if less than the required amount
  if (nrow(df) < n_calibration_pairs) {
    cat(crayon::red(paste0("Note: Unable to generate the number of negative control pairs (", n_calibration_pairs, ") requested. Generating as many negative control pairs (", nrow(df),") as possible.\n")))
  }
  return(df)
}


compute_calibration_group_size <- function(grna_target_data_frame) {
  median_group_size <- grna_target_data_frame |>
    dplyr::filter(grna_group != "non-targeting") |>
    dplyr::pull(grna_group) |> table() |> stats::median() |> round()
  n_ntc <- grna_target_data_frame |> dplyr::filter(grna_group == "non-targeting") |> nrow()
  calibration_group_size <- as.integer(min(median_group_size, n_ntc))
  return(calibration_group_size)
}


get_undercover_group_names <- function(possible_groups_m, nt_grna_names) {
  apply(X = possible_groups_m, MARGIN = 1, FUN = function(r) {
    paste0(nt_grna_names[r], collapse = "&")
  })
}
