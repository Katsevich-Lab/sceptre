construct_negative_control_grna_groups <- function(sceptre_object) {
  # extract relevant fields from sceptre_object
  grna_assignments <- sceptre_object@grna_assignments
  grna_target_data_frame <- sceptre_object@grna_target_data_frame
  discovery_pairs <- sceptre_object@discovery_pairs

  # 1. set a few variables
  N_POSSIBLE_GROUPS_THRESHOLD <- 100L
  nt_grna_names <- names(grna_assignments[["indiv_nt_grna_idxs"]])
  n_nt_grnas <- length(nt_grna_names)

  # 2. determine the group size
  calibration_group_size <- compute_calibration_group_size(grna_target_data_frame)

  # 3. compute the number of possible groups; compute the possible groups
  n_possible_groups <- choose(n_nt_grnas, calibration_group_size)
  if (n_possible_groups <= N_POSSIBLE_GROUPS_THRESHOLD) {
    # iterate over the entire set of combinations
    n_possible_groups <- as.integer(n_possible_groups)
    possible_groups_m <- iterate_over_combinations(n_nt_grnas, calibration_group_size, n_possible_groups)
  } else {
    # sample from the set of combinations
    possible_groups_m <- sample_combinations_v2(calibration_group_size, n_calibration_pairs, n_possible_groups,
                                                n_nt_grnas, as.numeric(n_genes), N_POSSIBLE_GROUPS_THRESHOLD, p_hat)
  }
}
