# The possible fields in the output are grna_group_idxs, all_nt_idxs, indiv_nt_grna_idxs.
# If in high moi, we always return the gRNA-group-to-idx map (grna_group_idxs) for the non-targeting gRNA groups; if we are running a calibration check, we additionally return the individual NT gRNA idxs (indiv_nt_grna_idxs). This latter vector is absolute (i.e., not relative to any other vector).
# If in low moi, we likewise always return the gRNA-group-to-idx map for the non-targeting gRNA groups. If the control group is the NT cells, we additionally return all_nt_idxs, which is the set of NT cell idxs. Finally, if we are running a calibration check, we return the indices of the individual NT gRNAs. If the control group is the NT cells, then these indices are relative to the NT cells. If the control group is the complement set, then these indices are absolute.
assign_grnas_to_cells <- function(sceptre_object) {
  # extract pieces from sceptre_object
  grna_matrix <- sceptre_object@grna_matrix
  grna_group_data_frame <- sceptre_object@grna_group_data_frame
  low_moi <- sceptre_object@low_moi
  grna_assignment_method <- sceptre_object@grna_assignment_method
  cell_covariate_data_frame <- sceptre_object@covariate_data_frame
  grna_assignment_hyperparameters <- sceptre_object@grna_assignment_hyperparameters
  n_cells <- ncol(grna_matrix)

  # assign grnas via the selected strategy; obtain the grna assignments and cells containing multiple grnas
  if (grna_assignment_method == "mixture") {
    initial_assignment_list <- assign_grnas_to_cells_mixture(grna_matrix = grna_matrix,
                                                             cell_covariate_data_frame = cell_covariate_data_frame,
                                                             grna_assignment_hyperparameters = grna_assignment_hyperparameters)
    run_process_initial_assignment_list <- TRUE
  }
  if (grna_assignment_method == "thresholding") {
    initial_assignment_list <- assign_grnas_to_cells_thresholding(grna_matrix = grna_matrix,
                                                                  grna_assign_threshold = grna_assignment_hyperparameters$threshold)
    run_process_initial_assignment_list <- TRUE
  }
  if (grna_assignment_method == "maximum") {
    out_list <- assign_grnas_to_cells_maximum(grna_matrix = grna_matrix,
                                              grna_group_data_frame = grna_group_data_frame,
                                              grna_lib_size = cell_covariate_data_frame$grna_n_umis)
    grna_assignments <- out_list$grna_assignments
    cells_w_multiple_grnas <- which(out_list$frac_umis < sceptre_object@grna_assignment_hyperparameters$umi_fraction_threshold)
  }

  # process the initial assignment list
  if (run_process_initial_assignment_list) {
    processed_assignment_out <- process_initial_assignment_list(initial_assignment_list = initial_assignment_list,
                                                                grna_group_data_frame = grna_group_data_frame,
                                                                n_cells = n_cells, low_moi = low_moi)
    sceptre_object@grna_assignments_raw <- processed_assignment_out$grna_assignments_raw
    sceptre_object@cells_w_multiple_grnas <- processed_assignment_out$cells_w_multiple_grnas

  }

  # update sceptre object with grna_assignments and cells_w_multiple_grnas
  #sceptre_object@grna_assignments_raw <- grna_assignments
  #sceptre_object@cells_w_multiple_grnas <- cells_w_multiple_grnas
  return(sceptre_object)
}


process_initial_assignment_list <- function(initial_assignment_list, grna_group_data_frame, n_cells, low_moi) {
  # 0. compute the number of cells per grna
  cells_per_grna <- sapply(initial_assignment_list, length)
  # 1. compute the vector of grnas per cell
  grnas_per_cell <- compute_n_grnas_per_cell_vector(initial_assignment_list, n_cells)
  # 2. determine the cells that contain multiple grnas (if in low MOI)
  cells_w_multiple_grnas <- if (low_moi) which(grnas_per_cell >= 1L) else integer()
  # 3. pool together targeting gRNAs via the or operation
  targeting_grna_group_data_table <- grna_group_data_frame |>
    dplyr::filter(grna_group != "non-targeting") |> data.table::as.data.table()
  targeting_grna_groups <- targeting_grna_group_data_table$grna_group |> unique()
  grna_group_idxs <- lapply(targeting_grna_groups, function(targeting_grna_group) {
    curr_grna_ids <- targeting_grna_group_data_table[
      targeting_grna_group_data_table$grna_group == targeting_grna_group,]$grna_id
    initial_assignment_list[curr_grna_ids] |> unlist() |> unique()
  }) |> stats::setNames(targeting_grna_groups)
  # 4. obtain the individual non-targeting grna idxs
  nontargeting_grna_ids <- grna_group_data_frame |>
    dplyr::filter(grna_group == "non-targeting") |> dplyr::pull(grna_id)
  indiv_nt_grna_idxs <- initial_assignment_list[nontargeting_grna_ids]
  # 5. construct the grna_group_idxs list
  grna_assignments_raw <- list(grna_group_idxs = grna_group_idxs,
                               indiv_nt_grna_idxs = indiv_nt_grna_idxs)
  # 6. compute the number of cells per targeting grna group
  cells_per_targeting_grna_group <- sapply(grna_group_idxs, length)
  # 7. compute the number of targeting grna groups per cell
  out <- list(grna_assignments_raw = grna_assignments_raw, cells_per_grna = cells_per_grna,
              grnas_per_cell = grnas_per_cell, cells_w_multiple_grnas = cells_w_multiple_grnas,
              cells_per_targeting_grna_group = cells_per_targeting_grna_group)
  return(out)
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
