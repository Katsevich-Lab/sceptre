
assign_grnas_to_cells_thresholding <- function(grna_matrix, grna_assign_threshold, grna_group_data_frame) {
  out <- list()

  # 1. make the matrix row-accessible
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
  grna_ids <- rownames(grna_matrix)

  # 2. obtain the grna ids and groups
  grna_groups <- as.character(unique(grna_group_data_frame$grna_group))
  grna_groups <- grna_groups[grna_groups != "non-targeting"]

  # 3. loop over the targeting grna groups, obtaining the cell assignments for each
  grna_group_idxs <- sapply(grna_groups, function(grna_group) {
    l <- grna_group_data_frame$grna_group == grna_group
    curr_grna_ids <- grna_group_data_frame$grna_id[l]
    row_idxs <- match(x = curr_grna_ids, grna_ids)
    cell_idxs <- group_and_threshold(j = grna_matrix@j, p = grna_matrix@p, x = grna_matrix@x,
                                     row_idxs = row_idxs, threshold = grna_assign_threshold)
    return(cell_idxs)
  })
  out$grna_group_idxs <- grna_group_idxs


  # 4. when running a calibration check, also the individual NT gRNAs
  nt_grnas <- grna_group_data_frame |>
    dplyr::filter(grna_group == "non-targeting") |> dplyr::pull("grna_id")
  indiv_nt_grna_idxs <- sapply(nt_grnas, function(nt_grna) {
    row_idx <- match(x = nt_grna, grna_ids)
    cell_idxs <- group_and_threshold(j = grna_matrix@j, p = grna_matrix@p, x = grna_matrix@x,
                                     row_idxs = row_idx, threshold = grna_assign_threshold)
    return(cell_idxs)
  })
  out$indiv_nt_grna_idxs <- indiv_nt_grna_idxs
  return(out)
}


assign_grnas_to_cells_maximum <- function(grna_matrix, grna_group_data_frame, grna_lib_size) {
  grna_assignments <- list()

  # 1. make grna matrix column accessible
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = FALSE)

  # 2. get the individual grna assignments
  l <- compute_colwise_max(i = grna_matrix@i, p = grna_matrix@p,
                           x = grna_matrix@x, n_cells = ncol(grna_matrix),
                           grna_lib_size = grna_lib_size)
  grna_idx <- l$assignment_vect
  grna_rownames <- factor(rownames(grna_matrix))
  indiv_grna_id_assignments <- grna_rownames[grna_idx]

  # 3. obtain the grna group assignments
  grna_group_assignments <- grna_group_data_frame$grna_group[match(x = indiv_grna_id_assignments,
                                                                   table = grna_group_data_frame$grna_id)]

  # 4. convert the grna group assignments into a map
  unique_grna_groups <- unique(grna_group_data_frame$grna_group)
  unique_grna_groups <- unique_grna_groups[unique_grna_groups != "non-targeting"]
  grna_group_idxs <- lapply(unique_grna_groups, function(unique_grna_group) {
    which(grna_group_assignments == unique_grna_group)
  }) |> stats::setNames(unique_grna_groups)
  grna_assignments$grna_group_idxs <- grna_group_idxs

  # 5. return the indices of the individual nt grnas
  nt_grnas <- grna_group_data_frame |>
    dplyr::filter(grna_group == "non-targeting") |> dplyr::pull("grna_id")
  indiv_nt_grna_idxs <- lapply(nt_grnas, function(nt_grna) {
    which(nt_grna == indiv_grna_id_assignments)
  }) |> stats::setNames(nt_grnas)
  grna_assignments$indiv_nt_grna_idxs <- indiv_nt_grna_idxs

  return(list(grna_assignments = grna_assignments, frac_umis = l$frac_umis))
}
