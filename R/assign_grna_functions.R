# The possible fields in the output are grna_group_idxs, all_nt_idxs, indiv_nt_grna_idxs.
# If in high moi, we always return the gRNA-group-to-idx map (grna_group_idxs) for the non-targeting gRNA groups; if we are running a calibration check, we additionally return the individual NT gRNA idxs (indiv_nt_grna_idxs). This latter vector is absolute (i.e., not relative to any other vector).
# If in low moi, we likewise always return the gRNA-group-to-idx map for the non-targeting gRNA groups. If the control group is the NT cells, we additionally return all_nt_idxs, which is the set of NT cell idxs. Finally, if we are running a calibration check, we return the indices of the individual NT gRNAs. If the control group is the NT cells, then these indices are relative to the NT cells. If the control group is the complement set, then these indices are absolute.
assign_grnas_to_cells <- function(grna_matrix, grna_group_data_frame, grna_assign_threshold, low_moi, control_group_complement, calibration_check) {
  if (low_moi) {
    grna_assignments <- assign_grnas_to_cells_lowmoi(grna_matrix, grna_group_data_frame, calibration_check, control_group_complement)
  } else {
    grna_assignments <- assign_grnas_to_cells_highmoi(grna_matrix, grna_assign_threshold, grna_group_data_frame, calibration_check)
  }
 return(grna_assignments)
}


assign_grnas_to_cells_highmoi <- function(grna_matrix, threshold, grna_group_data_frame, calibration_check) {
  out <- list()

  # 1. make the matrix row-accessible
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)

  # 2. obtain the grna ids and groups
  grna_ids <- rownames(grna_matrix)
  grna_groups <- as.character(unique(grna_group_data_frame$grna_group))
  grna_groups <- grna_groups[grna_groups != "non-targeting"]

  # 3. loop over the targeting grna groups, obtaining the cell assignments for each
  grna_group_idxs <- sapply(grna_groups, function(grna_group) {
    l <- grna_group_data_frame$grna_group == grna_group
    curr_grna_ids <- grna_group_data_frame$grna_id[l]
    row_idxs <- match(x = curr_grna_ids, grna_ids)
    cell_idxs <- group_and_threshold(j = grna_matrix@j, p = grna_matrix@p, x = grna_matrix@x,
                                     row_idxs = row_idxs, threshold = threshold)
    return(cell_idxs)
  })
  out$grna_group_idxs <- grna_group_idxs

  # 4. if running a calibration check, also the individual NT gRNAs
  if (calibration_check) {
    nt_grnas <- grna_group_data_frame |>
      dplyr::filter(grna_group == "non-targeting") |> dplyr::pull("grna_id")
    indiv_nt_grna_idxs <- sapply(nt_grnas, function(nt_grna) {
      row_idx <- match(x = nt_grna, grna_ids)
      cell_idxs <- group_and_threshold(j = grna_matrix@j, p = grna_matrix@p, x = grna_matrix@x,
                                       row_idxs = row_idx, threshold = threshold)
      return(cell_idxs)
    })
    out$indiv_nt_grna_idxs <- indiv_nt_grna_idxs
  }

  return(out)
}


assign_grnas_to_cells_lowmoi <- function(grna_matrix, grna_group_data_frame, calibration_check, control_group_complement) {
  out <- list()
  # 1. make grna matrix column accessible
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = FALSE)

  # 2. get the individual grna assignments
  grna_idx <- compute_colwise_max(i = grna_matrix@i, p = grna_matrix@p,
                                  x = grna_matrix@x, n_cells = ncol(grna_matrix))
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
  out$grna_group_idxs <- grna_group_idxs

  # 5. if we are running a calibration check, return the indices of the individual nt grnas
  if (calibration_check) {
    nt_grnas <- grna_group_data_frame |>
      dplyr::filter(grna_group == "non-targeting") |> dplyr::pull("grna_id")
    indiv_nt_grna_idxs <- lapply(nt_grnas, function(nt_grna) {
      which(nt_grna == indiv_grna_id_assignments)
    }) |> stats::setNames(nt_grnas)
    out$indiv_nt_grna_idxs <- indiv_nt_grna_idxs
  }

  # 6. handle NT cells as the control group
  if (!control_group_complement) { # NT cell control group
    if (!calibration_check) { # discovery analysis; just return the indices of the non-targeting cells
      all_nt_idxs <- which(grna_group_assignments == "non-targeting")
    } else { # calibration check; update indiv_nt_grna_idxs so that they are relative to all_nt_idxs
      all_nt_idxs <- stats::setNames(unlist(indiv_nt_grna_idxs), NULL)
      n_cells_per_nt <- sapply(indiv_nt_grna_idxs, length)
      stop <- cumsum(n_cells_per_nt)
      start <- c(0L, stop[-length(stop)]) + 1L
      indiv_nt_grna_idxs <- lapply(seq(1, length(nt_grnas)), function(i) {
        seq(start[i], stop[i])
      }) |> stats::setNames(nt_grnas)
      out$indiv_nt_grna_idxs <- indiv_nt_grna_idxs
    }
    out$all_nt_idxs <- all_nt_idxs
  }
  return(out)
}
