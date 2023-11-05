assign_grnas_to_cells <- function(sceptre_object, print_progress, parallel, n_processors) {
  # extract pieces from sceptre_object
  grna_matrix <- sceptre_object@grna_matrix
  grna_target_data_frame <- sceptre_object@grna_target_data_frame
  low_moi <- sceptre_object@low_moi
  grna_assignment_method <- sceptre_object@grna_assignment_method
  cell_covariate_data_frame <- sceptre_object@covariate_data_frame
  grna_assignment_hyperparameters <- sceptre_object@grna_assignment_hyperparameters
  n_cells <- ncol(grna_matrix)
  maximum_assignment <- grna_assignment_method == "maximum"

  # assign grnas via the selected strategy; obtain the grna assignments and cells containing multiple grnas
  if (grna_assignment_method == "mixture") {
    initial_assignment_list <- assign_grnas_to_cells_mixture(grna_matrix = grna_matrix, cell_covariate_data_frame = cell_covariate_data_frame,
                                                             grna_assignment_hyperparameters = grna_assignment_hyperparameters,
                                                             print_progress = print_progress, parallel = parallel, n_processors = n_processors)
  }
  if (grna_assignment_method == "thresholding") {
    initial_assignment_list <- assign_grnas_to_cells_thresholding(grna_matrix = grna_matrix,
                                                                  grna_assign_threshold = grna_assignment_hyperparameters$threshold)
  }
  if (grna_assignment_method == "maximum") {
    max_result <- assign_grnas_to_cells_maximum(grna_matrix = grna_matrix,
                                                grna_lib_size = cell_covariate_data_frame$grna_n_umis,
                                                umi_fraction_threshold = grna_assignment_hyperparameters$umi_fraction_threshold)
    initial_assignment_list <- max_result$initial_assignment_list
  }

  # process the initial assignment list
  processed_assignment_out <- process_initial_assignment_list(initial_assignment_list = initial_assignment_list,
                                                              grna_target_data_frame = grna_target_data_frame,
                                                              n_cells = n_cells, low_moi = low_moi,
                                                              maximum_assignment = maximum_assignment)
  sceptre_object@grna_assignments_raw <- processed_assignment_out$grna_assignments_raw
  sceptre_object@grnas_per_cell <- processed_assignment_out$grnas_per_cell
  sceptre_object@cells_w_multiple_grnas <- if (!maximum_assignment) processed_assignment_out$cells_w_multiple_grnas else max_result$cells_w_multiple_grnas
  sceptre_object@initial_grna_assignment_list <- initial_assignment_list
  return(sceptre_object)
}


assign_grnas_to_cells_thresholding <- function(grna_matrix, grna_assign_threshold) {
  # 1. make the grna expression matrix row-accessible; ensure threshold is numeric
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
  grna_ids <- rownames(grna_matrix)
  grna_assign_threshold <- as.numeric(grna_assign_threshold)

  # 2. perform the assignments
  initial_assignment_list <- sapply(grna_ids, function(grna_id) {
    threshold_count_matrix(j = grna_matrix@j, p = grna_matrix@p, x = grna_matrix@x,
                           row_idx = which(grna_id == grna_ids), n_cells = ncol(grna_matrix),
                           threshold = grna_assign_threshold)
  })

  return(initial_assignment_list)
}


assign_grnas_to_cells_maximum <- function(grna_matrix, grna_lib_size, umi_fraction_threshold) {
 # 1. make grna matrix column accessible
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = FALSE)

  # 2. compute the column-wise max and fraction of reads coming from the assigned grna
  l <- compute_colwise_max(i = grna_matrix@i, p = grna_matrix@p, x = grna_matrix@x,
                           n_cells = ncol(grna_matrix), grna_lib_size = grna_lib_size)

  # 3. construct a list in which each individual grna is mapped to the set of cells containing it
  indiv_grna_id_assignments <- l$assignment_vect
  grna_rownames <- rownames(grna_matrix)
  initial_assignment_list <- lapply(X = seq_along(grna_rownames), FUN = function(i) {
    which(i == indiv_grna_id_assignments)
  }) |> stats::setNames(grna_rownames)

  # 4. obtain the cells containing multiple grnas
  cells_w_multiple_grnas <- which(l$frac_umis <= umi_fraction_threshold)

  return(list(initial_assignment_list = initial_assignment_list,
              cells_w_multiple_grnas = cells_w_multiple_grnas))
}


process_initial_assignment_list <- function(initial_assignment_list, grna_target_data_frame, n_cells, low_moi, maximum_assignment) {
  # 1. compute the vector of grnas per cell
  grnas_per_cell <- if (!maximum_assignment) {
    compute_n_grnas_per_cell_vector(initial_assignment_list, n_cells)
  } else {
    integer()
  }
  # 2. determine the cells that contain multiple grnas (if in low MOI and not using maximum_assignment)
  if (low_moi && !maximum_assignment) {
    cells_w_multiple_grnas <- which(grnas_per_cell >= 2L)
  } else {
    cells_w_multiple_grnas <- integer()
  }
  # 3. pool together targeting gRNAs via the or operation
  targeting_grna_group_data_table <- grna_target_data_frame |>
    dplyr::filter(grna_group != "non-targeting") |> data.table::as.data.table()
  targeting_grna_groups <- targeting_grna_group_data_table$grna_group |> unique()
  grna_group_idxs <- lapply(targeting_grna_groups, function(targeting_grna_group) {
    curr_grna_ids <- targeting_grna_group_data_table[
      targeting_grna_group_data_table$grna_group == targeting_grna_group,]$grna_id
    initial_assignment_list[curr_grna_ids] |> unlist() |> unique()
  }) |> stats::setNames(targeting_grna_groups)
  # 4. obtain the individual non-targeting grna idxs
  nontargeting_grna_ids <- grna_target_data_frame |>
    dplyr::filter(grna_group == "non-targeting") |> dplyr::pull(grna_id)
  indiv_nt_grna_idxs <- initial_assignment_list[nontargeting_grna_ids]
  # 5. construct the grna_group_idxs list
  grna_assignments_raw <- list(grna_group_idxs = grna_group_idxs,
                               indiv_nt_grna_idxs = indiv_nt_grna_idxs)
  # 6. compute the number of targeting grna groups per cell
  out <- list(grna_assignments_raw = grna_assignments_raw,
              grnas_per_cell = grnas_per_cell,
              cells_w_multiple_grnas = cells_w_multiple_grnas)
  return(out)
}
