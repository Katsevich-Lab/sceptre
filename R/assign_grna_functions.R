###########################################
# TOP-LEVEL FUNCTION: ASSIGN GRNAS TO CELLS
###########################################
assign_grnas_to_cells <- function(sceptre_object, print_progress, parallel, n_processors, log_dir) {
  # extract pieces from sceptre_object
  grna_matrix <- sceptre_object@grna_matrix
  grna_assignment_method <- sceptre_object@grna_assignment_method
  cell_covariate_data_frame <- sceptre_object@covariate_data_frame
  grna_assignment_hyperparameters <- sceptre_object@grna_assignment_hyperparameters

  # determine the gRNAs to analyze
  if (sceptre_object@nf_pipeline) {
    grnas_in_use <- sceptre_object@elements_to_analyze
  } else {
    grnas_in_use <- determine_grnas_in_use(sceptre_object)
  }

  # assign grnas via the selected strategy; obtain the grna assignments and cells containing multiple grnas
  if (grna_assignment_method == "mixture") {
    sceptre_object@initial_grna_assignment_list <- assign_grnas_to_cells_mixture(
      grna_matrix = grna_matrix, cell_covariate_data_frame = cell_covariate_data_frame,
      grna_assignment_hyperparameters = grna_assignment_hyperparameters,
      print_progress = print_progress, parallel = parallel, n_processors = n_processors,
      log_dir = log_dir, grna_ids = grnas_in_use
    )
  }
  if (grna_assignment_method == "thresholding") {
    sceptre_object@initial_grna_assignment_list <- assign_grnas_to_cells_thresholding(
      grna_matrix = grna_matrix,
      grna_assign_threshold = grna_assignment_hyperparameters$threshold,
      grna_ids = grnas_in_use
    )
  }
  if (grna_assignment_method == "maximum") {
    max_result <- assign_grnas_to_cells_maximum(
      grna_matrix = grna_matrix,
      grna_lib_size = cell_covariate_data_frame$grna_n_umis,
      umi_fraction_threshold = grna_assignment_hyperparameters$umi_fraction_threshold
    )
    sceptre_object@initial_grna_assignment_list <- max_result$initial_assignment_list
    sceptre_object@cells_w_multiple_grnas <- max_result$cells_w_multiple_grnas # set cells w/ multiple gRNAs for max method
  }

  # process the initial assignment list
  if (!sceptre_object@nf_pipeline) sceptre_object <- process_initial_assignment_list(sceptre_object)
  return(sceptre_object)
}


#######################
# THRESHOLDING FUNCTION
#######################
assign_grnas_to_cells_thresholding <- function(grna_matrix, grna_assign_threshold, grna_ids) {
  # take cases on the class of grna_matrix
  if (methods::is(grna_matrix, "odm")) {
    initial_assignment_list <- sapply(grna_ids, function(grna_id) {
      ondisc:::threshold_count_matrix_cpp(file_name_in = grna_matrix@h5_file,
                                          f_row_ptr = grna_matrix@ptr,
                                          row_idx = which(grna_id == rownames(grna_matrix)),
                                          threshold = grna_assign_threshold)
    })
  } else {
    # 1. make the grna expression matrix row-accessible; ensure threshold is numeric
    grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
    grna_assign_threshold <- as.numeric(grna_assign_threshold)

    # 2. perform the assignments
    initial_assignment_list <- sapply(grna_ids, function(grna_id) {
      threshold_count_matrix(j = grna_matrix@j, p = grna_matrix@p, x = grna_matrix@x,
                             row_idx = which(grna_id == rownames(grna_matrix)),
                             threshold = grna_assign_threshold)
    })
  }

  return(initial_assignment_list)
}

#########
# MAXIMUM
#########
assign_grnas_to_cells_maximum <- function(grna_matrix, grna_lib_size, umi_fraction_threshold) {
  # take cases on the class of grna_matrix
  if (methods::is(grna_matrix, "odm")) {
    grna_ids <- unique(sceptre_object@ondisc_grna_assignment_info$max_grna)
    initial_assignment_list <- lapply(grna_ids, function(grna_id) {
      which(sceptre_object@ondisc_grna_assignment_info$max_grna == grna_id)
    }) |> stats::setNames(grna_ids)
    cells_w_multiple_grnas <- which(sceptre_object@ondisc_grna_assignment_info$max_grna_frac_umis <= umi_fraction_threshold)
  } else {
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
  }

  return(list(initial_assignment_list = initial_assignment_list,
              cells_w_multiple_grnas = cells_w_multiple_grnas))
}


######################
# AUXILLIARY FUNCTIONS
######################
preprocess_initial_assignment_list_vector_supplied <- function(sceptre_object) {
  # update the initial grna assignment list
  initial_grna_assignment_list <- sceptre_object@initial_grna_assignment_list
  grna_target_data_table <- sceptre_object@grna_target_data_frame_with_vector |> data.table::as.data.table()
  vector_ids <- unique(grna_target_data_table$vector_id)
  vector_groups <- lapply(vector_ids, function(curr_vector_id) {
    grna_target_data_table[grna_target_data_table$vector_id == curr_vector_id,]$grna_id
  }) |> stats::setNames(vector_ids)
  initial_grna_assignment_list_modified <- lapply(vector_groups, function(elem) {
    l <- initial_grna_assignment_list[elem]
    unique(unlist(l))
  })
  # update the grna target data table
  grna_target_data_table_modified <- grna_target_data_table |>
    dplyr::select(grna_id = vector_id, grna_target) |>
    dplyr::mutate(grna_group = grna_target) |>
    dplyr::distinct()
  # update the fields of the sceptre_object and return
  sceptre_object@grna_target_data_frame <- grna_target_data_table_modified
  sceptre_object@initial_grna_assignment_list <- initial_grna_assignment_list_modified
  return(sceptre_object)
}


# add three fields to sceptre_object:
# 1. grnas_per_cell, 2. cells_w_multiple_grnas, 3. grna_assignments_raw
process_initial_assignment_list <- function(sceptre_object) {
  # -1. preprocessing step: if vector_id has been supplied, combine grnas within the same vector, and treat vectors as if they are individual grnas
  if (nrow(sceptre_object@grna_target_data_frame_with_vector) >= 1L) {
    sceptre_object <- preprocess_initial_assignment_list_vector_supplied(sceptre_object)
  }
  # 0. set variables
  initial_assignment_list <- sceptre_object@initial_grna_assignment_list
  grna_target_data_frame <- sceptre_object@grna_target_data_frame
  low_moi <- sceptre_object@low_moi
  n_cells <- ncol(sceptre_object@grna_matrix)
  maximum_assignment <- sceptre_object@grna_assignment_method == "maximum"
  # 1. if not using maximum assignment, compute n grnas per cell and cells with multiple grnas
  if (!maximum_assignment) {
    sceptre_object@grnas_per_cell <- compute_n_grnas_per_cell_vector(initial_assignment_list, n_cells)
    if (low_moi) {
      sceptre_object@cells_w_multiple_grnas <- which(sceptre_object@grnas_per_cell >= 2L)
    }
  }
  # 2. pool together targeting gRNAs via the or operation
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
  # 6. initialize output
  sceptre_object@grna_assignments_raw <- grna_assignments_raw
  return(sceptre_object)
}


determine_grnas_in_use <- function(sceptre_object) {
  all_grna_targets <- unique(c(sceptre_object@positive_control_pairs$grna_target,
                               sceptre_object@discovery_pairs$grna_target, "non-targeting"))
  if (nrow(sceptre_object@grna_target_data_frame_with_vector) >= 1L) {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame_with_vector
  } else {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame
  }
  grnas_in_use <- dplyr::filter(grna_target_data_frame, grna_target %in% all_grna_targets) |>
    dplyr::pull(grna_id)
  return(grnas_in_use)
}
