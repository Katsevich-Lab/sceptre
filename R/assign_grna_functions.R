###########################################
# TOP-LEVEL FUNCTION: ASSIGN GRNAS TO CELLS
###########################################
assign_grnas_to_cells <- function(sceptre_object, print_progress,
                                  parallel, n_processors, log_dir) {
  # extract pieces from sceptre_object
  grna_matrix <- get_grna_matrix(sceptre_object)
  grna_assignment_method <- sceptre_object@grna_assignment_method
  cell_covariate_data_frame <- sceptre_object@covariate_data_frame
  grna_assignment_hyperparameters <- sceptre_object@grna_assignment_hyperparameters

  # determine the gRNAs to analyze
  if (sceptre_object@nf_pipeline) {
    grnas_in_use <- sceptre_object@elements_to_analyze
  } else {
    grnas_in_use <- determine_grnas_in_use(sceptre_object)
  }

  # assign grnas via the selected strategy; obtain the grna
  # assignments and cells containing multiple grnas
  if (grna_assignment_method == "mixture") {
    sceptre_object@initial_grna_assignment_list <- assign_grnas_to_cells_mixture(
      grna_matrix = grna_matrix,
      cell_covariate_data_frame = cell_covariate_data_frame,
      grna_assignment_hyperparameters = grna_assignment_hyperparameters,
      print_progress = print_progress,
      parallel = parallel,
      n_processors = n_processors,
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
      sceptre_object = sceptre_object,
      umi_fraction_threshold = grna_assignment_hyperparameters$umi_fraction_threshold,
      grna_n_umis = cell_covariate_data_frame$grna_n_umis,
      min_grna_n_umis_threshold = grna_assignment_hyperparameters$min_grna_n_umis_threshold,
      grna_ids = grnas_in_use
    )
    sceptre_object@initial_grna_assignment_list <- max_result$initial_assignment_list
    # set cells w/ multiple gRNAs for max method
    sceptre_object@cells_w_zero_or_twoplus_grnas <- max_result$cells_w_zero_or_twoplus_grnas
  }
  # check that at least some gRNAs were assigned
  if (all(vapply(sceptre_object@initial_grna_assignment_list,
                 length, FUN.VALUE = integer(1)) == 0L)) {
    warning("No gRNA was assigned to any cell. Consider setting `method` to 'thresholding' and `threshold` to a small, positive number.")
  }

  # process the initial assignment list
  if (!sceptre_object@nf_pipeline) {
    sceptre_object <- process_initial_assignment_list(sceptre_object)
  }
  return(sceptre_object)
}


#######################
# THRESHOLDING FUNCTION
#######################
assign_grnas_to_cells_thresholding <- function(grna_matrix, grna_assign_threshold, grna_ids) {
  # take cases on the class of grna_matrix
  if (methods::is(grna_matrix, "odm")) {
    initial_assignment_list <- lapply(grna_ids, function(grna_id) {
      ondisc:::threshold_count_matrix_cpp(
        file_name_in = grna_matrix@h5_file,
        f_row_ptr = grna_matrix@ptr,
        row_idx = which(grna_id == rownames(grna_matrix)),
        threshold = grna_assign_threshold
      )
    }) |> stats::setNames(grna_ids)
  } else {
    # 1. make the grna expression matrix row-accessible; ensure threshold is numeric
    grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
    grna_assign_threshold <- as.numeric(grna_assign_threshold)

    # 2. perform the assignments
    initial_assignment_list <- lapply(grna_ids, function(grna_id) {
      threshold_count_matrix(
        j = grna_matrix@j, p = grna_matrix@p, x = grna_matrix@x,
        row_idx = which(grna_id == rownames(grna_matrix)),
        threshold = grna_assign_threshold
      )
    }) |> stats::setNames(grna_ids)
  }

  return(initial_assignment_list)
}

#########
# MAXIMUM
#########
assign_grnas_to_cells_maximum <- function(sceptre_object, umi_fraction_threshold,
                                          grna_n_umis, min_grna_n_umis_threshold, grna_ids) {
  initial_assignment_list <- lapply(grna_ids, function(grna_id) {
    which(sceptre_object@import_grna_assignment_info$max_grna == grna_id)
  }) |> stats::setNames(grna_ids)
  cells_w_multiple_grnas <- which(sceptre_object@import_grna_assignment_info$max_grna_frac_umis <= umi_fraction_threshold)
  cells_w_zero_grnas <- which(grna_n_umis < min_grna_n_umis_threshold)
  cells_w_zero_or_twoplus_grnas <- union(cells_w_multiple_grnas, cells_w_zero_grnas)
  return(list(
    initial_assignment_list = initial_assignment_list,
    cells_w_zero_or_twoplus_grnas = cells_w_zero_or_twoplus_grnas
  ))
}


######################
# AUXILLIARY FUNCTIONS
######################
# add three fields to sceptre_object:
# 1. grnas_per_cell, 2. cells_w_zero_or_twoplus_grnas, 3. grna_assignments_raw
process_initial_assignment_list <- function(sceptre_object) {
  # 0. set variables
  initial_assignment_list <- sceptre_object@initial_grna_assignment_list
  grna_target_data_frame <- sceptre_object@grna_target_data_frame
  low_moi <- sceptre_object@low_moi
  n_cells <- ncol(get_grna_matrix(sceptre_object))
  maximum_assignment <- sceptre_object@grna_assignment_method == "maximum"
  # 1. if not using maximum assignment, compute n grnas per cell and cells with multiple grnas
  if (!maximum_assignment) {
    sceptre_object@grnas_per_cell <- compute_n_grnas_per_cell_vector(initial_assignment_list, n_cells)
    if (low_moi) {
      cells_w_multiple_grnas <- which(sceptre_object@grnas_per_cell >= 2L)
      cells_w_zero_grnas <- seq(1, n_cells)[-sort(unique(unlist(initial_assignment_list)))]
      sceptre_object@cells_w_zero_or_twoplus_grnas <- union(cells_w_multiple_grnas, cells_w_zero_grnas)
    }
  }
  # 2. pool together targeting gRNAs via the or operation
  targeting_grna_group_data_table <- grna_target_data_frame |>
    dplyr::filter(grna_group != "non-targeting") |>
    data.table::as.data.table()
  targeting_grna_groups <- targeting_grna_group_data_table$grna_group |> unique()
  grna_group_idxs <- lapply(targeting_grna_groups, function(targeting_grna_group) {
    curr_grna_ids <- targeting_grna_group_data_table[
      targeting_grna_group_data_table$grna_group == targeting_grna_group,
    ]$grna_id
    initial_assignment_list[curr_grna_ids] |>
      unlist() |>
      unique()
  }) |> stats::setNames(targeting_grna_groups)
  # 3. For control_group = nt_cells, keep only cells for each group that do not 
  # have gRNAs in any other targeting gRNA groups
  all_targeting_cells <- unlist(grna_group_idxs, use.names = FALSE)
  if(!sceptre_object@control_group_complement){
    index_counts <- table(all_targeting_cells)
    grna_group_idxs <- lapply(grna_group_idxs, function(vec) {
      vec[index_counts[as.character(vec)] == 1]
    }) 
  }
  # 4. obtain the individual non-targeting grna idxs and all non-targeting idxs
  nontargeting_grna_ids <- grna_target_data_frame |>
    dplyr::filter(grna_group == "non-targeting") |>
    dplyr::pull(grna_id)
  indiv_nt_grna_idxs <- initial_assignment_list[nontargeting_grna_ids]
  # For control_group = nt_cells, keep only cells for each group that do not 
  # have any targeting gRNAs
  if(!sceptre_object@control_group_complement){
    indiv_nt_grna_idxs <- lapply(indiv_nt_grna_idxs,
                                 setdiff,
                                 all_targeting_cells)
  }
  # 5. construct the grna_group_idxs list
  grna_assignments_raw <- list(
    grna_group_idxs = grna_group_idxs,
    indiv_nt_grna_idxs = indiv_nt_grna_idxs
  )
  # 6. initialize output
  sceptre_object@grna_assignments_raw <- grna_assignments_raw
  # 7. save mean cells per grna
  sceptre_object@mean_cells_per_grna <- vapply(sceptre_object@initial_grna_assignment_list,
    length,
    FUN.VALUE = integer(1)
  ) |> mean()
  return(sceptre_object)
}

determine_grnas_in_use <- function(sceptre_object, restricted_grnas = FALSE) {
  grna_target_data_frame <- sceptre_object@grna_target_data_frame
  if (restricted_grnas) {
    if (!sceptre_object@nuclear) {
      all_grna_targets <- unique(c(
        sceptre_object@positive_control_pairs$grna_target,
        sceptre_object@discovery_pairs$grna_target, "non-targeting"
      ))
    } else {
      set.seed(4)
      all_grna_targets <- c(grna_target_data_frame |>
        dplyr::filter(grna_target != "non-targeting") |>
        dplyr::sample_n(min(dplyr::n(), 30)) |> dplyr::pull(grna_target), "non-targeting")
    }
    grnas_in_use <- dplyr::filter(grna_target_data_frame, grna_target %in% all_grna_targets) |>
      dplyr::pull(grna_id)
  } else {
    grnas_in_use <- grna_target_data_frame$grna_id
  }
  return(grnas_in_use)
}
