check_create_sceptre_object_inputs <- function(response_matrix, grna_matrix, covariate_data_frame,
                                               grna_group_data_frame, moi) {
  # 1. check column names of grna_group_data_frame
  colnames_present <- all(c("grna_id", "grna_group") %in% colnames(grna_group_data_frame))
  if (!colnames_present) {
    stop("The data frame `grna_group_data_frame` must have columns `grna_id` and `grna_group`. The `grna_group` column should specify the group to which each `grna_id` belongs.")
  }

  # 3. verify that the row names are unique for both response and grna modalities
  response_ids <- rownames(response_matrix)
  grna_ids <- rownames(grna_matrix)
  if (length(response_ids) != length(unique(response_ids))) stop("The rownames of the `response_matrix` must be unique.")
  if (length(grna_ids) != length(unique(grna_ids))) stop("The rownames of the `grna_matrix` must be unique.")

  # 4. ensure that the ampersand symbol (&) is absent from the grna ids; ensure that no gRNA is named "non-targeting"
  problematic_grna_ids <- grep(pattern = "&", x = grna_ids)
  if (length(problematic_grna_ids) >= 1) {
    stop(paste0("The ampersand character (&) cannot be present in the gRNA IDs. The following gRNA IDs contain an ampersand: ", paste0(grna_ids[problematic_grna_ids], collapse = ", ")))
  }
  if (any(grna_ids == "non-targeting")) {
    stop("No individual gRNA can have the ID `non-targeting`. The string `non-targeting` is reserved for the `grna_group` column of the `grna_group_data_frame`.")
  }

  # 5. check that the ids in the grna group data frame are a subset of the ids in the grna matrix
  if (!all(grna_group_data_frame$grna_id %in% rownames(grna_matrix))) {
    stop("The column `grna_id` of the `response_grna_group_pairs` data frame must be a subset of the row names of the grna expression matrix.")
  }

  # 7. check the covariate data frame has at least one column
  if (ncol(covariate_data_frame) == 0) {
    stop("The global cell covariate matrix must contain at least one column.")
  }

  # 9. check type of input matrices
  check_matrix_class <- function(input_matrix, input_matrix_name, allowed_matrix_classes) {
    ok_class <- sapply(X = allowed_matrix_classes, function(mat_class) methods::is(input_matrix, mat_class)) |> any()
    if (!ok_class) {
      stop(paste0("`", input_matrix_name, "` must be an object of class ", paste0(allowed_matrix_classes, collapse = ", "), "."))
    }
  }
  check_matrix_class(response_matrix, "response_matrix", c("matrix", "dgTMatrix", "dgCMatrix", "dgRMatrix"))
  check_matrix_class(grna_matrix, "grna_matrix", c("matrix", "dgTMatrix", "dgCMatrix", "dgRMatrix", "lgTMatrix", "lgCMatrix", "lgRMatrix"))

  # 10. check for agreement in number of cells
  check_ncells <- (ncol(response_matrix) == ncol(grna_matrix)) && (ncol(response_matrix) == nrow(covariate_data_frame))
  if (!check_ncells) {
    stop("The number of cells in the `response_matrix`, `grna_matrix`, and `covariate_data_frame` must coincide.")
  }

  # 12. verify that moi is specified
  if (!(moi %in% c("low", "high"))) {
    stop("`moi` should be either `low` or `high`.")
  }

  return(NULL)
}


check_prepare_analysis_inputs <- function(response_matrix, grna_matrix, covariate_data_frame,
                                          grna_group_data_frame, formula_object, response_grna_group_pairs,
                                          control_group, resampling_mechanism, side, low_moi) {
  # 5. if response_grna_group_pairs has been supplied, check its characteristics
  if (!is.null(response_grna_group_pairs)) {
    # i. verify that `grna_group` and `response_id` are columns
    all(c("grna_group", "response_id") %in% colnames(response_grna_group_pairs))
    # ii. check that the response ids in the `response_grna_group_pairs` data frame are a subset of the response ids
    if (!all(response_grna_group_pairs$response_id %in% rownames(response_matrix))) {
      stop("The column `response_id` of the `response_grna_group_pairs` data frame must be a subset of the row names of the response expression matrix.")
    }
    # iii. check that the `grna_group` column of the `response_grna_group_pairs` data frame is a subset of the `grna_group` column of the `grna_group_data_frame`
    if (!all(response_grna_group_pairs$grna_group %in% grna_group_data_frame$grna_group)) {
      stop("The column `grna_group` of the `response_grna_group_pairs` data frame must be a subset of the colummn `grna_group` of the `grna_group_data_frame`.")
    }
    # iv. ensure that "non-targeting" is not a group in the pairs to analyze data frame
    if ("non-targeting" %in% unique(response_grna_group_pairs$grna_group)) {
      stop("The `response_grna_group_pairs` data frame cannot contain the gRNA group `non-targeting`. To test non-targeting gRNAs, run a calibration check.")
    }
  }

  # 6. check that there are no offsets in the formula object
  if (grepl("offset", as.character(formula_object)[2])) stop("Offsets are not currently supported in formula objects.")

  # 8. check that the variables in the formula object are a subset of the column names of the covariate data frame
  formula_object_vars <- all.vars(formula_object)
  check_var <- formula_object_vars %in% colnames(covariate_data_frame)
  if (!all(check_var)) {
    stop(paste0("The variables in the `formula_object` must be a subset of the columns of the `covariate_data_frame`. Check the following variables: ", paste0(formula_object_vars[!check_var], collapse = ", ")))
  }

  # 13. verify that control group is either nt_cells, complement, or default
  if (!control_group %in% c("nt_cells", "complement", "default")) {
    stop("`control_group` should set to `nt_cells`, `complement`, or `default`.")
  }

  # 14. verify that resampling_mechanism is one of "permutations", "crt", or "default"
  if (!(resampling_mechanism %in% c("permutations", "crt", "default"))) {
    stop("`resampling_mechanism` should set to `permutations`, `crt`, or `default`.")
  }

  # 15. verify that the moi is consistent with the control group
  if (!low_moi && control_group == "nt_cells") {
    stop("The control group cannot be the NT cells in high MOI.")
  }

  # 16. verify that, if the control_group is NT cells, there are NT gRNAs present
  if (control_group == "nt_cells") {
    nt_present <- "non-targeting" %in% grna_group_data_frame$grna_group
    if (!nt_present) {
      stop(paste0("The string 'non-targeting' must be present in the `grna_group` column of the `grna_group_data_frame`."))
    }
  }

  # 17. verify that "side" is among "both", "left", or "right"
  if (!(side %in% c("both", "left", "right"))) {
    stop("'side' must be one of 'both', 'left', or 'right'.")
  }

  return(NULL)
}


check_calibration_check_inputs <- function(grna_group_data_frame, control_group_complement) {
  n_nt_grnas <- grna_group_data_frame |>
    dplyr::filter(grna_group == "non-targeting") |> nrow()
  if (control_group_complement) {
    if (n_nt_grnas < 2) stop("Two or more non-targeting gRNAs must be present to run the calibration check. gRNAs that are non-targeting should be assigned a gRNA group label of 'non-targeting' in the `grna_group_data_frame`.")
  } else {
    if (n_nt_grnas < 1) stop("At least one non-targeting gRNA must be present to run the calibration check. gRNAs that are non-targeting should be assigned a gRNA group label of 'non-targeting' in the `grna_group_data_frame`.")
  }

  return(NULL)
}
