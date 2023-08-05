check_create_sceptre_object_inputs <- function(response_matrix, grna_matrix, grna_group_data_frame, moi, extra_covariates) {
  # 1. check column names of grna_group_data_frame
  colnames_present <- all(c("grna_id", "grna_group") %in% colnames(grna_group_data_frame))
  if (!colnames_present) {
    stop("The data frame `grna_group_data_frame` must have columns `grna_id` and `grna_group`. The `grna_group` column should specify the group to which each `grna_id` belongs.")
  }

  # 2. verify that the row names are unique for both response and grna modalities
  response_ids <- rownames(response_matrix)
  grna_ids <- rownames(grna_matrix)
  if (length(response_ids) != length(unique(response_ids))) stop("The rownames of the `response_matrix` must be unique.")
  if (length(grna_ids) != length(unique(grna_ids))) stop("The rownames of the `grna_matrix` must be unique.")

  # 5. ensure that the ampersand symbol (&) is absent from the grna ids; ensure that no gRNA is named "non-targeting"
  problematic_grna_ids <- grep(pattern = "&", x = grna_ids)
  if (length(problematic_grna_ids) >= 1) {
    stop(paste0("The ampersand character (&) cannot be present in the gRNA IDs. The following gRNA IDs contain an ampersand: ", paste0(grna_ids[problematic_grna_ids], collapse = ", ")))
  }
  if (any(grna_ids == "non-targeting")) {
    stop("No individual gRNA can have the ID `non-targeting`. The string `non-targeting` is reserved for the `grna_group` column of the `grna_group_data_frame`.")
  }

  # 6. check that the ids in the grna group data frame are a subset of the ids in the grna matrix
  if (!all(grna_group_data_frame$grna_id %in% rownames(grna_matrix))) {
    stop("The column `grna_id` of the `response_grna_group_pairs` data frame must be a subset of the row names of the grna expression matrix.")
  }

  # 7. check type of input matrices
  check_matrix_class <- function(input_matrix, input_matrix_name, allowed_matrix_classes) {
    ok_class <- sapply(X = allowed_matrix_classes, function(mat_class) methods::is(input_matrix, mat_class)) |> any()
    if (!ok_class) {
      stop(paste0("`", input_matrix_name, "` must be an object of class ", paste0(allowed_matrix_classes, collapse = ", "), "."))
    }
  }
  check_matrix_class(response_matrix, "response_matrix", c("matrix", "dgTMatrix", "dgCMatrix", "dgRMatrix"))
  check_matrix_class(grna_matrix, "grna_matrix", c("matrix", "dgTMatrix", "dgCMatrix", "dgRMatrix", "lgTMatrix", "lgCMatrix", "lgRMatrix"))

  # 8. check for agreement in number of cells
  check_ncells <- ncol(response_matrix) == ncol(grna_matrix)
  if (!is.null(extra_covariates)) {
    check_ncells <- check_ncells && (ncol(response_matrix) == nrow(extra_covariates))
  }
  if (!check_ncells) {
    stop("The number of cells in the `response_matrix`, `grna_matrix`, and `extra_covariates` (if supplied) must coincide.")
  }

  # 9. check that column names of extra_covariates are not already taken
  reserved_covariate_names <- c("response_n_nonzero", "response_n_umis", "response_p_mito", "grna_n_nonzero", "grna_n_umis")
  extra_covariate_names <- colnames(extra_covariates)
  if (any(extra_covariate_names %in% reserved_covariate_names)) {
    stop("The covariate names `response_n_nonzero`, `response_n_umis`, `response_p_mito`, `grna_n_nonzero`, and `grna_n_umis` are reserved. Change the column names of the `extra_covariates` data frame.")
  }

  # 10. verify that the types of the extra covariates are acceptable
  for (extra_covariate_name in extra_covariate_names) {
    v <- extra_covariates[,extra_covariate_name]
    accept_type <- is(v, "numeric") || is(v, "character") || is(v, "factor")
    if (!accept_type) {
      stop(paste0("The column `", extra_covariate_name, "` of the `extra_covariates` data frame should be of type numeric, character, or factor."))
    }
  }

  # 11. verify that moi is specified
  if (!(moi %in% c("low", "high"))) {
    stop("`moi` should be either `low` or `high`.")
  }

  return(NULL)
}


check_prepare_analysis_inputs <- function(response_matrix, grna_matrix, covariate_data_frame,
                                          grna_group_data_frame, formula_object, response_grna_group_pairs_list,
                                          control_group, resampling_mechanism, side, low_moi) {
  # 5. if response_grna_group_pairs has been supplied, check its characteristics
  for (response_grna_group_pairs in response_grna_group_pairs_list) {
    if (nrow(response_grna_group_pairs) >= 1L) {
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


check_calibration_check_inputs <- function(analysis_prepared, grna_group_data_frame, control_group_complement) {
  if (!analysis_prepared) {
    stop("Prepare the analysis via prepare_analysis() before running a calibration check.")
  }
  n_nt_grnas <- grna_group_data_frame |>
    dplyr::filter(grna_group == "non-targeting") |> nrow()
  if (control_group_complement) {
    if (n_nt_grnas < 2) stop("Two or more non-targeting gRNAs must be present to run the calibration check. gRNAs that are non-targeting should be assigned a gRNA group label of 'non-targeting' in the `grna_group_data_frame`.")
  } else {
    if (n_nt_grnas < 1) stop("At least one non-targeting gRNA must be present to run the calibration check. gRNAs that are non-targeting should be assigned a gRNA group label of 'non-targeting' in the `grna_group_data_frame`.")
  }

  return(NULL)
}


check_discovery_analysis_inputs <- function(response_grna_group_pairs,
                                            control_group_complement,
                                            grna_group_data_frame,
                                            calibration_check_run,
                                            pc_analysis) {
  # 1. check that positive control pairs are available
  if (nrow(response_grna_group_pairs) == 0L) {
    stop(paste0(ifelse(pc_analysis, "Positive control", "Discovery"), " pairs have not been supplied. Thus, the ", ifelse(pc_analysis, "power check", "discovery analysis"), " cannot be run. You can supply ", ifelse(pc_analysis, "positive control", "discovery"), " pairs in the function prepare_analysis()."))
  }

  # 2. check that negative control gRNAs are present (if the control group is the complement set)
  if (control_group_complement) {
    nt_present <- "non-targeting" %in% grna_group_data_frame$grna_group
    if (!nt_present) {
      stop(paste0("At least one non-targeting gRNA must be present to run a ", ifelse(pc_analysis, "power check", "discovery analysis"), " when the control group is the complement set."))
    }
  }

  # 3. check that the calibration check has been run
  if (!sceptre_object@calibration_check_run) {
    warning("Either the calibration check has not yet been run, or the analysis parameters have been updated since the calibration check was last run. We recommend running the calibration check (via the calibration_check() function) before continuing with the ", ifelse(pc_analysis, "power check", "discovery analysis"), ".")
  }

  return(NULL)
}
