run_calibration_check_v2 <- function(sceptre_object, output_amount = 1, n_calibration_pairs = NULL,
                                     calibration_group_size = NULL, print_progress = TRUE, parallel = FALSE) {
  # 0. verify that function called in correct order
  check_function_call(sceptre_object, "run_calibration_check")

  # 1. handle the default arguments
  if (is.null(calibration_group_size)) calibration_group_size <- compute_calibration_group_size(sceptre_object@grna_group_data_frame)
  if (is.null(n_calibration_pairs)) n_calibration_pairs <- sceptre_object@n_ok_discovery_pairs

  # 2. check inputs
  check_calibration_check_inputs(sceptre_object) |> invisible()

  # 3. construct the negative control pairs
  cat("Constructing negative control pairs.")
  response_grna_group_pairs <- construct_negative_control_pairs(n_calibration_pairs = n_calibration_pairs,
                                                                calibration_group_size = calibration_group_size,
                                                                grna_assignments = sceptre_object@grna_assignments,
                                                                response_matrix = sceptre_object@response_matrix,
                                                                grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                                                low_moi = sceptre_object@low_moi,
                                                                n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                                                n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                                                n_nonzero_m = sceptre_object@M_matrix,
                                                                n_nonzero_tot = sceptre_object@n_nonzero_tot_vector,
                                                                cells_in_use = sceptre_object@cells_in_use)
  cat(crayon::green(' \u2713\n'))



}














run_calibration_check <- function(sceptre_object, output_amount = 1, n_calibration_pairs = NULL,
                                  calibration_group_size = NULL, print_progress = TRUE, parallel = FALSE) {
  # 1. do argument check
  check_calibration_check_inputs(analysis_prepared = sceptre_object@analysis_prepared,
                                 grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                 control_group_complement = sceptre_object@control_group_complement) |> invisible()

  # 2. get grna assignments
  grna_assignments <- sceptre_object@grna_assignments

  # 3. construct the negative control pairs
  if (is.na(sceptre_object@calibration_group_size)) {
    calibration_group_size <- compute_calibration_group_size(sceptre_object@grna_group_data_frame)
  }  else {
    calibration_group_size <- sceptre_object@calibration_group_size
  }

  cat("Constructing negative control pairs.")
  if (nrow(sceptre_object@negative_control_pairs) == 0L) {
    response_grna_group_pairs <- construct_negative_control_pairs(n_calibration_pairs = sceptre_object@n_calibration_pairs,
                                                                  calibration_group_size = calibration_group_size,
                                                                  grna_assignments = grna_assignments,
                                                                  response_matrix = sceptre_object@response_matrix,
                                                                  grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                                                  low_moi = sceptre_object@low_moi,
                                                                  n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                                                  n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                                                  n_nonzero_m = sceptre_object@M_matrix,
                                                                  n_nonzero_tot = sceptre_object@n_nonzero_tot_vector)
  } else {
    response_grna_group_pairs <- sceptre_object@negative_control_pairs
  }
  cat(crayon::green(' \u2713\n'))

  # 4. generate the set of synthetic indicator idxs (if running permutations)
  if (sceptre_object@run_permutations) {
    cat("Generating permutation resamples.")
    n_cells <- nrow(sceptre_object@covariate_matrix)
    synthetic_idxs <- get_synthetic_permutation_idxs(grna_assignments = grna_assignments,
                                                     B = sceptre_object@B1 + sceptre_object@B2 + sceptre_object@B3,
                                                     calibration_check = TRUE,
                                                     control_group_complement = sceptre_object@control_group_complement,
                                                     calibration_group_size = calibration_group_size,
                                                     n_cells = n_cells)
    cat(crayon::green(' \u2713\n'))
  }
  gc() |> invisible()

  # 5. run the method
  if (sceptre_object@run_permutations) {
    out <- run_perm_test_in_memory(response_matrix = sceptre_object@response_matrix,
                                   grna_assignments = grna_assignments,
                                   covariate_matrix = sceptre_object@covariate_matrix,
                                   response_grna_group_pairs = response_grna_group_pairs,
                                   synthetic_idxs = synthetic_idxs,
                                   output_amount = output_amount,
                                   fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                   B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                   B3 = sceptre_object@B3, calibration_check = TRUE,
                                   control_group_complement = sceptre_object@control_group_complement,
                                   n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                   n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                   side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                   response_precomputations = sceptre_object@response_precomputations,
                                   print_progress = print_progress, parallel = parallel)
  } else {
    stop("CRT not yet implemented.")
    ret <- run_crt_in_memory_v2(response_matrix = sceptre_object@response_matrix,
                                grna_assignments = grna_assignments,
                                covariate_matrix = sceptre_object@covariate_matrix,
                                response_grna_group_pairs = response_grna_group_pairs,
                                output_amount = output_amount,
                                fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                B3 = sceptre_object@B3, calibration_check = TRUE,
                                control_group_complement = sceptre_object@control_group_complement,
                                n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                print_progress = print_progress)
  }

  # update fields of sceptre object
  sceptre_object@calibration_result <- out$ret_pass_qc
  sceptre_object@negative_control_pairs <- response_grna_group_pairs
  sceptre_object@response_precomputations <- out$response_precomputations
  sceptre_object@calibration_check_run <- TRUE
  sceptre_object@last_function_called <- "run_calibration_check"
  return(sceptre_object)
}


run_power_check <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE) {
  # 1. do argument check
  check_discovery_analysis_inputs(response_grna_group_pairs = sceptre_object@positive_control_pairs_with_info,
                                  control_group_complement = sceptre_object@control_group_complement,
                                  grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                  calibration_check_run = sceptre_object@calibration_check_run,
                                  pc_analysis = TRUE) |> invisible()

  # 2. get grna assignments
  grna_assignments <- sceptre_object@grna_assignments

  # 4. generate the set of synthetic indicator idxs (if running permutations)
  if (sceptre_object@run_permutations) {
    cat("Generating permutation resamples.")
    n_cells <- nrow(sceptre_object@covariate_matrix)
    synthetic_idxs <- get_synthetic_permutation_idxs(grna_assignments = grna_assignments,
                                                     B = sceptre_object@B1 + sceptre_object@B2 + sceptre_object@B3,
                                                     calibration_check = FALSE,
                                                     control_group_complement = sceptre_object@control_group_complement,
                                                     calibration_group_size = NULL,
                                                     n_cells = n_cells)
    cat(crayon::green(' \u2713\n'))
  }
  gc() |> invisible()

  # 5. run the method
  if (sceptre_object@run_permutations) {
    out <- run_perm_test_in_memory(response_matrix = sceptre_object@response_matrix,
                                   grna_assignments = grna_assignments,
                                   covariate_matrix = sceptre_object@covariate_matrix,
                                   response_grna_group_pairs = sceptre_object@positive_control_pairs_with_info |> dplyr::filter(pass_qc),
                                   synthetic_idxs = synthetic_idxs,
                                   output_amount = output_amount,
                                   fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                   B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                   B3 = sceptre_object@B3, calibration_check = FALSE,
                                   control_group_complement = sceptre_object@control_group_complement,
                                   n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                   n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                   side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                   response_precomputations = sceptre_object@response_precomputations,
                                   print_progress = print_progress, parallel = parallel)
  } else {
    stop("CRT not yet implemented.")
    ret <- run_crt_in_memory_v2(response_matrix = sceptre_object@response_matrix,
                                grna_assignments = grna_assignments,
                                covariate_matrix = sceptre_object@covariate_matrix,
                                response_grna_group_pairs = sceptre_object@positive_control_pairs_with_info,
                                output_amount = output_amount,
                                fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                B3 = sceptre_object@B3, calibration_check = FALSE,
                                control_group_complement = sceptre_object@control_group_complement,
                                n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                print_progress = print_progress)
  }
  # 6. create the final results data frame by combining the pairs that pass qc with those that do not
  final_result <- data.table::rbindlist(list(out$ret_pass_qc,
                                             sceptre_object@positive_control_pairs_with_info |> dplyr::filter(!pass_qc)), fill = TRUE)
  data.table::setorderv(final_result, cols = c("p_value", "response_id"), na.last = TRUE)

  # update fields of sceptre object
  sceptre_object@power_result <- out$ret
  sceptre_object@response_precomputations <- out$response_precomputations
  sceptre_object@power_check_run <- TRUE
  sceptre_object@last_function_called <- "run_power_check"
  return(sceptre_object)
}


run_discovery_analysis <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE) {
  # 1. do argument check
  check_discovery_analysis_inputs(response_grna_group_pairs = sceptre_object@discovery_pairs_with_info,
                                  control_group_complement = sceptre_object@control_group_complement,
                                  grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                  calibration_check_run = sceptre_object@calibration_check_run,
                                  pc_analysis = FALSE) |> invisible()

  # 2. get grna assignments
  grna_assignments <- sceptre_object@grna_assignments

  # 4. generate the set of synthetic indicator idxs (if running permutations)
  if (sceptre_object@run_permutations) {
    cat("Generating permutation resamples.")
    n_cells <- nrow(sceptre_object@covariate_matrix)
    synthetic_idxs <- get_synthetic_permutation_idxs(grna_assignments = grna_assignments,
                                                     B = sceptre_object@B1 + sceptre_object@B2 + sceptre_object@B3,
                                                     calibration_check = FALSE,
                                                     control_group_complement = sceptre_object@control_group_complement,
                                                     calibration_group_size = NULL,
                                                     n_cells = n_cells)
    cat(crayon::green(' \u2713\n'))
  }
  gc() |> invisible()

  # 5. run the method on the pairs passing qc
  if (sceptre_object@run_permutations) {
    out <- run_perm_test_in_memory(response_matrix = sceptre_object@response_matrix,
                                   grna_assignments = grna_assignments,
                                   covariate_matrix = sceptre_object@covariate_matrix,
                                   response_grna_group_pairs = sceptre_object@discovery_pairs_with_info |> dplyr::filter(pass_qc),
                                   synthetic_idxs = synthetic_idxs,
                                   output_amount = output_amount,
                                   fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                   B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                   B3 = sceptre_object@B3, calibration_check = FALSE,
                                   control_group_complement = sceptre_object@control_group_complement,
                                   n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                   n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                   side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                   response_precomputations = sceptre_object@response_precomputations,
                                   print_progress = print_progress, parallel = parallel)
  } else {
    stop("CRT not yet implemented.")
    ret <- run_crt_in_memory_v2(response_matrix = sceptre_object@response_matrix,
                                grna_assignments = grna_assignments,
                                covariate_matrix = sceptre_object@covariate_matrix,
                                response_grna_group_pairs = sceptre_object@discovery_pairs_with_info,
                                output_amount = output_amount,
                                fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                B3 = sceptre_object@B3, calibration_check = FALSE,
                                control_group_complement = sceptre_object@control_group_complement,
                                n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                print_progress = print_progress)
  }

  # 6. create the final results data frame by combining the pairs that pass qc with those that do not
  final_result <- data.table::rbindlist(list(out$ret_pass_qc,
                                             sceptre_object@discovery_pairs_with_info |> dplyr::filter(!pass_qc)), fill = TRUE)
  data.table::setorderv(final_result, cols = c("p_value", "response_id"), na.last = TRUE)

  # update fields of sceptre object
  sceptre_object@discovery_result <- final_result
  sceptre_object@response_precomputations <- out$response_precomputations
  sceptre_object@discovery_analysis_run <- TRUE
  sceptre_object@last_function_called <- "run_discovery_analysis"
  return(sceptre_object)
}


get_result <- function(sceptre_object, analysis_type, alpha = 0.1, multiple_testing_correction = "BH") {
  if (!(analysis_type %in% c("calibration", "power", "discovery"))) {
    stop("`analysis_type` must be one of `calibration`, `power`, or `discovery`.")
  }
  if (analysis_type == "calibration" && !sceptre_object@calibration_check_run) {
    stop("Calibration check has not yet been run.")
  }
  if (analysis_type == "power" && !sceptre_object@power_check_run) {
    stop("Power check has not yet been run.")
  }
  if (analysis_type == "discovery" && !sceptre_object@discovery_analysis_run) {
    stop("Discovery analysis not yet run.")
  }
  field_to_extract <- switch(EXPR = analysis_type,
                             calibration = "calibration_result",
                             power = "power_result",
                             discovery = "discovery_result")
  out <- slot(sceptre_object, field_to_extract)
  # if calibration check or discovery analysis, apply multiplicity correction
  if (analysis_type %in% c("calibration", "discovery")) {
    out <- out |>
      dplyr::mutate(reject = stats::p.adjust(p_value, method = multiple_testing_correction) < alpha) |>
      dplyr::select(-pass_qc)
  }
  return(out)
}
