run_calibration_check <- function(sceptre_object, output_amount = 1, n_calibration_pairs = NULL,
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
  response_grna_group_pairs <- construct_negative_control_pairs_v2(sceptre_object = sceptre_object,
                                                                   n_calibration_pairs = n_calibration_pairs,
                                                                   calibration_group_size = calibration_group_size)
  cat(crayon::green(' \u2713\n'))

  # 4. update uncached objects
  sceptre_object@calibration_group_size <- calibration_group_size
  sceptre_object@n_calibration_pairs <- n_calibration_pairs

  # 5. run the sceptre analysis (high-level function call)
  out <- run_sceptre_analysis_high_level(sceptre_object = sceptre_object,
                                         response_grna_group_pairs = response_grna_group_pairs,
                                         calibration_check = TRUE,
                                         output_amount = output_amount,
                                         print_progress = print_progress,
                                         parallel = parallel)

  # 6. update fields of sceptre object with results
  sceptre_object@last_function_called <- "run_calibration_check"
  sceptre_object@calibration_result <- out$result
  sceptre_object@negative_control_pairs <- response_grna_group_pairs
  sceptre_object@response_precomputations <- out$response_precomputations
  return(sceptre_object)
}


run_power_check <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE) {
  # 0. verify that function called in correct order
  check_function_call(sceptre_object, "run_power_check")

  # 1. extract relevant arguments
  response_grna_group_pairs <- sceptre_object@positive_control_pairs_with_info

  # 2. check inputs
  check_discovery_analysis_inputs(response_grna_group_pairs = response_grna_group_pairs,
                                  control_group_complement = sceptre_object@control_group_complement,
                                  grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                  pc_analysis = TRUE) |> invisible()

  # 3.  run the sceptre analysis (high-level function call)
  out <- run_sceptre_analysis_high_level(sceptre_object = sceptre_object,
                                         response_grna_group_pairs = response_grna_group_pairs,
                                         calibration_check = FALSE,
                                         output_amount = output_amount,
                                         print_progress = print_progress,
                                         parallel = parallel)

  # 4. update fields of sceptre object with results
  sceptre_object@last_function_called <- "run_power_check"
  sceptre_object@power_result <- out$result
  sceptre_object@response_precomputations <- out$response_precomputations
  return(sceptre_object)
}


run_discovery_analysis <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE) {
  # 0. verify that function called in correct order
  check_function_call(sceptre_object, "run_discovery_analysis")

  # 1. extract relevant arguments
  response_grna_group_pairs <- sceptre_object@discovery_pairs_with_info

  # 2. check inputs
  check_discovery_analysis_inputs(response_grna_group_pairs = response_grna_group_pairs,
                                  control_group_complement = sceptre_object@control_group_complement,
                                  grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                  pc_analysis = FALSE) |> invisible()

  # 3.  run the sceptre analysis (high-level function call)
  out <- run_sceptre_analysis_high_level(sceptre_object = sceptre_object,
                                         response_grna_group_pairs = response_grna_group_pairs,
                                         calibration_check = FALSE,
                                         output_amount = output_amount,
                                         print_progress = print_progress,
                                         parallel = parallel)

  # 4. update fields of sceptre object with results
  sceptre_object@last_function_called <- "run_discovery_analysis"
  sceptre_object@discovery_result <- out$result
  sceptre_object@response_precomputations <- out$response_precomputations
  return(sceptre_object)
}


run_sceptre_analysis_high_level <- function(sceptre_object, response_grna_group_pairs, calibration_check, output_amount, print_progress, parallel) {
  # if running permutations, generate the permutation idxs
  if (sceptre_object@run_permutations) {
    cat("Generating permutation resamples.")
    synthetic_idxs <- get_synthetic_permutation_idxs(grna_assignments = sceptre_object@grna_assignments,
                                                     B = sceptre_object@B1 + sceptre_object@B2 + sceptre_object@B3,
                                                     calibration_check = calibration_check,
                                                     control_group_complement = sceptre_object@control_group_complement,
                                                     calibration_group_size = sceptre_object@calibration_group_size,
                                                     n_cells = length(sceptre_object@cells_in_use))
    cat(crayon::green(' \u2713\n'))
  }
  gc() |> invisible()

  # initialize the args to pass
  args_to_pass <- list(response_matrix = sceptre_object@response_matrix,
                       grna_assignments = sceptre_object@grna_assignments,
                       covariate_matrix = sceptre_object@covariate_matrix,
                       response_grna_group_pairs = response_grna_group_pairs,
                       output_amount = output_amount,
                       fit_parametric_curve = sceptre_object@fit_parametric_curve,
                       B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                       B3 = sceptre_object@B3, calibration_check = calibration_check,
                       control_group_complement = sceptre_object@control_group_complement,
                       n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                       n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                       side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                       response_precomputations = sceptre_object@response_precomputations,
                       cells_in_use = sceptre_object@cells_in_use,
                       print_progress = print_progress, parallel = parallel)

  # run the method
  out <- if (sceptre_object@run_permutations) {
    args_to_pass$synthetic_idxs <- synthetic_idxs
    do.call(what = "run_perm_test_in_memory", args = args_to_pass)
  } else {
    do.call(what = "run_crt_in_memory_v2", args = args_to_pass)
  }

  # return the output
  return(out)
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
