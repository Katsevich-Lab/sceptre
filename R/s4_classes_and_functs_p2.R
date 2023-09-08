#' @export
run_calibration_check <- function(sceptre_object, output_amount = 1, n_calibration_pairs = NULL,
                                  calibration_group_size = NULL, print_progress = TRUE, parallel = FALSE) {
  if (!parallel) cat(crayon::red("Note: Set `parallel = TRUE` in the function call to improve speed.\n\n"))
  # 0. advance function (if necessary), and check function call
  sceptre_object <- advance_set_analysis_parameters_by_two(sceptre_object, parallel)
  check_function_call(sceptre_object, "run_calibration_check")

  # 1. handle the default arguments
  if (is.null(calibration_group_size)) calibration_group_size <- compute_calibration_group_size(sceptre_object@grna_group_data_frame)
  if (is.null(n_calibration_pairs)) n_calibration_pairs <- sceptre_object@n_ok_discovery_pairs

  # 2. check inputs
  check_calibration_check_inputs(sceptre_object, n_calibration_pairs) |> invisible()

  # 3. construct the negative control pairs
  cat("Constructing negative control pairs.")
  response_grna_group_pairs <- construct_negative_control_pairs_v2(sceptre_object = sceptre_object,
                                                                   n_calibration_pairs = n_calibration_pairs,
                                                                   calibration_group_size = calibration_group_size)
  cat(crayon::green(' \u2713\n'))

  # 4. update uncached objects
  sceptre_object@calibration_group_size <- as.integer(calibration_group_size)
  sceptre_object@n_calibration_pairs <- as.integer(n_calibration_pairs)

  # 5. run the sceptre analysis (high-level function call)
  out <- run_sceptre_analysis_high_level(sceptre_object = sceptre_object,
                                         response_grna_group_pairs = response_grna_group_pairs,
                                         calibration_check = TRUE,
                                         analysis_type = "calibration_check",
                                         output_amount = output_amount,
                                         print_progress = print_progress,
                                         parallel = parallel)

  # 6. update fields of sceptre object with results
  sceptre_object@last_function_called <- "run_calibration_check"
  sceptre_object@calibration_result <- out$result |>
    dplyr::mutate(reject = stats::p.adjust(p_value, method = sceptre_object@multiple_testing_method) < sceptre_object@multiple_testing_alpha)
  sceptre_object@negative_control_pairs <- response_grna_group_pairs
  sceptre_object@response_precomputations <- out$response_precomputations
  return(sceptre_object)
}

#' @export
run_power_check <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE) {
  if (!parallel) cat(crayon::red("Note: Set `parallel = TRUE` in the function call to improve speed.\n\n"))
  # 0. verify that function called in correct order
  sceptre_object <- advance_set_analysis_parameters_by_two(sceptre_object, parallel)
  check_function_call(sceptre_object, "run_power_check")

  # 1. extract relevant arguments
  response_grna_group_pairs <- sceptre_object@positive_control_pairs_with_info

  # 2. check inputs
  check_discovery_analysis_inputs(response_grna_group_pairs = response_grna_group_pairs,
                                  control_group_complement = sceptre_object@control_group_complement,
                                  grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                  pc_analysis = TRUE,
                                  calibration_result = sceptre_object@calibration_result,
                                  n_ok_pairs = sceptre_object@n_ok_positive_control_pairs) |> invisible()

  # 3.  run the sceptre analysis (high-level function call)
  out <- run_sceptre_analysis_high_level(sceptre_object = sceptre_object,
                                         response_grna_group_pairs = response_grna_group_pairs,
                                         calibration_check = FALSE,
                                         analysis_type = "power_check",
                                         output_amount = output_amount,
                                         print_progress = print_progress,
                                         parallel = parallel)

  # 4. update fields of sceptre object with results
  sceptre_object@last_function_called <- "run_power_check"
  sceptre_object@power_result <- out$result
  sceptre_object@response_precomputations <- out$response_precomputations
  return(sceptre_object)
}

#' @export
run_discovery_analysis <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE) {
  if (!parallel) cat(crayon::red("Note: Set `parallel = TRUE` in the function call to improve speed.\n\n"))
  # 0. verify that function called in correct order
  sceptre_object <- advance_set_analysis_parameters_by_two(sceptre_object, parallel)
  check_function_call(sceptre_object, "run_discovery_analysis")

  # 1. extract relevant arguments
  response_grna_group_pairs <- sceptre_object@discovery_pairs_with_info

  # 2. check inputs
  check_discovery_analysis_inputs(response_grna_group_pairs = response_grna_group_pairs,
                                  control_group_complement = sceptre_object@control_group_complement,
                                  grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                  pc_analysis = FALSE,
                                  calibration_result = sceptre_object@calibration_result,
                                  n_ok_pairs = sceptre_object@n_ok_discovery_pairs) |> invisible()

  # 3.  run the sceptre analysis (high-level function call)
  out <- run_sceptre_analysis_high_level(sceptre_object = sceptre_object,
                                         response_grna_group_pairs = response_grna_group_pairs,
                                         calibration_check = FALSE,
                                         analysis_type = "discovery_analysis",
                                         output_amount = output_amount,
                                         print_progress = print_progress,
                                         parallel = parallel)

  # 4. update fields of sceptre object with results
  sceptre_object@last_function_called <- "run_discovery_analysis"
  sceptre_object@discovery_result <- out$result |>
    dplyr::mutate(significant = stats::p.adjust(p_value, method = sceptre_object@multiple_testing_method) < sceptre_object@multiple_testing_alpha)
  sceptre_object@response_precomputations <- out$response_precomputations
  return(sceptre_object)
}


run_sceptre_analysis_high_level <- function(sceptre_object, response_grna_group_pairs, calibration_check, analysis_type, output_amount, print_progress, parallel) {
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
                       response_grna_group_pairs = response_grna_group_pairs |> dplyr::filter(pass_qc),
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
                       print_progress = print_progress, parallel = parallel, analysis_type = analysis_type)

  # run the method
  out <- if (sceptre_object@run_permutations) {
    args_to_pass$synthetic_idxs <- synthetic_idxs
    do.call(what = "run_perm_test_in_memory", args = args_to_pass)
  } else {
    do.call(what = "run_crt_in_memory_v2", args = args_to_pass)
  }

  # wrangle the result and return
  final_result <- data.table::rbindlist(list(out$result, response_grna_group_pairs |> dplyr::filter(!pass_qc)), fill = TRUE)
  data.table::setorderv(final_result, cols = c("p_value", "response_id"), na.last = TRUE)
  out$result <- final_result
  return(out)
}


#' @export
get_result <- function(sceptre_object, analysis) {
  if (!(analysis %in% c("run_calibration_check", "run_power_check", "run_discovery_analysis"))) {
    stop("`analysis` must be one of `run_calibration_check`, `run_power_check`, or `run_discovery_analysis`.")
  }
  funts_run <- get_funct_run_vect(sceptre_object)
  if (!funts_run[[analysis]]) stop(paste0(analysis, " has not yet been run."))
  field_to_extract <- switch(EXPR = analysis,
                             run_calibration_check = "calibration_result",
                             run_power_check = "power_result",
                             run_discovery_analysis = "discovery_result")
  out <- slot(sceptre_object, field_to_extract)
  return(out)
}


advance_set_analysis_parameters_by_two <- function(sceptre_object, parallel) {
  if (sceptre_object@last_function_called == "set_analysis_parameters") {
    cat(crayon::red("Note: Automatically running `assign_grnas()` and `run_qc()` with default options.\n\n"))
    sceptre_object <- sceptre_object |> assign_grnas(parallel = parallel) |> run_qc()
  }
  return(sceptre_object)
}
