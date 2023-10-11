#' Run calibration check
#'
#' @param sceptre_object TBD
#' @param output_amount TBD
#' @param n_calibration_pairs TBD
#' @param calibration_group_size TBD
#' @param print_progress TBD
#' @param parallel TBD
#'
#' @export
run_calibration_check <- function(sceptre_object, output_amount = 1, n_calibration_pairs = NULL,
                                  calibration_group_size = NULL, print_progress = TRUE, parallel = FALSE) {
  # 0. advance function (if necessary), and check function call
  sceptre_object <- skip_assign_grnas_and_run_qc(sceptre_object, parallel)
  sceptre_object <- perform_status_check_and_update(sceptre_object, "run_calibration_check")
  if (!parallel) cat(crayon::red("Note: Set `parallel = TRUE` in the function call to improve speed.\n\n"))

  # 1. handle the default arguments
  if (is.null(calibration_group_size)) calibration_group_size <- compute_calibration_group_size(sceptre_object@grna_target_data_frame)
  if (is.null(n_calibration_pairs)) n_calibration_pairs <- sceptre_object@n_ok_discovery_pairs

  # 2. check inputs
  check_calibration_check_inputs(sceptre_object, n_calibration_pairs) |> invisible()

  # 3. construct the negative control pairs
  response_grna_group_pairs <- construct_negative_control_pairs_v2(sceptre_object = sceptre_object,
                                                                   n_calibration_pairs = n_calibration_pairs,
                                                                   calibration_group_size = calibration_group_size)

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
  sceptre_object@calibration_result <- out$result |>
    dplyr::mutate(significant = stats::p.adjust(p_value, method = sceptre_object@multiple_testing_method) <
                    sceptre_object@multiple_testing_alpha) |>
    apply_grouping_to_result(sceptre_object)
  sceptre_object@negative_control_pairs <- response_grna_group_pairs
  sceptre_object@response_precomputations <- out$response_precomputations
  return(sceptre_object)
}


#' Run power check
#'
#' @param sceptre_object TBD
#' @param output_amount TBD
#' @param print_progress TBD
#' @param parallel TBD
#'
#' @export
run_power_check <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE) {
  # 0. verify that function called in correct order
  sceptre_object <- skip_assign_grnas_and_run_qc(sceptre_object, parallel)
  sceptre_object <- perform_status_check_and_update(sceptre_object, "run_power_check")
  if (!parallel) cat(crayon::red("Note: Set `parallel = TRUE` in the function call to improve speed.\n\n"))

  # 1. extract relevant arguments
  response_grna_group_pairs <- sceptre_object@positive_control_pairs_with_info

  # 2. check inputs
  check_discovery_analysis_inputs(response_grna_group_pairs = response_grna_group_pairs,
                                  control_group_complement = sceptre_object@control_group_complement,
                                  grna_target_data_frame = sceptre_object@grna_target_data_frame,
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
  sceptre_object@power_result <- out$result |> apply_grouping_to_result(sceptre_object)
  sceptre_object@response_precomputations <- out$response_precomputations
  return(sceptre_object)
}


#' Run discovery analysis
#'
#' @param sceptre_object TBD
#' @param output_amount TBD
#' @param print_progress TBD
#' @param parallel TBD
#'
#' @export
run_discovery_analysis <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE) {
  # 0. verify that function called in correct order
  sceptre_object <- skip_assign_grnas_and_run_qc(sceptre_object, parallel)
  sceptre_object <- perform_status_check_and_update(sceptre_object, "run_discovery_analysis")
  if (!parallel) cat(crayon::red("Note: Set `parallel = TRUE` in the function call to improve speed.\n\n"))

  # 1. extract relevant arguments
  response_grna_group_pairs <- sceptre_object@discovery_pairs_with_info

  # 2. check inputs
  check_discovery_analysis_inputs(response_grna_group_pairs = response_grna_group_pairs,
                                  control_group_complement = sceptre_object@control_group_complement,
                                  grna_target_data_frame = sceptre_object@grna_target_data_frame,
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
  sceptre_object@discovery_result <- out$result |>
    dplyr::mutate(significant = stats::p.adjust(p_value, method = sceptre_object@multiple_testing_method) <
                    sceptre_object@multiple_testing_alpha) |>
    apply_grouping_to_result(sceptre_object)
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


apply_grouping_to_result <- function(result, sceptre_object) {
  grna_grouping_strategy <- sceptre_object@grna_grouping_strategy
  if (grna_grouping_strategy == "union") {
    result <- result |> dplyr::rename("grna_target" = "grna_group")
  }
  if (grna_grouping_strategy == "singleton") {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame |>
      dplyr::select(grna_id, grna_target)
    result <- result |> dplyr::rename("grna_id" = "grna_group") |>
      dplyr::left_join(grna_target_data_frame, by = "grna_id") |>
      dplyr::relocate(response_id, grna_id, grna_target)
  }
  return (result)
}


#' Get result
#'
#' @param sceptre_object TBD
#' @param analysis TBD
#'
#' @export
get_result <- function(sceptre_object, analysis) {
  if (!(analysis %in% c("run_calibration_check", "run_power_check", "run_discovery_analysis"))) {
    stop("`analysis` must be one of `run_calibration_check`, `run_power_check`, or `run_discovery_analysis`.")
  }
  if (!sceptre_object@functs_called[[analysis]]) stop(paste0(analysis, " has not yet been run."))
  field_to_extract <- switch(EXPR = analysis,
                             run_calibration_check = "calibration_result",
                             run_power_check = "power_result",
                             run_discovery_analysis = "discovery_result")
  out <- methods::slot(sceptre_object, field_to_extract)
  return(out)
}


#' Write outputs to directory
#'
#' @param sceptre_object TBD
#' @param directory TBD
#'
#' @export
write_outputs_to_directory <- function(sceptre_object, directory) {
  # 0. create directory
  if (!dir.exists(directory)) dir.create(path = directory, recursive = TRUE)
  fs <- list.files(path = directory, full.names = TRUE)
  for (f in fs) file.remove(f)

  # 1. create analysis_summary.txt file
  summary_file_fp <- paste0(directory, "/analysis_summary.txt")
  sink(file = summary_file_fp, append = FALSE)
  print(sceptre_object)
  sink(NULL)

  # 2. determine the functions that have been called
  plotting_params <- list(device = "png", scale = 1.1, width = 5, height = 4, dpi = 330)
  functs_run <- sceptre_object@functs_called
  functs_run_plots <- functs_run[names(functs_run) %in% c("assign_grnas", "run_qc", "run_calibration_check", "run_power_check", "run_discovery_analysis")]
  functs_run_plots <- c(functs_run_plots, "grna_count_distributions" = TRUE)
  for (funct in names(functs_run_plots)) {
    if (functs_run_plots[[funct]]) {
      plotting_params$plot <- do.call(what = paste0("plot_", funct), args = list(sceptre_object))
      plotting_params$filename <- paste0(directory, "/plot_", funct, ".png")
      do.call(what = ggplot2::ggsave, args = plotting_params)
    }
  }

  # 3. save results
  for (analysis in c("run_calibration_check", "run_power_check", "run_discovery_analysis")) {
    if (functs_run_plots[[analysis]]) {
      res <- get_result(sceptre_object = sceptre_object, analysis = analysis)
      saveRDS(object = res, file = paste0(directory, "/results_", analysis, ".rds"))
    }
  }

  return(NULL)
}


skip_assign_grnas_and_run_qc <- function(sceptre_object, parallel) {
  functs_run <- sceptre_object@functs_called
  if (functs_run[["import_data"]] && functs_run[["set_analysis_parameters"]] &&
      !functs_run[["assign_grnas"]] && !functs_run[["run_qc"]]) { # advance by two
    cat(crayon::red("Note: Automatically running `assign_grnas()` and `run_qc()` with default arguments.\n\n"))
    sceptre_object <- assign_grnas(sceptre_object, parallel = parallel) |> run_qc()
  }
  return(sceptre_object)
}
