#' Run calibration check
#'
#' `run_calibration_check()` runs the calibration check. See \href{https://timothy-barry.github.io/sceptre-book/run-calibration-check.html}{Chapter 5 of the manual} for detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param output_amount (optional; default `1`) an integer taking values 1, 2, or 3 specifying the amount of information to return
#' @param n_calibration_pairs (optional) the number of negative control pairs to construct and test for association
#' @param calibration_group_size (optional) the number of negative control gRNAs to put into each negative control target
#' @param print_progress (optional; default `TRUE`) a logical indicating whether to print progress updates
#' @param parallel (optional; default `FALSE`) a logical indicating whether to run the function in parallel
#' @param n_processors (optional; default "auto") an integer specifying the number of processors to use if `parallel` is set to `TRUE`. The default, "auto," automatically detects the number of processors available on the machine.
#' @param log_dir (optional; default `tempdir()`) a string indicating the directory in which to write the log files (ignored if `parallel = FALSE`)
#' @return an updated `sceptre_object` in which the calibration check has been carried out
#'
#' @export
#' @examples
#' library(sceptredata)
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#'
#' # set analysis parameters, assign grnas, run qc, run calibration check
#' sceptre_object <- sceptre_object |>
#'   set_analysis_parameters(
#'     side = "left",
#'     resampling_mechanism = "permutations"
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc() |>
#'   run_calibration_check(
#'     n_calibration_pairs = 500,
#'     calibration_group_size = 2,
#'     parallel = TRUE,
#'     n_processors = 2
#'   )
run_calibration_check <- function(sceptre_object, output_amount = 1, n_calibration_pairs = NULL,
                                  calibration_group_size = NULL, print_progress = TRUE, parallel = FALSE,
                                  n_processors = "auto", log_dir = tempdir()) {
  sceptre_object <- sceptre_object |>
    run_calibration_check_pt_1(
      n_calibration_pairs = n_calibration_pairs,
      calibration_group_size = calibration_group_size,
      parallel = parallel,
      n_processors = n_processors
    ) |>
    run_calibration_check_pt_2(
      output_amount = output_amount,
      print_progress = print_progress,
      parallel = parallel,
      n_processors = n_processors,
      log_dir = log_dir
    )
  return(sceptre_object)
}


run_calibration_check_pt_1 <- function(sceptre_object, n_calibration_pairs = NULL, calibration_group_size = NULL,
                                       parallel = FALSE, n_processors = "auto") {
  # 0. advance function (if necessary), and check function call
  sceptre_object <- skip_assign_grnas_and_run_qc(sceptre_object, parallel, n_processors)
  sceptre_object <- perform_status_check_and_update(sceptre_object, "run_calibration_check")
  if (!parallel) cat(crayon::red("Note: Set `parallel = TRUE` in the function call to improve speed.\n\n"))

  # 1. handle the default arguments
  if (sceptre_object@n_discovery_pairs == 0L && (is.null(n_calibration_pairs) || is.null(calibration_group_size))) {
    stop("`calibration_group_size` and `n_calibration_pairs` must be supplied when discovery pairs are not specified.")
  }
  if (is.null(calibration_group_size)) calibration_group_size <- compute_calibration_group_size(sceptre_object@grna_target_data_frame)
  if (is.null(n_calibration_pairs)) {
    n_calibration_pairs <- sceptre_object@n_ok_discovery_pairs
  }
  if (!is.null(n_calibration_pairs) && sceptre_object@resampling_approximation == "no_approximation") {
    mult_fact <- if (sceptre_object@side_code == 0L) 10 else 5
    sceptre_object@B3 <- ceiling(mult_fact * n_calibration_pairs / sceptre_object@multiple_testing_alpha) |> as.integer()
  }

  # 2. check inputs
  check_calibration_check_inputs(sceptre_object, n_calibration_pairs, n_processors) |> invisible()

  # 3. construct the negative control pairs
  negative_control_pairs <- construct_negative_control_pairs_v2(
    sceptre_object = sceptre_object,
    n_calibration_pairs = n_calibration_pairs,
    calibration_group_size = calibration_group_size
  )

  # 4. update uncached objects
  sceptre_object@calibration_group_size <- as.integer(calibration_group_size)
  sceptre_object@n_calibration_pairs <- as.integer(n_calibration_pairs)
  sceptre_object@negative_control_pairs <- negative_control_pairs
  return(sceptre_object)
}


run_calibration_check_pt_2 <- function(sceptre_object, output_amount = 1, print_progress = TRUE,
                                       parallel = FALSE, n_processors = "auto", log_dir = tempdir()) {
  # 5. run the sceptre analysis (high-level function call)
  response_grna_group_pairs <- sceptre_object@negative_control_pairs
  out <- run_sceptre_analysis_high_level(
    sceptre_object = sceptre_object,
    response_grna_group_pairs = response_grna_group_pairs,
    calibration_check = TRUE,
    analysis_type = "calibration_check",
    output_amount = output_amount,
    print_progress = print_progress,
    parallel = parallel,
    n_processors = n_processors,
    log_dir = log_dir
  )

  # 6. update fields of sceptre object with results
  sceptre_object@calibration_result <- if (!sceptre_object@nf_pipeline) {
    process_calibration_result(out$result, sceptre_object)
  } else {
    out$result
  }
  sceptre_object@response_precomputations <- out$response_precomputations
  return(sceptre_object)
}


process_calibration_result <- function(result, sceptre_object) {
  result |>
    apply_grouping_to_result(sceptre_object, TRUE) |>
    dplyr::mutate(significant = stats::p.adjust(p_value, method = sceptre_object@multiple_testing_method) <
      sceptre_object@multiple_testing_alpha)
}


#' Run power check
#'
#' `run_power_check()` runs the power check. See \href{https://timothy-barry.github.io/sceptre-book/run-power-check-and-discovery-analysis.html}{Chapter 6 of the manual} for detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param output_amount (optional; default `1`) an integer taking values 1, 2, or 3 specifying the amount of information to return
#' @param print_progress (optional; default `TRUE`) a logical indicating whether to print progress updates
#' @param parallel (optional; default `FALSE`) a logical indicating whether to run the function in parallel
#' @param n_processors (optional; default "auto") an integer specifying the number of processors to use if `parallel` is set to `TRUE`. The default, "auto," automatically detects the number of processors available on the machine.
#' @param log_dir (optional; default `tempdir()`) a string indicating the directory in which to write the log files (ignored if `parallel = FALSE`)
#' @return an updated `sceptre_object` in which the power check has been carried out
#'
#' @export
#' @examples
#' library(sceptredata)
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#'
#' # set analysis parameters, assign grnas, run qc
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' sceptre_object <- sceptre_object |>
#'   set_analysis_parameters(
#'     side = "left",
#'     resampling_mechanism = "permutations",
#'     positive_control_pairs = positive_control_pairs
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc() |>
#'   run_power_check()
run_power_check <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE,
                            n_processors = "auto", log_dir = tempdir()) {
  # 0. verify that function called in correct order
  sceptre_object <- skip_assign_grnas_and_run_qc(sceptre_object, parallel, n_processors)
  sceptre_object <- perform_status_check_and_update(sceptre_object, "run_power_check")
  if (!parallel) cat(crayon::red("Note: Set `parallel = TRUE` in the function call to improve speed.\n\n"))

  # 1. extract relevant arguments
  response_grna_group_pairs <- sceptre_object@positive_control_pairs_with_info

  # 2. check inputs
  check_discovery_analysis_inputs(
    response_grna_group_pairs = response_grna_group_pairs,
    control_group_complement = sceptre_object@control_group_complement,
    grna_target_data_frame = sceptre_object@grna_target_data_frame,
    pc_analysis = TRUE,
    calibration_result = sceptre_object@calibration_result,
    n_ok_pairs = sceptre_object@n_ok_positive_control_pairs,
    n_processors = n_processors
  ) |> invisible()

  # 3.  run the sceptre analysis (high-level function call)
  out <- run_sceptre_analysis_high_level(
    sceptre_object = sceptre_object,
    response_grna_group_pairs = response_grna_group_pairs,
    calibration_check = FALSE,
    analysis_type = "power_check",
    output_amount = output_amount,
    print_progress = print_progress,
    parallel = parallel,
    n_processors = n_processors,
    log_dir = log_dir
  )

  # 4. update fields of sceptre object with results
  sceptre_object@power_result <- if (!sceptre_object@nf_pipeline) {
    out$result |> apply_grouping_to_result(sceptre_object)
  } else {
    out$result
  }
  sceptre_object@response_precomputations <- out$response_precomputations
  return(sceptre_object)
}


#' Run discovery analysis
#'
#' `run_discovery_analysis()` runs the discovery analysis. See \href{https://timothy-barry.github.io/sceptre-book/run-power-check-and-discovery-analysis.html}{Chapter 6 of the manual} for detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param output_amount (optional; default `1`) an integer taking values 1, 2, or 3 specifying the amount of information to return
#' @param print_progress (optional; default `TRUE`) a logical indicating whether to print progress updates
#' @param parallel (optional; default `FALSE`) a logical indicating whether to run the function in parallel
#' @param n_processors (optional; default "auto") an integer specifying the number of processors to use if `parallel` is set to `TRUE`. The default, "auto," automatically detects the number of processors available on the machine.
#' @param log_dir (optional; default `tempdir()`) a string indicating the directory in which to write the log files (ignored if `parallel = FALSE`)
#' @return an updated `sceptre_object` in which the discovery analysis has been carried out
#'
#' @export
#' @examples
#' library(sceptredata)
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' # set analysis parameters, assign grnas, run qc
#' discovery_pairs <- construct_cis_pairs(sceptre_object)
#' sceptre_object <- sceptre_object |>
#'   set_analysis_parameters(
#'     side = "left",
#'     resampling_mechanism = "permutations",
#'     discovery_pairs = discovery_pairs
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc() |>
#'   run_discovery_analysis(
#'     parallel = TRUE,
#'     n_processors = 2
#'   )
run_discovery_analysis <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE,
                                   n_processors = "auto", log_dir = tempdir()) {
  # 0. verify that function called in correct order
  sceptre_object <- skip_assign_grnas_and_run_qc(sceptre_object, parallel, n_processors)
  sceptre_object <- perform_status_check_and_update(sceptre_object, "run_discovery_analysis")
  if (!parallel) cat(crayon::red("Note: Set `parallel = TRUE` in the function call to improve speed.\n\n"))

  # 1. extract relevant arguments
  response_grna_group_pairs <- sceptre_object@discovery_pairs_with_info

  # 2. check inputs
  check_discovery_analysis_inputs(
    response_grna_group_pairs = response_grna_group_pairs,
    control_group_complement = sceptre_object@control_group_complement,
    grna_target_data_frame = sceptre_object@grna_target_data_frame,
    pc_analysis = FALSE,
    calibration_result = sceptre_object@calibration_result,
    n_ok_pairs = sceptre_object@n_ok_discovery_pairs,
    n_processors = n_processors
  ) |> invisible()

  # 3.  run the sceptre analysis (high-level function call)
  out <- run_sceptre_analysis_high_level(
    sceptre_object = sceptre_object,
    response_grna_group_pairs = response_grna_group_pairs,
    calibration_check = FALSE,
    analysis_type = "discovery_analysis",
    output_amount = output_amount,
    print_progress = print_progress,
    parallel = parallel,
    n_processors = n_processors,
    log_dir = log_dir
  )

  # 4. update fields of sceptre object with results
  sceptre_object@discovery_result <- if (!sceptre_object@nf_pipeline) {
    process_discovery_result(out$result, sceptre_object)
  } else {
    out$result
  }
  sceptre_object@response_precomputations <- out$response_precomputations
  return(sceptre_object)
}


process_discovery_result <- function(result, sceptre_object) {
  result |>
    apply_grouping_to_result(sceptre_object) |>
    dplyr::mutate(significant = stats::p.adjust(p_value, method = sceptre_object@multiple_testing_method) <
      sceptre_object@multiple_testing_alpha)
}


run_sceptre_analysis_high_level <- function(sceptre_object, response_grna_group_pairs, calibration_check, analysis_type,
                                            output_amount, print_progress, parallel, n_processors, log_dir) {
  # if running permutations, generate the permutation idxs
  if (sceptre_object@resampling_mechanism == "permutations") {
    cat("Generating permutation resamples.")
    synthetic_idxs <- get_synthetic_permutation_idxs(
      grna_assignments = sceptre_object@grna_assignments,
      B = sceptre_object@B1 + sceptre_object@B2 + sceptre_object@B3,
      calibration_check = calibration_check,
      control_group_complement = sceptre_object@control_group_complement,
      calibration_group_size = sceptre_object@calibration_group_size,
      n_cells = length(sceptre_object@cells_in_use)
    )
    cat(crayon::green(" \u2713\n"))
  }
  if(sceptre_object@resampling_mechanism == "asymptotic_normality"){
    synthetic_idxs <- integer(0)
  }
  gc() |> invisible()

  # initialize the args to pass
  args_to_pass <- list(
    response_matrix = get_response_matrix(sceptre_object),
    grna_assignments = sceptre_object@grna_assignments,
    covariate_matrix = sceptre_object@covariate_matrix,
    response_grna_group_pairs = response_grna_group_pairs |> dplyr::filter(pass_qc),
    output_amount = output_amount,
    resampling_mechanism = sceptre_object@resampling_mechanism,
    resampling_approximation = sceptre_object@resampling_approximation,
    B1 = sceptre_object@B1, B2 = sceptre_object@B2,
    B3 = sceptre_object@B3, calibration_check = calibration_check,
    control_group_complement = sceptre_object@control_group_complement,
    n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
    n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
    side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
    response_precomputations = sceptre_object@response_precomputations,
    response_regression_method = sceptre_object@response_regression_method,
    cells_in_use = sceptre_object@cells_in_use, print_progress = print_progress,
    parallel = parallel, n_processors = n_processors, log_dir = log_dir,
    analysis_type = analysis_type
  )

  # run the method
  out <- if (sceptre_object@resampling_mechanism %in% c("permutations", "asymptotic_normality")) {
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


apply_grouping_to_result <- function(result, sceptre_object, is_calibration_check = FALSE) {
  grna_integration_strategy <- sceptre_object@grna_integration_strategy
  if (grna_integration_strategy == "union") {
    new_result <- result |> dplyr::rename("grna_target" = "grna_group")
  }
  if (grna_integration_strategy %in% c("singleton", "bonferroni")) {
    grna_target_data_frame <- sceptre_object@grna_target_data_frame |>
      dplyr::select(grna_id, grna_target)
    new_result <- result |>
      dplyr::rename("grna_id" = "grna_group") |>
      dplyr::left_join(grna_target_data_frame, by = "grna_id") |>
      dplyr::relocate(response_id, grna_id, grna_target)
    if (grna_integration_strategy == "bonferroni" && !is_calibration_check) {
      new_result <- new_result |>
        dplyr::group_by(response_id, grna_target) |>
        dplyr::group_modify(.f = function(tbl, key) {
          if (!any(tbl$pass_qc)) {
            n_nonzero_trt_out <- max(tbl$n_nonzero_trt)
            n_nonzero_cntrl_out <- max(tbl$n_nonzero_cntrl)
            log_2_fold_change_out <- p_value_out <- NA
            pass_qc_out <- FALSE
          } else {
            min_p_idx <- which.min(tbl$p_value)
            n_nonzero_trt_out <- tbl$n_nonzero_trt[min_p_idx]
            n_nonzero_cntrl_out <- tbl$n_nonzero_cntrl[min_p_idx]
            log_2_fold_change_out <- tbl$log_2_fold_change[min_p_idx]
            p_value_out <- min(sum(tbl$pass_qc) * tbl$p_value[min_p_idx], 1)
            pass_qc_out <- TRUE
          }
          data.frame(
            p_value = p_value_out,
            n_nonzero_trt = n_nonzero_trt_out,
            n_nonzero_cntrl = n_nonzero_cntrl_out,
            log_2_fold_change = log_2_fold_change_out,
            pass_qc = pass_qc_out
          )
        }) |>
        dplyr::ungroup()
    }
    data.table::setorderv(new_result, c("p_value", "response_id"), na.last = TRUE)
  }
  return(new_result)
}


#' Get result
#'
#' `get_result()` returns a data frame containing the result of a calibration check, power analysis, or discovery analysis. See \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_write_outputs_to_directory}{Section 8 of the introductory chapter in the manual} for more information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param analysis a string indicating the name of the analysis whose results we are querying, one of `"run_calibration_check"`, `"run_power_check"`, or `"run_discovery_analysis"`.
#'
#' @returns a data frame containing the results of the analysis
#' @export
#' @examples
#' library(sceptredata)
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' pc_result <- sceptre_object |>
#'   set_analysis_parameters(
#'     side = "left",
#'     resampling_mechanism = "permutations",
#'     positive_control_pairs = positive_control_pairs
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc() |>
#'   run_power_check() |>
#'   get_result("run_power_check")
get_result <- function(sceptre_object, analysis) {
  if (!(analysis %in% c("run_calibration_check", "run_power_check", "run_discovery_analysis"))) {
    stop("`analysis` must be one of `run_calibration_check`, `run_power_check`, or `run_discovery_analysis`.")
  }
  if (!sceptre_object@functs_called[[analysis]]) stop(analysis, " has not yet been run.")
  field_to_extract <- switch(EXPR = analysis,
    run_calibration_check = "calibration_result",
    run_power_check = "power_result",
    run_discovery_analysis = "discovery_result"
  )
  out <- methods::slot(sceptre_object, field_to_extract)
  return(out)
}


#' Write outputs to directory
#'
#' `write_outputs_to_directory()` writes the outputs of a `sceptre` analysis to a directory on disk. See \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_write_outputs_to_directory}{Section 8 of the introductory chapter in the manual} for more information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param directory a string giving the file path to a directory on disk in which to write the results
#'
#' @return the value NULL
#' @export
#' @examples
#' library(sceptredata)
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_cis_pairs(sceptre_object,
#'   positive_control_pairs = positive_control_pairs,
#'   distance_threshold = 5e6
#' )
#' sceptre_object |>
#'   set_analysis_parameters(
#'     side = "left",
#'     resampling_mechanism = "permutations",
#'     discovery_pairs = discovery_pairs,
#'     positive_control_pairs = positive_control_pairs
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc() |>
#'   run_calibration_check(
#'     parallel = TRUE,
#'     n_processors = 2
#'   ) |>
#'   run_power_check() |>
#'   run_discovery_analysis(
#'     parallel = TRUE,
#'     n_processors = 2
#'   ) |>
#'   write_outputs_to_directory(paste0(tempdir(), "/sceptre_outputs"))
#' # files written to "sceptre_outputs" in tempdir()
write_outputs_to_directory <- function(sceptre_object, directory) {
  # 0. create directory
  if (!dir.exists(directory)) dir.create(path = directory, recursive = TRUE)
  fs <- paste0(directory, "/", c(
    "analysis_summary.txt", "plot_assign_grnas.png", "plot_run_qc.png",
    "plot_run_calibration_check.png", "plot_run_power_check.png", "plot_run_discovery_analysis.png",
    "results_run_calibration_check.rds", "results_run_power_check.rds", "results_run_discovery_analysis.rds"
  ))
  for (f in fs) if (file.exists(f)) file.remove(f)

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

  # 4. save gRNA-to-cell assignments
  grna_assignments <- get_grna_assignments(sceptre_object)
  saveRDS(object = grna_assignments, file = paste0(directory, "/grna_assignment_matrix.rds"))
  return(NULL)
}


skip_assign_grnas_and_run_qc <- function(sceptre_object, parallel, n_processors) {
  functs_run <- sceptre_object@functs_called
  if (functs_run[["import_data"]] && functs_run[["set_analysis_parameters"]]) {
    if (!functs_run[["assign_grnas"]] && !functs_run[["run_qc"]]) {
      cat(crayon::red("Note: Automatically running `assign_grnas()` and `run_qc()` with default arguments.\n\n"))
      sceptre_object <- assign_grnas(sceptre_object, parallel = parallel, n_processors = n_processors) |> run_qc() # advance by two
    }
    if (functs_run[["assign_grnas"]] && !functs_run[["run_qc"]]) {
      cat(crayon::red("Note: Automatically running `run_qc()` with default arguments.\n\n"))
      sceptre_object <- sceptre_object |> run_qc() # advance by one
    }
  }
  return(sceptre_object)
}
