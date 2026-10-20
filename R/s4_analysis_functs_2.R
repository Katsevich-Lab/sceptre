#' Run calibration check
#'
#' `run_calibration_check()` runs the calibration check. The calibration check
#' involves applying sceptre to analyze negative control target-response pairs
#' --- pairs for which we know there is no association between the target and
#' response --- to ensure control of the false discovery rate. The calibration
#' check enables us to verify that the discovery set that sceptre ultimately
#' produces is not contaminated by excess false positives. See
#' \href{https://timothy-barry.github.io/sceptre-book/run-calibration-check.html}{Chapter 5 of the manual}
#' for more detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param output_amount (optional; default `1`) an integer taking values 1, 2,
#' or 3 specifying the amount of information to return. `1` returns the least
#' amount of information and `3` the most.
#' @param n_calibration_pairs (optional) the number of negative control pairs to
#' construct and test for association
#' @param calibration_group_size (optional) the number of negative control gRNAs
#' to randomly assemble to form each negative control target
#' @param print_progress (optional; default `TRUE`) a logical indicating whether
#' to print progress updates
#' @param parallel (optional; default `FALSE`) a logical indicating whether to
#' run the function in parallel. `parallel = TRUE` is recommended only on Mac;
#' it is not supported on Windows and may behave unreliably on Linux clusters.
#' @param n_processors (optional; default `"auto"`) an integer specifying the
#' number of processors to use if `parallel` is set to `TRUE`. The default,
#' `"auto"`, uses half the physical cores. The fraction may be tuned via the
#' `parallelly.availableCores.fraction` R option.
#' @param log_dir (optional; default `tempdir()`) a string indicating the
#' directory in which to write the log files (ignored if `parallel = FALSE`)
#' @return an updated `sceptre_object` in which the calibration check has been
#' carried out
#'
#' @export
#' @examples
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
#'     calibration_group_size = 2
#'   )
run_calibration_check <- function(
    sceptre_object,
    n_calibration_pairs = NULL,
    calibration_group_size = NULL,
    print_progress = TRUE,
    parallel = FALSE,
    n_processors = "auto",
    log_dir = tempdir(),
    output_amount = 1
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    sceptre_object <- invalidate_other_inference_results(sceptre_object)
    check_parallel_supported(parallel) |> invisible()
    prepared <- run_calibration_check_pt_1(
        sceptre_object,
        n_calibration_pairs = n_calibration_pairs,
        calibration_group_size = calibration_group_size,
        parallel = parallel,
        n_processors = n_processors
    )
    run_calibration_check_pt_2(
        prepared$object,
        prepared$plan,
        output_amount = output_amount,
        print_progress = print_progress,
        parallel = parallel,
        n_processors = n_processors,
        log_dir = log_dir
    )
}


#' Prepare calibration settings and choose negative-control pairs
#'
#' Complete missing assignment/QC steps, record workflow status, choose defaults
#' and check calibration sizes. For empirical p-values, calculate B3 for the
#' requested calibration size before sampling the pairs.
#' @return List with the object and a transient validated run plan. The plan
#'   owns
#'   calibration settings and eligible pairs; no result is committed yet.
#' @noRd
run_calibration_check_pt_1 <- function(
    sceptre_object,
    n_calibration_pairs = NULL,
    calibration_group_size = NULL,
    parallel = FALSE,
    n_processors = "auto"
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    sceptre_object <- skip_assign_grnas_and_run_qc(
        sceptre_object,
        parallel,
        n_processors
    )
    require_stage(sceptre_object, "run_qc")
    if (!parallel) {
        message(crayon::red(
            "Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`."
        ))
    }
    if (
        pair_count(sceptre_object, "discovery") == 0L &&
            (is.null(n_calibration_pairs) || is.null(calibration_group_size))
    ) {
        stop(
            "`calibration_group_size` and `n_calibration_pairs` must be supplied when discovery pairs are not specified."
        )
    }
    if (is.null(calibration_group_size)) {
        calibration_group_size <-
            compute_calibration_group_size(guide_group_metadata(
                sceptre_object
            ))
    }
    if (is.null(n_calibration_pairs)) {
        n_calibration_pairs <- pair_count(
            sceptre_object,
            "discovery",
            passing = TRUE
        )
    }
    check_calibration_check_inputs(
        sceptre_object,
        n_calibration_pairs,
        n_processors,
        calibration_group_size
    ) |>
        invisible()
    budgets <- run_resampling_budget(
        sceptre_object,
        "run_calibration_check",
        n_calibration_pairs
    )
    pairs <- construct_negative_control_pairs(
        sceptre_object,
        n_calibration_pairs = n_calibration_pairs,
        calibration_group_size = calibration_group_size
    )
    plan <- make_run_plan(
        sceptre_object,
        "run_calibration_check",
        pairs,
        n_calibration_pairs = n_calibration_pairs,
        calibration_group_size = calibration_group_size,
        budgets = budgets
    )
    list(object = sceptre_object, plan = plan)
}


#' Test previously prepared calibration pairs and save their results
#'
#' Use the pairs chosen in phase one and the current cells and settings.
#' Analysis
#' helpers handle random draws, progress and serial/parallel work.
#' @param plan Validated transient plan produced by calibration preparation.
#' @return Object with a complete-family calibration record and reusable fits,
#'   committed only after successful execution and result finalization.
#' @noRd
run_calibration_check_pt_2 <- function(
    sceptre_object,
    plan,
    output_amount = 1,
    print_progress = TRUE,
    parallel = FALSE,
    n_processors = "auto",
    log_dir = tempdir()
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    require_stage(sceptre_object, "run_qc")
    plan <- validate_run_plan(sceptre_object, plan)
    if (!identical(plan@analysis, "run_calibration_check")) {
        stop("Expected a calibration run plan.")
    }
    out <- analyze_association_pairs(
        sceptre_object,
        response_grna_group_pairs = plan@pairs,
        analysis_type = "calibration_check",
        output_amount = output_amount,
        print_progress = print_progress,
        parallel = parallel,
        n_processors = n_processors,
        log_dir = log_dir,
        plan = plan
    )
    result <- process_calibration_result(out$result, sceptre_object)
    sceptre_object <- commit_result_record(
        sceptre_object,
        "run_calibration_check",
        result,
        plan@pairs,
        output_amount,
        plan = plan
    )
    commit_analysis_fits(sceptre_object, out)
}


#' Run power check
#'
#' `run_power_check()` runs the power check. The power check entails applying
#' `sceptre` to analyze positive control pairs --- pairs for which we know there
#' is an association between the target and the response --- to ensure that
#' `sceptre` is capable of detecting true associations. See
#' \href{https://timothy-barry.github.io/sceptre-book/run-power-check-and-discovery-analysis.html}{Chapter 6 of the manual}
#' for more detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @inheritParams run_calibration_check
#' @return an updated `sceptre_object` in which the power check has been carried
#' out
#'
#' @export
#' @examples
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
run_power_check <- function(
    sceptre_object,
    output_amount = 1,
    print_progress = TRUE,
    parallel = FALSE,
    n_processors = "auto",
    log_dir = tempdir()
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    sceptre_object <- invalidate_other_inference_results(sceptre_object)
    check_parallel_supported(parallel) |> invisible()
    sceptre_object <- skip_assign_grnas_and_run_qc(
        sceptre_object,
        parallel,
        n_processors
    )
    require_stage(sceptre_object, "run_qc")
    if (!parallel) {
        message(crayon::red(
            "Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`."
        ))
    }
    response_grna_group_pairs <-
        sceptre_object@preparation@positive_control_pairs_with_info
    check_discovery_analysis_inputs(
        response_grna_group_pairs = response_grna_group_pairs,
        control_group_complement = sceptre_object@spec@control_group_complement,
        grna_target_data_frame = guide_group_metadata(sceptre_object),
        pc_analysis = TRUE,
        calibration_result = current_result(
            sceptre_object,
            "run_calibration_check",
            required = FALSE
        ),
        n_ok_pairs = pair_count(
            sceptre_object,
            "positive_control",
            passing = TRUE
        ),
        n_processors = n_processors
    ) |>
        invisible()
    plan <- make_run_plan(
        sceptre_object,
        "run_power_check",
        response_grna_group_pairs
    )
    out <- analyze_association_pairs(
        sceptre_object,
        response_grna_group_pairs = response_grna_group_pairs,
        analysis_type = "power_check",
        output_amount = output_amount,
        print_progress = print_progress,
        parallel = parallel,
        n_processors = n_processors,
        log_dir = log_dir,
        plan = plan
    )
    result <- apply_grouping_to_result(out$result, sceptre_object)
    sceptre_object <- commit_result_record(
        sceptre_object,
        "run_power_check",
        result,
        plan@pairs,
        output_amount,
        plan = plan
    )
    commit_analysis_fits(sceptre_object, out)
}


#' Run discovery analysis
#'
#' `run_discovery_analysis()` runs the discovery analysis. The discovery
#' analysis involves applying `sceptre` to analyze discovery pairs, or
#' target-response pairs whose association status we do not know but seek to
#' learn. Identifying associations among the discovery pairs is the primary
#' objective of the single-cell CRISPR screen analysis. See
#' \href{https://timothy-barry.github.io/sceptre-book/run-power-check-and-discovery-analysis.html}{Chapter 6 of the manual}
#' for more detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @inheritParams run_calibration_check
#' @return an updated `sceptre_object` in which the discovery analysis has been
#' carried out
#'
#' @export
#' @examples
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
#'   run_discovery_analysis()
run_discovery_analysis <- function(
    sceptre_object,
    output_amount = 1,
    print_progress = TRUE,
    parallel = FALSE,
    n_processors = "auto",
    log_dir = tempdir()
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    sceptre_object <- invalidate_other_inference_results(sceptre_object)
    check_parallel_supported(parallel) |> invisible()
    sceptre_object <- skip_assign_grnas_and_run_qc(
        sceptre_object,
        parallel,
        n_processors
    )
    require_stage(sceptre_object, "run_qc")
    if (!parallel) {
        message(crayon::red(
            "Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`."
        ))
    }
    response_grna_group_pairs <-
        sceptre_object@preparation@discovery_pairs_with_info
    check_discovery_analysis_inputs(
        response_grna_group_pairs = response_grna_group_pairs,
        control_group_complement = sceptre_object@spec@control_group_complement,
        grna_target_data_frame = guide_group_metadata(sceptre_object),
        pc_analysis = FALSE,
        calibration_result = current_result(
            sceptre_object,
            "run_calibration_check",
            required = FALSE
        ),
        n_ok_pairs = pair_count(sceptre_object, "discovery", passing = TRUE),
        n_processors = n_processors
    ) |>
        invisible()
    plan <- make_run_plan(
        sceptre_object,
        "run_discovery_analysis",
        response_grna_group_pairs
    )
    out <- analyze_association_pairs(
        sceptre_object,
        response_grna_group_pairs = response_grna_group_pairs,
        analysis_type = "discovery_analysis",
        output_amount = output_amount,
        print_progress = print_progress,
        parallel = parallel,
        n_processors = n_processors,
        log_dir = log_dir,
        plan = plan
    )
    result <- process_discovery_result(out$result, sceptre_object)
    sceptre_object <- commit_result_record(
        sceptre_object,
        "run_discovery_analysis",
        result,
        plan@pairs,
        output_amount,
        plan = plan
    )
    commit_analysis_fits(sceptre_object, out)
}


#' Get result
#'
#' `get_result()` returns a data frame containing the result of a calibration
#' check, power check, or discovery analysis. We pass as arguments
#' `sceptre_object` and `analysis`, where the latter is a string indicating the
#' function whose results we are querying. The output is a data frame, the rows
#' of which correspond to target-response pairs, and the columns of which are as
#' follows: `response_id`, `grna_target`, `n_nonzero_trt`, `n_nonzero_cntrl`,
#' `pass_qc` (a `TRUE`/`FALSE` value indicating whether the pair passes pairwise
#' QC), `p_value`, `fold_change`, `se_fold_change` (standard error for the fold
#' change estimate), `log_2_fold_change`, and `significant` (a `TRUE`/`FALSE`
#' value indicating whether the pair is called as significant). The p-value
#' contained within the `p_value` column is a raw (i.e.,
#' non-multiplicity-adjusted) p-value. See
#' \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_write_outputs_to_directory}{Section 8 of the introductory chapter in the manual}
#' for more information about this function.
#'
#'
#' @note If `output_amount` is set to `2` or `3` in `run_calibration_check()`,
#' `run_power_check()`, or `run_discovery_analysis()`, then the result data
#' frame contains additional columns; see
#' \href{https://timothy-barry.github.io/sceptre-book/run-calibration-check.html#sec-run_calibration_check_output_amount}{Chapter 6 in the manual}
#' for more information.
#' @param sceptre_object a `sceptre_object`
#' @param analysis a string indicating the name of the analysis whose results we
#' are querying, one of `"run_calibration_check"`, `"run_power_check"`, or
#' `"run_discovery_analysis"`.
#'
#' @param historical Return the saved table even when its provenance is
#'   historical.
#'   The returned table carries its saved record as an attribute.
#' @returns An owned data frame containing the results of the analysis. Editing
#'   it, including nested historical provenance, does not change the object.
#' @export
#' @examples
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
get_result <- function(sceptre_object, analysis, historical = FALSE) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    if (
        !is.character(analysis) ||
            length(analysis) != 1L ||
            is.na(analysis) ||
            !analysis %in% analysis_names()
    ) {
        stop(
            "`analysis` must be one of `run_calibration_check`, `run_power_check`, or `run_discovery_analysis`."
        )
    }
    if (!isTRUE(historical)) {
        return(owned_result_value(current_result(sceptre_object, analysis)))
    }
    record <- owned_result_value(sceptre_object@results@records[[analysis]])
    result <- record$table
    if (is.null(record) || !is.data.frame(result)) {
        stop(analysis, " has no saved result.")
    }
    attr(result, "sceptre_result_record") <- record
    if (!result_record_is_current(sceptre_object, analysis)) {
        warning(
            "Returning historical results; they do not represent the current ",
            "analysis settings.",
            call. = FALSE
        )
    }
    result
}


#' Write outputs to directory
#'
#' `write_outputs_to_directory()` writes the outputs of a `sceptre` analysis to
#' a directory on disk. `write_outputs_to_directory()` writes several files to
#' the specified directory: a text-based summary of the analysis
#' (`analysis_summary.txt`), the various plots (`*.png`), the calibration check,
#' power check, discovery analysis results (`results_run_calibration_check.rds`,
#' `results_run_power_check.rds`, and `results_run_discovery_analysis.rds`,
#' respectively), and the binary gRNA-to-cell assignment matrix
#' (`grna_assignment_matrix.rds`). See
#' \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_write_outputs_to_directory}{Section 8 of the introductory chapter in the manual}
#' for more information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param directory a string giving the file path to a directory on disk in
#' which to write the results
#'
#' @return No return value; called for its side effect of writing files to
#'   `directory`.
#' @export
#' @examples
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
#'   run_calibration_check() |>
#'   run_power_check() |>
#'   run_discovery_analysis() |>
#'   write_outputs_to_directory(paste0(tempdir(), "/sceptre_outputs"))
#' # files written to "sceptre_outputs" in tempdir()
write_outputs_to_directory <- function(sceptre_object, directory) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    for (analysis in names(result_analysis_fields())) {
        if (result_record_is_current(sceptre_object, analysis)) {
            current_result(sceptre_object, analysis) |> invisible()
        }
    }
    # 0. create directory
    if (!dir.exists(directory)) {
        dir.create(path = directory, recursive = TRUE)
    }
    fs <- paste0(
        directory,
        "/",
        c(
            "analysis_summary.txt",
            "plot_assign_grnas.png",
            "plot_run_qc.png",
            "plot_run_calibration_check.png",
            "plot_run_power_check.png",
            "plot_run_discovery_analysis.png",
            "results_run_calibration_check.rds",
            "results_run_power_check.rds",
            "results_run_discovery_analysis.rds"
        )
    )
    for (f in fs) {
        if (file.exists(f)) file.remove(f)
    }

    # 1. create analysis_summary.txt file
    summary_file_fp <- paste0(directory, "/analysis_summary.txt")
    write_analysis_summary(sceptre_object, summary_file_fp)

    # 2. determine the functions that have been called
    plotting_params <- list(
        device = "png",
        scale = 1.1,
        width = 5,
        height = 4,
        dpi = 330
    )
    functs_run <- workflow_status(sceptre_object)
    functs_run_plots <- functs_run[
        names(functs_run) %in%
            c(
                "assign_grnas",
                "run_qc",
                "run_calibration_check",
                "run_power_check",
                "run_discovery_analysis"
            )
    ]
    functs_run_plots <- c(functs_run_plots, "grna_count_distributions" = TRUE)
    for (funct in names(functs_run_plots)) {
        if (functs_run_plots[[funct]]) {
            plotting_params$plot <- do.call(
                what = paste0("plot_", funct),
                args = list(sceptre_object)
            )
            plotting_params$filename <- paste0(
                directory,
                "/plot_",
                funct,
                ".png"
            )
            do.call(what = ggplot2::ggsave, args = plotting_params)
        }
    }

    # 3. save results
    for (analysis in c(
        "run_calibration_check",
        "run_power_check",
        "run_discovery_analysis"
    )) {
        if (functs_run_plots[[analysis]]) {
            res <- get_result(
                sceptre_object = sceptre_object,
                analysis = analysis
            )
            saveRDS(
                object = res,
                file = paste0(directory, "/results_", analysis, ".rds")
            )
        }
    }

    # 4. save gRNA-to-cell assignments
    grna_assignments <- get_grna_assignments(sceptre_object)
    saveRDS(
        object = grna_assignments,
        file = paste0(directory, "/grna_assignment_matrix.rds")
    )
    return(invisible(NULL))
}


#' Run missing assignment and QC steps with their public defaults
#'
#' Calibration, power and discovery call this after checking parallel support.
#' If data and parameters are ready, run assignment then QC when both are
#' missing,
#' or just QC when assignment is complete. Otherwise return the object
#' unchanged;
#' the later workflow check reports missing prerequisites. Forward parallel and
#' n_processors to assignment. Its messages, randomness and logs, and QC's
#' removal
#' of unusable saved fits, come from those normal public calls.
#' @noRd
skip_assign_grnas_and_run_qc <- function(
    sceptre_object,
    parallel,
    n_processors
) {
    functs_run <- workflow_status(sceptre_object)
    if (
        functs_run[["import_data"]] && functs_run[["set_analysis_parameters"]]
    ) {
        if (!functs_run[["assign_grnas"]] && !functs_run[["run_qc"]]) {
            message(crayon::red(
                "Note: Automatically running `assign_grnas()` and `run_qc()` with default arguments."
            ))
            sceptre_object <- assign_grnas(
                sceptre_object,
                parallel = parallel,
                n_processors = n_processors
            ) |>
                run_qc() # advance by two
        }
        if (functs_run[["assign_grnas"]] && !functs_run[["run_qc"]]) {
            message(crayon::red(
                "Note: Automatically running `run_qc()` with default arguments."
            ))
            sceptre_object <- sceptre_object |> run_qc() # advance by one
        }
    }
    return(sceptre_object)
}
