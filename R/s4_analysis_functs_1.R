#' Set analysis parameters
#'
#' `set_analysis_parameters()` sets the analysis parameters that control how the
#' statistical analysis is to be conducted. See
#' \href{https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html}{Chapter 2 of the manual}
#' for more detailed information about this function.
#'
#' @note Every argument to this function is optional, but typically, users want
#' to specify `discovery_pairs` at minimum.
#'
#' @param sceptre_object a `sceptre_object`
#' @param discovery_pairs (optional) a data frame with columns `grna_target` and
#' `response_id` specifying the discovery pairs to analyze
#' @param positive_control_pairs (optional) a data frame with columns
#' `grna_target` and `response_id` specifying the positive control pairs to
#' analyze
#' @param formula_object (optional) a formula object specifying how to adjust
#' for the covariates in the model
#' @param side (optional; default `"both"`) the sidedness of the test, one of
#' `"left"`, `"right"`, or `"both"`
#' @param grna_integration_strategy (optional; default `"union"`) a string
#' specifying the gRNA integration strategy, either `"singleton"`, `"union"`, or
#' `"bonferroni"`
#' @param resampling_approximation (optional; default `"skew_normal"`) a string
#' indicating the resampling approximation to make to the null distribution of
#' test statistics, either `"skew_normal"` or `"no_approximation"`
#' @param control_group (optional) a string specifying the control group to use
#' in the differential expression analysis, either `"complement"` or
#' `"nt_cells"`
#' @param resampling_mechanism (optional) a string specifying the resampling
#' mechanism to use, either `"permutations"` or `"crt"`
#' @param multiple_testing_method (optional; default `"BH"`) a string specifying
#' the multiple testing correction method to use; see `p.adjust.methods` for
#' options
#' @param multiple_testing_alpha (optional; default `0.1`) a numeric specifying
#' the nominal level of the multiple testing correction method
#'
#' @return an updated `sceptre_object` in which the analysis parameters have
#' been set
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
#' # set analysis parameters
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_cis_pairs(sceptre_object,
#'   positive_control_pairs = positive_control_pairs,
#'   distance_threshold = 5e6
#' )
#' sceptre_object <- sceptre_object |>
#'   set_analysis_parameters(
#'     discovery_pairs = discovery_pairs,
#'     positive_control_pairs = positive_control_pairs,
#'     side = "left"
#'   )
set_analysis_parameters <- function(
    sceptre_object,
    discovery_pairs = data.frame(
        grna_target = character(0),
        response_id = character(0)
    ),
    positive_control_pairs = data.frame(
        grna_target = character(0),
        response_id = character(0)
    ),
    side = "both",
    grna_integration_strategy = "union",
    formula_object = "default",
    resampling_approximation = "skew_normal",
    control_group = "default",
    resampling_mechanism = "default",
    multiple_testing_method = "BH",
    multiple_testing_alpha = 0.1
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    # Requested inputs survive calibration's changes to the working budget.
    requested <- as.list(environment())
    requested$sceptre_object <- NULL
    # verify that function called in correct order
    sceptre_object <- perform_status_check_and_update(
        sceptre_object,
        "set_analysis_parameters"
    )

    # handle default arguments
    if (!sceptre_object@data@low_moi) {
        control_group <- "complement"
        if (resampling_mechanism == "default") resampling_mechanism <- "crt"
    }
    if (sceptre_object@data@low_moi) {
        if (control_group == "default") {
            control_group <- "nt_cells"
        }
        if (resampling_mechanism == "default") {
            resampling_mechanism <- "permutations"
        }
    }
    if (identical(formula_object, "default")) {
        formula_object <- auto_construct_formula_object(
            cell_covariates = sceptre_object@data@covariate_data_frame,
            include_grna_covariates = !sceptre_object@data@low_moi
        )
    }
    B1 <- 499L
    if (resampling_approximation == "skew_normal") {
        B2 <- 4999L
        B3 <- if (resampling_mechanism == "permutations") 24999L else 0L
    } else if (resampling_approximation == "no_approximation") {
        B2 <- 0L # no curve fitting; thus, B2 = 0L
        B3 <- 0L # to be updated in the run_qc step
    }

    # check inputs
    check_set_analysis_parameters(
        sceptre_object = sceptre_object,
        formula_object = formula_object,
        response_grna_target_pairs_list = list(
            discovery_pairs = discovery_pairs,
            positive_control_pairs = positive_control_pairs
        ),
        control_group = control_group,
        resampling_mechanism = resampling_mechanism,
        side = side,
        low_moi = sceptre_object@data@low_moi,
        grna_integration_strategy = grna_integration_strategy,
        resampling_approximation = resampling_approximation
    ) |>
        invisible()

    # determine whether to reset response precomputations
    design <- convert_covariate_df_to_design_matrix(
        covariate_data_frame = sceptre_object@data@covariate_data_frame,
        formula_object = formula_object
    )
    reset_response_precomps <- response_cache_parameters_changed(
        sceptre_object,
        design,
        control_group == "complement"
    )

    # update uncached fields of the sceptre object
    side_code <- which(side == c("left", "both", "right")) - 2L
    control_group_complement <- control_group == "complement"
    run_permutations <- resampling_mechanism == "permutations"
    sceptre_object@spec@discovery_pairs <- discovery_pairs |>
        dplyr::mutate(
            grna_target = as.character(grna_target),
            response_id = as.character(response_id)
        )
    sceptre_object@spec@positive_control_pairs <- positive_control_pairs |>
        dplyr::mutate(
            grna_target = as.character(grna_target),
            response_id = as.character(response_id)
        )
    sceptre_object@preparation@n_discovery_pairs <- nrow(sceptre_object@spec@discovery_pairs)
    sceptre_object@preparation@n_positive_control_pairs <- nrow(
        sceptre_object@spec@positive_control_pairs
    )
    sceptre_object@spec@formula_object <- formula_object
    sceptre_object@spec@side_code <- side_code
    sceptre_object@spec@resampling_approximation <- resampling_approximation
    sceptre_object@spec@control_group_complement <- control_group_complement
    sceptre_object@spec@run_permutations <- run_permutations
    sceptre_object@spec@B1 <- B1
    sceptre_object@spec@B2 <- B2
    sceptre_object@spec@B3 <- B3
    sceptre_object@spec@multiple_testing_alpha <- multiple_testing_alpha
    sceptre_object@spec@multiple_testing_method <- multiple_testing_method
    sceptre_object@spec@grna_integration_strategy <- grna_integration_strategy
    sceptre_object@preparation@covariate_matrix <- design

    # modify the grna target df and response grna target dfs
    sceptre_object <- update_dfs_based_on_grouping_strategy(sceptre_object)

    # update cached fields
    sceptre_object <- invalidate_response_cache(
        sceptre_object,
        "parameters",
        reset_response_precomps
    )

    sceptre_object@functs_called["set_analysis_parameters"] <- TRUE
    sceptre_object@last_function_called <- "set_analysis_parameters"
    sceptre_object@spec@requested <- requested
    sceptre_object <- invalidate_result_records(sceptre_object, "analysis parameters changed")
    sceptre_object@preparation@stages$design <- list(complete = TRUE,
        identity = digest::digest(sceptre_object@preparation@covariate_matrix, algo = "sha256"))
    sceptre_object@preparation@stages$assignment <- list(complete = FALSE)
    sceptre_object@preparation@stages$cell_qc <- list(complete = FALSE)
    sceptre_object@preparation@stages$pair_qc <- list(complete = FALSE)
    validObject(sceptre_object)
    return(sceptre_object)
}


#' Assign gRNAs to cells
#'
#' `assign_grnas()` performs the gRNA-to-cell assignments. `sceptre` provides
#' three gRNA-to-cell assignment strategies: the mixture method, the
#' thresholding method, and the maximum method. The mixture method assigns gRNAs
#' by thresholding fitted posterior membership probabilities. Next, the
#' thresholding method assigns a gRNA to a cell if the UMI count of the gRNA in
#' the cell is greater than or equal to some integer threshold. Finally, the
#' maximum method assigns the gRNA that accounts for the greatest number of UMIs
#' in a given cell to that cell. The maximum method is available only in low
#' MOI. See
#' \href{https://timothy-barry.github.io/sceptre-book/assign-grnas.html}{Chapter 3 of the manual}
#' for more detailed information about `assign_grnas()`.
#'
#' @note The `"mixture"` method retains the historical fitting algorithm used
#' in `sceptre` version 0.9.0 and later. Owing to an implementation
#' inconsistency, the perturbation parameter is estimated under a multiplicative
#' parameterization but applied additively in the subsequent E-step. The
#' resulting parameter estimates therefore do not maximize the likelihood of
#' the mixture model used to compute the fitted probabilities. The historical
#' implementation is retained because it has performed well in benchmarks. This
#' limitation affects only `"mixture"`;
#' `"maximum"` and `"thresholding"` are unaffected. We plan to replace the
#' historical mixture implementation after a
#' statistically coherent and empirically validated alternative is available.
#' Any change to the public API or high-MOI default will follow that validation.
#' Progress is tracked in
#' [GitHub issue #220](https://github.com/Katsevich-Lab/sceptre/issues/220).
#' See the manual for information about the method-specific additional
#' arguments.
#' @param sceptre_object a `sceptre_object`
#' @param method (optional) a string indicating the method to use to assign the
#' gRNAs to cells, one of `"mixture"`, `"thresholding"`, or `"maximum"`. The
#' default is `"maximum"` in low MOI and `"mixture"` in high MOI.
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
#' @param ... optional method-specific additional arguments
#'
#' @return an updated `sceptre_object` in which the gRNA assignments have been
#' carried out
#' @export
#' @examples
#' data("lowmoi_example_data")
#' # 1. import data, set default analysis parameters
#' sceptre_object <- import_data(
#'   response_matrix = lowmoi_example_data$response_matrix,
#'   grna_matrix = lowmoi_example_data$grna_matrix,
#'   extra_covariates = lowmoi_example_data$extra_covariates,
#'   grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'   moi = "low"
#' ) |> set_analysis_parameters()
#'
#' # 2. assign gRNAs (three different methods)
#' sceptre_object <- sceptre_object |> assign_grnas(method = "thresholding")
#' sceptre_object <- sceptre_object |> assign_grnas(method = "maximum")
#' sceptre_object <- sceptre_object |> assign_grnas(method = "mixture")
assign_grnas <- function(
    sceptre_object,
    method = "default",
    print_progress = TRUE,
    parallel = FALSE,
    n_processors = "auto",
    log_dir = tempdir(),
    ...
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    # verify that function called in correct order
    check_parallel_supported(parallel) |> invisible()
    sceptre_object <- perform_status_check_and_update(
        sceptre_object,
        "assign_grnas"
    )

    # handle the default arguments
    if (identical(method, "default")) {
        method <- if (sceptre_object@data@low_moi) "maximum" else "mixture"
    }
    hyperparameters_default <- if (method == "maximum") {
        list(
            umi_fraction_threshold = 0.8,
            min_grna_n_umis_threshold = 5L
        )
    } else if (method == "thresholding") {
        list(threshold = 5)
    } else if (method == "mixture") {
        list(
            n_em_rep = 5L,
            pi_guess_range = c(1e-5, 0.1),
            g_pert_guess_range = log(c(10, 5000)),
            n_nonzero_cells_cutoff = 10L,
            backup_threshold = 5L,
            probability_threshold = 0.8,
            formula_object = auto_construct_formula_object(
                cell_covariates = sceptre_object@data@covariate_data_frame,
                include_grna_covariates = TRUE
            )
        )
    }
    hyperparameters <- list(...)
    if (length(hyperparameters) == 0L) {
        hyperparameters <- hyperparameters_default
    }
    for (hyperparam_name in names(hyperparameters)) {
        hyperparameters_default[[hyperparam_name]] <- hyperparameters[[
            hyperparam_name
        ]]
    }
    hyperparameters <- hyperparameters_default

    # check inputs
    check_assign_grna_inputs(
        sceptre_object,
        method,
        hyperparameters,
        n_processors
    ) |>
        invisible()

    # determine whether to reset response precomputations
    sceptre_object <- invalidate_response_cache(
        sceptre_object,
        "assignments",
        response_cache_assignments_changed(
            sceptre_object,
            method,
            hyperparameters
        )
    )

    # update uncached fields
    sceptre_object@spec@grna_assignment_method <- method
    sceptre_object@spec@grna_assignment_hyperparameters <- hyperparameters

    # assign the grnas
    sceptre_object <- assign_grnas_to_cells(
        sceptre_object,
        print_progress,
        parallel,
        n_processors,
        log_dir
    )

    sceptre_object <- invalidate_result_records(sceptre_object, "guide assignments changed")
    sceptre_object@functs_called["assign_grnas"] <- TRUE
    sceptre_object@last_function_called <- "assign_grnas"
    sceptre_object@preparation@stages$assignment <- list(complete = TRUE,
        coordinates = "source", identity = digest::digest(
            sceptre_object@preparation@initial_grna_assignment_list, algo = "sha256"))
    sceptre_object@preparation@stages$cell_qc <- list(complete = FALSE)
    sceptre_object@preparation@stages$pair_qc <- list(complete = FALSE)
    validObject(sceptre_object)
    return(sceptre_object)
}


#' Run QC
#'
#' `run_qc()` runs cellwise and pairwise QC on the data. Cellwise QC involves
#' filtering cells on the covariates `response_n_nonzero`, `response_n_umis`,
#' and `response_p_mito`. In low-MOI we additionally remove cells that contain
#' zero or multiple gRNAs. Next, pairwise QC involves filtering out
#' target-response pairs whose data are too sparse to be analyzed reliably. In
#' this context we define the “number of nonzero treatment cells” (resp., the
#' “number of nonzero control cells”) as the number of cells in the treatment
#' group (resp., control group) that contain nonzero expression of the response.
#' (We sometimes use the shorthand `n_nonzero_trt` and `n_nonzero_cntrl` to
#' refer to the number of nonzero treatment cells and control cells,
#' respectively.) Pairwise QC involves filtering target-response pairs on
#' `n_nonzero_trt` and `n_nonzero_cntrl`. See
#' \href{https://timothy-barry.github.io/sceptre-book/run-qc.html}{Chapter 4 of the manual}
#' for more detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param n_nonzero_trt_thresh (optional; default `7L`) an integer specifying
#' the number of nonzero *treatment* cells a pair must contain for it to be
#' retained
#' @param n_nonzero_cntrl_thresh (optional; default `7L`) an integer specifying
#' the number of nonzero *control* cells a pair must contain for it to be
#' retained
#' @param response_n_umis_range (optional; default `c(0.01, 0.99)`) a length-two
#' vector of percentiles specifying the location at which to clip the left and
#' right tails of the `response_n_umis` distribution
#' @param response_n_nonzero_range (optional; default `c(0.01, 0.99)`) a
#' length-two vector of percentiles specifying the location at which to clip the
#' left and right tails of the `response_n_nonzero` distribution
#' @param p_mito_threshold (optional; default `0.2`) a numeric value specifying
#' the location at which to clip the right tail of the `response_p_mito`
#' distribution
#' @param additional_cells_to_remove (optional) a vector of integer indices
#' specifying additional cells to remove
#'
#' @return an updated `sceptre_object` in which cellwise and pairwise QC have
#' been applied
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
#' discovery_pairs <- construct_cis_pairs(sceptre_object,
#'   positive_control_pairs = positive_control_pairs,
#'   distance_threshold = 5e6
#' )
#' sceptre_object <- sceptre_object |>
#'   set_analysis_parameters(
#'     discovery_pairs = discovery_pairs,
#'     positive_control_pairs = positive_control_pairs,
#'     side = "left"
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc()
run_qc <- function(
    sceptre_object,
    n_nonzero_trt_thresh = 7L,
    n_nonzero_cntrl_thresh = 7L,
    response_n_umis_range = c(0.01, 0.99),
    response_n_nonzero_range = c(0.01, 0.99),
    p_mito_threshold = 0.2,
    additional_cells_to_remove = integer()
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    run_qc_pt_1(
        sceptre_object,
        n_nonzero_trt_thresh,
        n_nonzero_cntrl_thresh,
        response_n_umis_range,
        response_n_nonzero_range,
        p_mito_threshold,
        additional_cells_to_remove
    ) |>
        run_qc_pt_2()
}


#' Run cellwise QC, the phase also used by nuclear Nextflow tasks
#' @param sceptre_object Object with initial guide assignments and cell covariates.
#' @param n_nonzero_trt_thresh,n_nonzero_cntrl_thresh Pairwise QC cutoffs to save.
#' @param response_n_umis_range,response_n_nonzero_range Quantile bounds for cell counts.
#' @param p_mito_threshold Maximum mitochondrial fraction, when that summary exists.
#' @param additional_cells_to_remove Original matrix column numbers to remove, starting at 1.
#' @return Object with QC status, thresholds, removal counts, original columns of
#'   retained cells and guide assignments renumbered for those cells. Clear saved
#'   response fits if retained cells changed; keep the full original design matrix.
#'   Warn about large low-MOI assignment exclusions and error if no cells remain.
#' @details Nextflow calls this phase separately before distributing pairwise QC.
#' @noRd
run_qc_pt_1 <- function(
    sceptre_object,
    n_nonzero_trt_thresh = 7L,
    n_nonzero_cntrl_thresh = 7L,
    response_n_umis_range = c(0.01, 0.99),
    response_n_nonzero_range = c(0.01, 0.99),
    p_mito_threshold = 0.2,
    additional_cells_to_remove = integer()
) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    # cellwise start
    # verify that function called in correct order
    sceptre_object <- perform_status_check_and_update(sceptre_object, "run_qc")

    # check inputs
    check_run_qc_inputs(
        n_nonzero_trt_thresh,
        n_nonzero_cntrl_thresh,
        response_n_umis_range,
        response_n_nonzero_range,
        sceptre_object@preparation@initial_grna_assignment_list
    ) |>
        invisible()

    # 5 update fields of sceptre_object
    sceptre_object@spec@cellwise_qc_thresholds <- list(
        response_n_umis_range = response_n_umis_range,
        response_n_nonzero_range = response_n_nonzero_range,
        p_mito_threshold = p_mito_threshold
    )

    # obtain previous cells_in_use for caching purposes
    current_cells_in_use <- sceptre_object@preparation@cells_in_use

    # update uncached fields of the sceptre object
    sceptre_object@spec@n_nonzero_trt_thresh <- as.integer(n_nonzero_trt_thresh)
    sceptre_object@spec@n_nonzero_cntrl_thresh <- as.integer(n_nonzero_cntrl_thresh)

    # determine the cells to retain after cellwise qc
    sceptre_object <- determine_cells_to_retain(
        sceptre_object,
        response_n_umis_range,
        response_n_nonzero_range,
        p_mito_threshold,
        additional_cells_to_remove
    )
    # informative error/warning messages
    frac_cells_removed_2p_grna <- sceptre_object@preparation@cell_removal_metrics[[
        "n_cells_rm_zero_twoplus_grnas"
    ]] /
        nrow(sceptre_object@data@covariate_data_frame)
    if (frac_cells_removed_2p_grna >= 0.5 && sceptre_object@data@low_moi) {
        warning(
            "A large percentage of cells (",
            round(100 * frac_cells_removed_2p_grna, 3),
            "%) contain two or more gRNAs. Consider switching from a low-MOI to a high-MOI analysis."
        )
    }
    if (length(sceptre_object@preparation@cells_in_use) == 0) {
        stop(
            "All cells were removed as part of cellwise QC. Consider relaxing the cellwise QC thresholds (`response_n_umis_range`, `response_n_nonzero_range`, `p_mito_threshold`).",
            if (sceptre_object@data@low_moi) {
                " Additionally, consider switching from a low-MOI to high-MOI analysis."
            } else {
                ""
            }
        )
    }

    # determine whether to reset response precomputation
    sceptre_object <- invalidate_response_cache(
        sceptre_object,
        "analysis_cells",
        !identical(current_cells_in_use, sceptre_object@preparation@cells_in_use)
    )

    # update the grna assignments given the cellwise qc
    sceptre_object <- update_grna_assignments_given_qc(sceptre_object)
    sceptre_object@functs_called["run_qc"] <- FALSE
    sceptre_object <- invalidate_result_records(sceptre_object, "QC selection or eligibility changed")
    sceptre_object@preparation@stages$cell_qc <- list(complete = TRUE,
        cells = sceptre_object@preparation@cells_in_use)
    sceptre_object@preparation@stages$pair_qc <- list(complete = FALSE)
    return(sceptre_object)
}


#' Run pairwise QC and set the empirical resampling size when needed
#' @param sceptre_object Object after cellwise QC, with updated guide assignments.
#' @return Object with pair counts, pass_qc flags, NT summaries and passing-pair
#'   counts. With no_approximation, also calculate and check B3. Count data through
#'   storage helpers; random test draws are generated later. Nuclear Nextflow
#'   tasks can call this phase separately.
#' @noRd
run_qc_pt_2 <- function(sceptre_object) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    # pairwise start
    # compute (i) the NT M matrix, (ii), n nonzero total vector, (iii)
    # n_nonzero_trt, and (iv) n_nonzero_cntrl vectors
    sceptre_object <- compute_pairwise_qc_information(sceptre_object)

    # compute the number of discovery pairs and (if applicable) pc pairs
    # passing qc
    sceptre_object <- compute_qc_metrics(sceptre_object)

    # update B3, the number of resamples to draw, if resampling_approximation is
    # no_approximation
    if (
        sceptre_object@spec@resampling_approximation == "no_approximation" &&
            (length(sceptre_object@preparation@n_ok_discovery_pairs) >= 1L ||
                length(sceptre_object@preparation@n_ok_positive_control_pairs) >= 1L)
    ) {
        sceptre_object@spec@B3 <- compute_empirical_budget(
            max(
                sceptre_object@preparation@n_ok_discovery_pairs,
                sceptre_object@preparation@n_ok_positive_control_pairs
            ),
            sceptre_object@spec@multiple_testing_alpha,
            sceptre_object@spec@side_code,
            sceptre_object@spec@B1,
            sceptre_object@spec@B2
        )
    }

    sceptre_object@preparation@stages$pair_qc <- list(complete = TRUE)
    sceptre_object@functs_called["run_qc"] <- TRUE
    sceptre_object@last_function_called <- "run_qc"
    validObject(sceptre_object)
    return(sceptre_object)
}
