#' Run the requested gene-guide tests that passed QC
#' @param sceptre_object Sceptre object supplying data and settings. This function
#'   returns results separately; the public workflow saves them into the object.
#' @param response_grna_group_pairs All requested pairs, including failed-QC rows.
#' @param analysis_type Calibration, power or discovery analysis name.
#' @param output_amount,print_progress,parallel,n_processors,log_dir Validated
#'   output, progress, worker-count and logging settings from the public function.
#' @param plan Optional validated transient run plan.
#' @return A list with the result table and reusable response fits in
#'   `response_precomputations`. Failed-QC requests are included in the table.
#' @details Prepare the cell selections and settings, then choose the permutation
#'   or CRT loop. Permutations share one full set of random draws across the run;
#'   CRT generates a full set per guide group. Keep that sharing and draw order:
#'   regenerating draws for each pair could change results. Guide integration and
#'   multiple-testing correction happen later, in the caller.
#' @noRd
analyze_association_pairs <- function(
    sceptre_object,
    response_grna_group_pairs,
    analysis_type,
    output_amount,
    print_progress,
    parallel,
    n_processors,
    log_dir,
    plan = NULL
) {
    if (is.null(plan)) plan <- make_run_plan(sceptre_object, analysis_type, response_grna_group_pairs)
    validate_run_plan(sceptre_object, plan)
    if (!identical(plan$analysis, analysis_record_name(analysis_type)) ||
        !identical(plan$pairs, response_grna_group_pairs)) {
        stop("Association inputs do not match the run plan.")
    }
    context <- prepare_analysis_context(
        sceptre_object,
        analysis_type,
        pairs = response_grna_group_pairs |> dplyr::filter(pass_qc),
        plan = plan
    )
    execution <- list(
        output_amount = output_amount,
        print_progress = print_progress,
        parallel = parallel,
        n_processors = n_processors,
        log_dir = log_dir
    )
    if (sceptre_object@spec@run_permutations) {
        message("Generating permutation resamples.")
        synthetic_idxs <- get_synthetic_permutation_idxs(
            grna_assignments = sceptre_object@preparation@analysis_groups,
            B = context$test_options$total_resamples,
            calibration_check = context$index$calibration_check,
            control_group_complement = context$index$control_group_complement,
            calibration_group_size = if (is.null(plan$calibration_group_size)) 1L else plan$calibration_group_size,
            n_cells = length(context$index$analysis_to_source)
        )
        message(crayon::green(" \u2713"))
    }
    gc() |> invisible()
    out <- if (sceptre_object@spec@run_permutations) {
        run_permutation_analysis(context, execution, synthetic_idxs)
    } else {
        run_crt_analysis(context, execution)
    }
    if (!context$shared_response_fit) {
        # Group-specific NT fits are transient. Preserve an existing shared
        # calibration cache with its own identity; do not relabel it all-QC.
        out$response_precomputations <- sceptre_object@fits@entries
    }
    attr(out$response_precomputations, "fit_identity") <- if (context$shared_response_fit) {
        context$fit_identity
    } else sceptre_object@fits@identity
    out$result <- assemble_pair_results(out$result, response_grna_group_pairs)
    out
}
