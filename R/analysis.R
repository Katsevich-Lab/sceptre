# S4-to-component compatibility adapter. Numerical workers receive explicit
# context arrays and execution options, and return values before facade commit.
run_sceptre_analysis_high_level <- function(
    sceptre_object,
    response_grna_group_pairs,
    calibration_check,
    analysis_type,
    output_amount,
    print_progress,
    parallel,
    n_processors,
    log_dir
) {
    total_resamples <- checked_resampling_total(sceptre_object@B1, sceptre_object@B2, sceptre_object@B3)
    context <- prepare_analysis_context(sceptre_object, analysis_type,
        pairs = response_grna_group_pairs |> dplyr::filter(pass_qc))
    # if running permutations, generate the permutation idxs
    if (sceptre_object@run_permutations) {
        message("Generating permutation resamples.")
        synthetic_idxs <- get_synthetic_permutation_idxs(
            grna_assignments = sceptre_object@grna_assignments,
            B = total_resamples,
            calibration_check = calibration_check,
            control_group_complement = sceptre_object@control_group_complement,
            calibration_group_size = sceptre_object@calibration_group_size,
            n_cells = length(sceptre_object@cells_in_use)
        )
        message(crayon::green(" \u2713"))
    }
    gc() |> invisible()

    # initialize the args to pass
    args_to_pass <- list(
        response_matrix = context$response_store,
        grna_assignments = sceptre_object@grna_assignments,
        covariate_matrix = context$design,
        response_grna_group_pairs = context$pairs,
        output_amount = output_amount,
        resampling_approximation = sceptre_object@resampling_approximation,
        B1 = sceptre_object@B1,
        B2 = sceptre_object@B2,
        B3 = sceptre_object@B3,
        calibration_check = calibration_check,
        control_group_complement = sceptre_object@control_group_complement,
        n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
        n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
        side_code = sceptre_object@side_code,
        low_moi = sceptre_object@low_moi,
        response_precomputations = context$caches,
        cells_in_use = sceptre_object@cells_in_use,
        print_progress = print_progress,
        parallel = parallel,
        n_processors = n_processors,
        log_dir = log_dir,
        analysis_type = analysis_type,
        context = context
    )

    # run the method
    out <- if (sceptre_object@run_permutations) {
        args_to_pass$synthetic_idxs <- synthetic_idxs
        do.call(what = "run_perm_test_in_memory", args = args_to_pass)
    } else {
        do.call(what = "run_crt_in_memory_v2", args = args_to_pass)
    }

    out$result <- assemble_pair_results(out$result, response_grna_group_pairs)
    return(out)
}
