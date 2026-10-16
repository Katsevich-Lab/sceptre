#' Run an eligible association family through its prepared permutation/CRT plan
#' @param sceptre_object Sceptre facade; no slot is committed here.
#' @param response_grna_group_pairs Full requested table, including failed QC rows.
#' @param analysis_type Calibration, power or discovery role string.
#' @param output_amount,print_progress,parallel,n_processors,log_dir Valid execution
#'   settings from the public wrapper. Output schema/log bytes remain unchanged.
#' @return List `result` (eligible native results plus failed rows) and compact
#'   `response_precomputations`, for the facade to commit after success.
#' @details Creates one full permutation bank before response traversal; CRT
#'   creates its full banks per guide. Performs the existing initial GC, preserves
#'   RNG/bank scope, and leaves guide integration/family correction to callers.
#' @noRd
analyze_association_pairs <- function(
    sceptre_object, response_grna_group_pairs, analysis_type,
    output_amount, print_progress, parallel, n_processors, log_dir
) {
    context <- prepare_analysis_context(sceptre_object, analysis_type,
        pairs = response_grna_group_pairs |> dplyr::filter(pass_qc))
    execution <- list(output_amount = output_amount, print_progress = print_progress,
        parallel = parallel, n_processors = n_processors, log_dir = log_dir)
    if (sceptre_object@run_permutations) {
        message("Generating permutation resamples.")
        synthetic_idxs <- get_synthetic_permutation_idxs(
            grna_assignments = sceptre_object@grna_assignments,
            B = context$test_options$total_resamples,
            calibration_check = context$index$calibration_check,
            control_group_complement = context$index$control_group_complement,
            calibration_group_size = sceptre_object@calibration_group_size,
            n_cells = length(context$index$analysis_to_source)
        )
        message(crayon::green(" \u2713"))
    }
    gc() |> invisible()
    out <- if (sceptre_object@run_permutations) {
        run_permutation_analysis(context, execution, synthetic_idxs)
    } else {
        run_crt_analysis(context, execution)
    }
    out$result <- assemble_pair_results(out$result, response_grna_group_pairs)
    out
}
