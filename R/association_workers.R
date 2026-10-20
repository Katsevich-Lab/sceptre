#' Test several guide groups using the same response fit and permutation draws
#' @param expression Gene's counts in the selected model-fitting cell order.
#' @param score Fitted means and score arrays in that same cell order.
#' @param treatment_lookup Function returning each group's treated-cell
#'   positions
#'   in those arrays, starting at 1, and their count.
#' @param groups Guide-group names in request order, including repeats.
#' @param synthetic_idxs Complete shared set of permutation draws. Draw
#'   positions
#'   start at 0 for C++; observed treatment positions start at 1 and are
#'   converted
#'   by `run_association_test()`.
#' @param test_options Stage sizes B1/B2/B3, test direction and approximation
#'   setting.
#' @param output_amount Public output level; level 3 includes resampled
#'   statistics.
#' @return C++ result lists in guide-group order. The response fit, score arrays
#'   and random draws are already prepared; this worker applies the tests.
#' @noRd
test_permutation_shared_fit <- function(
    expression,
    score,
    treatment_lookup,
    groups,
    synthetic_idxs,
    test_options,
    output_amount,
    response_id = NULL
) {
    results <- vector("list", length(groups))
    for (i in seq_along(groups)) {
        treatment <- treatment_lookup(groups[i])
        results[[i]] <- test_association_pair(
            y = expression,
            mu = score$mu,
            a = score$a,
            w = score$w,
            D = score$D,
            trt_idxs = treatment$trt_idxs,
            use_all_cells = FALSE,
            n_trt = treatment$n_trt,
            synthetic_idxs = synthetic_idxs,
            B1 = test_options$B1,
            B2 = test_options$B2,
            B3 = test_options$B3,
            fit_parametric_curve = test_options$fit_parametric_curve,
            return_resampling_dist = (output_amount == 3L),
            side_code = test_options$side_code,
            response_id = response_id,
            grna_group = groups[i]
        )
    }
    results
}

#' Fit and test one gene separately for each group's treatment-plus-NT cells
#' @param expression Gene's counts in the order of cells kept after QC.
#' @param design Design matrix with rows in that same order, even for an
#'   intercept-only formula.
#' @param index Cell lookups prepared for discovery with NT controls.
#' @param groups Guide-group names in request order.
#' @param synthetic_idxs Complete shared hybrid-permutation draws. Each test
#'   uses
#'   the first n_trt entries of each draw for its own treatment-cell count.
#' @param test_options,output_amount Validated test settings and output level.
#' @return C++ results in guide-group order. Each group gets its own response
#'   fit,
#'   with treatment cells first, then NT cells, keeping membership order and
#'   repeats.
#'   These group-specific fits are not saved in the shared response-fit list.
#'   The supplied draws are reused, rather than generated again here.
#' @noRd
test_permutation_group_fit <- function(
    expression,
    design,
    index,
    groups,
    synthetic_idxs,
    test_options,
    output_amount,
    response_id = NULL
) {
    results <- vector("list", length(groups))
    for (i in seq_along(groups)) {
        view <- make_fit_view(index, groups[i])
        group_expression <- expression[view$fit_to_analysis]
        group_design <- design[view$fit_to_analysis, , drop = FALSE]
        fit <- fit_response(group_expression, group_design)
        score <- prepare_score_components(group_expression, group_design, fit)
        results[[i]] <- test_association_pair(
            y = group_expression,
            mu = score$mu,
            a = score$a,
            w = score$w,
            D = score$D,
            trt_idxs = view$trt_fit_pos,
            use_all_cells = FALSE,
            n_trt = view$n_trt,
            synthetic_idxs = synthetic_idxs,
            B1 = test_options$B1,
            B2 = test_options$B2,
            B3 = test_options$B3,
            fit_parametric_curve = test_options$fit_parametric_curve,
            return_resampling_dist = (output_amount == 3L),
            side_code = test_options$side_code,
            response_id = response_id,
            grna_group = groups[i]
        )
    }
    results
}

#' Test one guide group's requested genes using shared response fits
#' @param response_store Response count matrix, readable one gene at a time.
#' @param response_rows Original matrix gene-row numbers, starting at 1, in
#'   request order.
#' @param response_ids Corresponding gene names used to find their saved fits.
#' @param design Design matrix whose rows follow `source_cells`.
#' @param source_cells Original matrix column numbers for this fit, in order;
#'   repeats are allowed.
#' @param treatment Treated cells' positions within this fit, starting at 1, and
#'   their count.
#' @param response_fits Previously prepared response fits for these same
#'   cells/design.
#' @param test_options,output_amount Validated test settings and output level.
#' @param probabilities Optional probabilities reconstructed from a compatible
#'   compact perturbation fit. NULL follows the original one-fit cold path.
#' @return C++ results in requested gene order. Fit assignment probabilities and
#'   generate B1+B2+B3 CRT draws once for this guide, then read each gene and
#'   use
#'   its saved response fit. This worker does not add new response fits.
#' @noRd
test_crt_shared_fit <- function(
    response_store,
    response_rows,
    response_ids,
    design,
    source_cells,
    treatment,
    response_fits,
    test_options,
    output_amount,
    probabilities = NULL,
    grna_group = NULL
) {
    if (is.null(probabilities)) {
        probabilities <- fit_assignment_probabilities(
            treatment$trt_idxs,
            design
        )
    }
    synthetic_idxs <- sample_crt_indices(
        probabilities,
        B = test_options$total_resamples
    )
    results <- vector("list", length(response_rows))
    for (i in seq_along(response_rows)) {
        expression <- storage_read_response(
            response_store,
            response_rows[i],
            source_cells
        )
        score <- prepare_score_components(
            expression,
            design,
            lookup_response_fit(response_fits, response_ids[i])
        )
        results[[i]] <- test_association_pair(
            y = expression,
            mu = score$mu,
            a = score$a,
            w = score$w,
            D = score$D,
            trt_idxs = treatment$trt_idxs,
            use_all_cells = TRUE,
            n_trt = treatment$n_trt,
            synthetic_idxs = synthetic_idxs,
            B1 = test_options$B1,
            B2 = test_options$B2,
            B3 = test_options$B3,
            fit_parametric_curve = test_options$fit_parametric_curve,
            return_resampling_dist = (output_amount == 3L),
            side_code = test_options$side_code,
            response_id = response_ids[i],
            grna_group = grna_group
        )
    }
    results
}

#' Fit and test one guide group's genes using its treatment-plus-NT cells
#' @param response_store Response count matrix, readable one gene at a time.
#' @param response_rows Original matrix gene-row numbers, starting at 1, in
#'   request order.
#' @param design Design matrix for all cells kept after QC, in that order.
#' @param fit_view Cell selection from `make_fit_view()`: `fit_to_analysis`
#'   selects
#'   design rows; `source_cells` selects original count-matrix columns;
#'   `trt_fit_pos` locates treatment cells within the selected data, starting at
#'   1.
#' @param test_options,output_amount Validated test settings and output level.
#' @param probabilities Optional probabilities reconstructed from a compatible
#'   compact perturbation fit. NULL follows the original one-fit cold path.
#' @param group_design Optional existing slice of design in fit_view order.
#'   NULL constructs the same slice locally.
#' @return C++ results in requested gene order. Fit assignment probabilities and
#'   generate a full set of CRT draws once for this guide. Then read, fit and
#'   test
#'   each gene using the selected cells. These response fits belong only to this
#'   group's cells and are not saved for other groups to reuse.
#' @noRd
test_crt_group_fit <- function(
    response_store,
    response_rows,
    design,
    fit_view,
    test_options,
    output_amount,
    probabilities = NULL,
    group_design = NULL,
    grna_group = NULL
) {
    if (is.null(group_design)) {
        group_design <- design[fit_view$fit_to_analysis, , drop = FALSE]
    }
    if (is.null(probabilities)) {
        probabilities <- fit_assignment_probabilities(
            fit_view$trt_fit_pos,
            group_design
        )
    }
    synthetic_idxs <- sample_crt_indices(
        probabilities,
        B = test_options$total_resamples
    )
    results <- vector("list", length(response_rows))
    for (i in seq_along(response_rows)) {
        expression <- storage_read_response(
            response_store,
            response_rows[i],
            fit_view$source_cells
        )
        fit <- fit_response(expression, group_design)
        score <- prepare_score_components(expression, group_design, fit)
        results[[i]] <- test_association_pair(
            y = expression,
            mu = score$mu,
            a = score$a,
            w = score$w,
            D = score$D,
            trt_idxs = fit_view$trt_fit_pos,
            use_all_cells = TRUE,
            n_trt = fit_view$n_trt,
            synthetic_idxs = synthetic_idxs,
            B1 = test_options$B1,
            B2 = test_options$B2,
            B3 = test_options$B3,
            fit_parametric_curve = test_options$fit_parametric_curve,
            return_resampling_dist = (output_amount == 3L),
            side_code = test_options$side_code,
            response_id = rownames(response_store)[response_rows[i]],
            grna_group = grna_group
        )
    }
    results
}

# Add pair identity after native cleanup, while retaining the original cause.
# Errors abort the complete family. No pair is silently filtered or
# assigned p=1.
test_association_pair <- function(..., response_id = NULL, grna_group = NULL) {
    tryCatch(run_association_test(...), error = function(cause) {
        response <- if (is.null(response_id)) "<unknown>" else response_id
        group <- if (is.null(grna_group)) "<unknown>" else grna_group
        stop(structure(list(
            message = paste0("Association test for response '", response,
                "' and guide group '", group, "' failed: ",
                conditionMessage(cause)),
            call = NULL, response_id = response_id, grna_group = grna_group,
            cause = cause, parent = cause),
            class = c("sceptre_pair_error", "error", "condition")))
    })
}
