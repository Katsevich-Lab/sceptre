#' Test guide groups against one shared permutation response fit.
#'
#' @param expression Numeric response in the common fit-view cell order.
#' @param score Full score components aligned to expression and the common design.
#' @param treatment_lookup Prepared lookup returning one-based fit-local treatment
#'   positions and their count for a guide-group ID.
#' @param groups Character guide-group IDs in request order, including repeats.
#' @param synthetic_idxs One complete shared permutation bank; native draws are
#'   zero-based, whereas observed positions enter native inference one-based.
#' @param test_options Checked stage budgets, side code and approximation flag.
#' @param output_amount Existing public output level (3 retains resampling values).
#' @return Native result lists in group order. Does not read storage, fit models,
#'   draw resamples, mutate inputs, or update the response cache.
#' @noRd
test_permutation_shared_fit <- function(expression, score, treatment_lookup,
                                        groups, synthetic_idxs, test_options,
                                        output_amount) {
    results <- vector("list", length(groups))
    for (i in seq_along(groups)) {
        treatment <- treatment_lookup(groups[i])
        results[[i]] <- run_low_level_test_full_v4(
            y = expression, mu = score$mu, a = score$a, w = score$w, D = score$D,
            trt_idxs = treatment$trt_idxs, use_all_cells = FALSE,
            n_trt = treatment$n_trt, synthetic_idxs = synthetic_idxs,
            B1 = test_options$B1, B2 = test_options$B2, B3 = test_options$B3,
            fit_parametric_curve = test_options$fit_parametric_curve,
            return_resampling_dist = (output_amount == 3L),
            side_code = test_options$side_code
        )
    }
    results
}

#' Fit and test each treatment-plus-NT permutation view of one response.
#'
#' @param expression Numeric response in analysis-cell order.
#' @param design Matrix in that same common analysis-cell order, including when
#'   it has only an intercept column.
#' @param index Prepared coordinate context for NT-only discovery.
#' @param groups Character requested groups in their existing order.
#' @param synthetic_idxs Complete shared hybrid bank; inference consumes the
#'   existing treatment-count prefix for each fit view.
#' @param test_options,output_amount Checked test settings and output level.
#' @return Native results in group order. Fits one response per group and keeps
#'   those guide-specific fits out of the shared response cache. No input mutation
#'   or RNG draws; each view retains treatment-first, then NT order and repeats.
#' @noRd
test_permutation_group_fit <- function(expression, design, index, groups,
                                       synthetic_idxs, test_options,
                                       output_amount) {
    results <- vector("list", length(groups))
    for (i in seq_along(groups)) {
        view <- make_fit_view(index, groups[i])
        group_expression <- expression[view$fit_to_analysis]
        group_design <- design[view$fit_to_analysis, , drop = FALSE]
        fit <- fit_response(group_expression, group_design)
        score <- prepare_score_components(group_expression, group_design, fit)
        results[[i]] <- run_low_level_test_full_v4(
            y = group_expression, mu = score$mu, a = score$a, w = score$w,
            D = score$D, trt_idxs = view$trt_fit_pos, use_all_cells = FALSE,
            n_trt = view$n_trt, synthetic_idxs = synthetic_idxs,
            B1 = test_options$B1, B2 = test_options$B2, B3 = test_options$B3,
            fit_parametric_curve = test_options$fit_parametric_curve,
            return_resampling_dist = (output_amount == 3L),
            side_code = test_options$side_code
        )
    }
    results
}

#' Test requested responses with a shared CRT assignment fit for one guide.
#'
#' @param response_store Row-readable response matrix or storage adapter.
#' @param response_rows One-based source response rows in request order.
#' @param response_ids Corresponding cache keys in that same order.
#' @param design Common fit-view matrix; rows align to source_cells.
#' @param source_cells One-based source-cell positions, ordered and repeatable.
#' @param treatment One-based fit-local observed positions and their count.
#' @param response_fits Complete compatible compact response-fit cache.
#' @param test_options,output_amount Checked test settings and output level.
#' @return Native results in response order. Fits the assignment model and draws
#'   one full B1+B2+B3 bank before reading responses. Reads each requested row,
#'   reuses cached response fits, and never inserts fits or mutates inputs.
#' @noRd
test_crt_shared_fit <- function(response_store, response_rows, response_ids,
                                design, source_cells, treatment, response_fits,
                                test_options, output_amount) {
    probabilities <- fit_assignment_probabilities(treatment$trt_idxs, design)
    synthetic_idxs <- crt_index_sampler_fast(probabilities,
                                             B = test_options$total_resamples)
    results <- vector("list", length(response_rows))
    for (i in seq_along(response_rows)) {
        expression <- storage_read_response(response_store, response_rows[i], source_cells)
        score <- prepare_score_components(expression, design,
            lookup_response_fit(response_fits, response_ids[i]))
        results[[i]] <- run_low_level_test_full_v4(
            y = expression, mu = score$mu, a = score$a, w = score$w, D = score$D,
            trt_idxs = treatment$trt_idxs, use_all_cells = TRUE,
            n_trt = treatment$n_trt, synthetic_idxs = synthetic_idxs,
            B1 = test_options$B1, B2 = test_options$B2, B3 = test_options$B3,
            fit_parametric_curve = test_options$fit_parametric_curve,
            return_resampling_dist = (output_amount == 3L),
            side_code = test_options$side_code
        )
    }
    results
}

#' Fit and test requested responses in one treatment-plus-NT CRT view.
#'
#' @param response_store Row-readable response storage.
#' @param response_rows One-based source response rows in request order.
#' @param design Common matrix in analysis-cell order.
#' @param fit_view Prepared guide view: fit_to_analysis selects design rows,
#'   source_cells selects storage cells, and trt_fit_pos is one-based fit-local.
#' @param test_options,output_amount Checked test settings and output level.
#' @return Native results in response order. Fits one assignment and draws one
#'   full bank per guide, then reads and fits every response in the ordered view.
#'   Guide-specific response fits are not cached. Does not mutate caller inputs.
#' @noRd
test_crt_group_fit <- function(response_store, response_rows, design, fit_view,
                               test_options, output_amount) {
    group_design <- design[fit_view$fit_to_analysis, , drop = FALSE]
    probabilities <- fit_assignment_probabilities(fit_view$trt_fit_pos, group_design)
    synthetic_idxs <- crt_index_sampler_fast(probabilities,
                                             B = test_options$total_resamples)
    results <- vector("list", length(response_rows))
    for (i in seq_along(response_rows)) {
        expression <- storage_read_response(response_store, response_rows[i],
                                            fit_view$source_cells)
        fit <- fit_response(expression, group_design)
        score <- prepare_score_components(expression, group_design, fit)
        results[[i]] <- run_low_level_test_full_v4(
            y = expression, mu = score$mu, a = score$a, w = score$w, D = score$D,
            trt_idxs = fit_view$trt_fit_pos, use_all_cells = TRUE,
            n_trt = fit_view$n_trt, synthetic_idxs = synthetic_idxs,
            B1 = test_options$B1, B2 = test_options$B2, B3 = test_options$B3,
            fit_parametric_curve = test_options$fit_parametric_curve,
            return_resampling_dist = (output_amount == 3L),
            side_code = test_options$side_code
        )
    }
    results
}
