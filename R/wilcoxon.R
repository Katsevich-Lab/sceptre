#' Run an experimental unadjusted Wilcoxon analysis
#'
#' Explicit opt-in alternative to covariate-adjusted SCEPTRE inference. Requires
#' completed assignment and QC, in-memory responses, low MOI, and
#' non-targeting controls. Tests
#' log1p(10000 * counts / full-library UMIs) with the asymptotic, tie-corrected
#' two-sample rank-sum test in `stats::wilcox.test`, with continuity correction.
#' The null is equality of the two cell-level response distributions. The rank
#' estimand is P(treatment > control) + 0.5 P(treatment = control); a location
#' shift interpretation requires additional distributional assumptions.
#'
#' @param sceptre_object A prepared `sceptre_object` after `run_qc()`.
#' @param analysis One of `"discovery_analysis"`, `"power_check"`, or
#'   `"calibration_check"` (the corresponding `run_` names are also accepted).
#' @param execution Must be `"direct"`. Distributed/Nextflow is unsupported.
#' @param parallel Must be `FALSE`. Parallel execution is unsupported.
#' @inheritParams run_calibration_check
#' @return An updated object with a complete-family result available through
#'   `get_result()`. This replaces the selected analysis record and invalidates
#'   current records from other inference methods. Keep a separate object for
#'   SCEPTRE comparisons. Fitted-model caches are not used or modified.
#' @details Uses the configured sidedness, QC, guide integration, and multiple
#'   testing policy. Union pools assigned cells; singleton tests each guide;
#'   Bonferroni uses the existing minimum-p guide reduction before correction
#'   over the family. Failed-QC rows remain missing and are not discoveries.
#'   An all-tied eligible comparison returns p = 1. Empty groups are errors.
#'
#'   `fold_change` and `log_2_fold_change` describe the ratio of arithmetic mean
#'   CP10K responses, not a Wilcoxon location estimate or adjusted effect.
#'   `se_fold_change` is unavailable. No covariates are fitted, regardless of the
#'   object's formula. Library-size scaling does not remove batch, cell-state,
#'   donor, clonality, or assignment confounding. Cells must be independent and
#'   exchangeable under the null; this is not donor-level inference. Asymptotic
#'   p-values can be inaccurate for small or highly discrete samples. The usual
#'   assumptions of the chosen multiple-testing procedure still apply. A null
#'   control check is a diagnostic, not proof of FDR control.
#'
#'   Resampling settings in the object are ignored: no permutations, CRT draws,
#'   regression, or tail fitting are performed. The saved run plan uses a zero
#'   resampling budget. Running a standard analysis function later selects
#'   SCEPTRE again and invalidates current Wilcoxon records.
#' @export
#' @examples
#' data(lowmoi_example_data)
#' d <- lowmoi_example_data
#' object <- import_data(
#'     response_matrix = d$response_matrix,
#'     grna_matrix = d$grna_matrix,
#'     grna_target_data_frame = d$grna_target_data_frame,
#'     moi = "low", extra_covariates = d$extra_covariates
#' )
#' # One arbitrary pair illustrates the API, not a known biological effect.
#' target <- setdiff(d$grna_target_data_frame$grna_target, "non-targeting")[1L]
#' pairs <- data.frame(
#'     grna_target = target, response_id = rownames(d$response_matrix)[1L]
#' )
#' object <- object |>
#'     set_analysis_parameters(discovery_pairs = pairs, control_group = "nt_cells") |>
#'     assign_grnas(method = "thresholding") |>
#'     run_qc()
#' alternative <- run_wilcoxon(object)
#' get_result(alternative, "run_discovery_analysis")

run_wilcoxon <- function(
    sceptre_object,
    analysis = "discovery_analysis",
    execution = "direct",
    parallel = FALSE,
    n_calibration_pairs = NULL,
    calibration_group_size = NULL
) {
    if (!identical(execution, "direct") || !identical(parallel, FALSE)) {
        stop("Experimental Wilcoxon supports direct serial execution only.")
    }
    object <- migrate_sceptre_object(sceptre_object)
    require_stage(object, "run_qc")
    if (
        !isTRUE(object@data@low_moi) ||
            !identical(object@spec@control_group_complement, FALSE)
    ) {
        stop("Experimental Wilcoxon requires low MOI and NT-cell controls.")
    }
    if (methods::is(object@data@response_matrix, "BackedCountStore")) {
        stop("Experimental Wilcoxon currently requires in-memory responses.")
    }
    analysis <- analysis_record_name(analysis)
    calibration <- identical(analysis, "run_calibration_check")
    if (
        !calibration &&
            (!is.null(n_calibration_pairs) || !is.null(calibration_group_size))
    ) {
        stop("Calibration sizes apply only to calibration_check.")
    }
    depth <- object@data@covariate_data_frame$response_n_umis
    cells <- object@preparation@analysis_cells
    if (
        length(depth) != ncol(get_response_matrix(object)) ||
            any(!is.finite(depth[cells]) | depth[cells] <= 0)
    ) {
        stop(
            "Wilcoxon requires positive finite full-library UMIs for QC cells."
        )
    }
    if (calibration) {
        prepared <- run_calibration_check_pt_1(
            object,
            n_calibration_pairs = n_calibration_pairs,
            calibration_group_size = calibration_group_size
        )
        pairs <- prepared$plan@pairs
        n_calibration_pairs <- prepared$plan@n_calibration_pairs
        calibration_group_size <- prepared$plan@calibration_group_size
    } else {
        pairs <- analysis_pairs(object, analysis)
    }
    if (!nrow(pairs)) {
        stop("Wilcoxon requires a nonempty prepared family.")
    }
    plan <- make_run_plan(
        object,
        analysis,
        pairs,
        n_calibration_pairs = n_calibration_pairs,
        calibration_group_size = calibration_group_size,
        budgets = c(B1 = 0L, B2 = 0L, B3 = 0L)
    )
    eligible <- pairs[pairs$pass_qc, , drop = FALSE]
    index <- prepare_index_context(
        cells,
        ncol(get_response_matrix(object)),
        object@preparation@analysis_groups,
        eligible,
        rownames(get_response_matrix(object)),
        calibration,
        FALSE,
        TRUE
    )
    raw <- eligible
    raw$p_value <- raw$fold_change <- raw$se_fold_change <-
        raw$log_2_fold_change <- rep(NA_real_, nrow(raw))
    alternative <- c("less", "two.sided", "greater")[
        object@spec@side_code + 2L
    ]
    analyze_block <- function(rows, proc_id = NULL) {
        for (response_idx in unique(index$pair_rows$response_idx[rows])) {
            y <- storage_read_response(
                get_response_matrix(object),
                response_idx,
                cells
            )
            normalized <- 10000 * y / depth[cells]
            request_rows <- rows[
                index$pair_rows$response_idx[rows] == response_idx
            ]
            for (i in request_rows) {
                view <- make_fit_view(index, raw$grna_group[i])
                response <- normalized[view@fit_to_analysis]
                trt <- view@trt_fit_pos
                control <- setdiff(seq_along(response), trt)
                if (!length(trt) || !length(control)) {
                    stop(
                        "Wilcoxon requires nonempty treatment and control groups."
                    )
                }
                raw$p_value[i] <- wilcoxon_probability(
                    log1p(response[trt]),
                    log1p(response[control]),
                    alternative
                )
                raw$fold_change[i] <- mean(response[trt]) /
                    mean(response[control])
                raw$log_2_fold_change[i] <- log2(raw$fold_change[i])
            }
        }
        raw
    }
    raw <- execute_analysis_blocks(
        list(seq_len(nrow(raw))),
        analyze_block,
        parallel = FALSE,
        analysis = analysis,
        phase = "wilcoxon_tests"
    )[[1L]]

    result <- assemble_pair_results(raw, pairs)
    result <- if (identical(analysis, "run_power_check")) {
        apply_grouping_to_result(result, object)
    } else if (calibration) {
        process_calibration_result(result, object)
    } else {
        process_discovery_result(result, object)
    }
    if (!"fold_change" %in% names(result)) {
        result$fold_change <- 2^result$log_2_fold_change
    }
    if (!"se_fold_change" %in% names(result)) {
        result$se_fold_change <- rep(NA_real_, nrow(result))
    }
    result$inference_method <- rep("wilcoxon_unadjusted", nrow(result))
    result$response_scale <- rep("log1p_CP10K_full_library", nrow(result))
    result$effect_estimator <- rep("mean_CP10K_ratio", nrow(result))
    result$test_status <- ifelse(result$pass_qc, "tested", "failed_qc")
    object <- invalidate_other_inference_results(object, "wilcoxon_unadjusted")
    object <- commit_result_record(object, analysis, result, pairs, 1L, plan)
    provenance <- object@results@records[[analysis]]$provenance
    provenance$inference_method <- "wilcoxon_unadjusted"
    provenance$response_fit_method <- provenance$perturbation_fit_method <-
        "not_used"
    provenance$normalization <- "log1p(10000 * count / full-library UMIs)"
    provenance$test <- "stats::wilcox.test; exact=FALSE; correct=TRUE"
    provenance$covariate_adjustment <- FALSE
    provenance$execution <- "direct_serial"
    object@results@records[[analysis]]$provenance <- provenance
    object
}

# Handle the degenerate variance explicitly; base R returns NaN for all ties.
wilcoxon_probability <- function(treatment, control, alternative) {
    if (length(unique(c(treatment, control))) == 1L) {
        return(1)
    }
    stats::wilcox.test(
        treatment,
        control,
        alternative = alternative,
        exact = FALSE,
        correct = TRUE
    )$p.value
}

# Absence means the established SCEPTRE method for older records. Invalidation
# preserves historical results and their provenance rather than relabelling them.
invalidate_other_inference_results <- function(object, method = "sceptre") {
    for (analysis in analysis_names()) {
        record <- object@results@records[[analysis]]
        if (is.null(record)) {
            next
        }
        saved_method <- record$provenance$inference_method
        if (is.null(saved_method)) {
            saved_method <- "sceptre"
        }
        if (!identical(saved_method, method)) {
            object@results@records[[analysis]]$invalidated_by <-
                paste("Inference method changed to", method)
        }
    }
    object
}

print_alternative_inference <- function(object) {
    for (analysis in analysis_names()) {
        if (!result_record_is_current(object, analysis)) {
            next
        }
        method <- object@results@records[[analysis]]$provenance$inference_method
        if (identical(method, "wilcoxon_unadjusted")) {
            cat(
                "\n",
                analysis,
                ": experimental unadjusted Wilcoxon; ",
                "log1p(CP10K), full-library depth; no covariate adjustment.\n",
                "Effect: descriptive mean CP10K ratio; no effect SE.\n",
                "Configured formula/resampling settings are ignored.\n",
                sep = ""
            )
        }
    }
}

has_current_wilcoxon_result <- function(object) {
    any(vapply(
        analysis_names(),
        function(analysis) {
            result_record_is_current(object, analysis) &&
                identical(
                    object@results@records[[
                        analysis
                    ]]$provenance$inference_method,
                    "wilcoxon_unadjusted"
                )
        },
        logical(1)
    ))
}

label_alternative_plot <- function(plot, object, analysis) {
    method <- object@results@records[[analysis]]$provenance$inference_method
    if (!identical(method, "wilcoxon_unadjusted")) {
        return(plot)
    }
    label <- function(p) {
        p +
            ggplot2::labs(
                caption = paste(
                    "Experimental unadjusted Wilcoxon: log1p(CP10K).",
                    "Effect: descriptive mean CP10K ratio; no covariate adjustment."
                )
            )
    }
    if (inherits(plot, "ggplot")) {
        return(label(plot))
    }
    lapply(plot, label)
}
