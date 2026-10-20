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
#'   `get_result_methods()` returns a named character vector of current
#'   completed families and method identifiers; absent or stale families
#'   are omitted.
#' @details Uses the configured sidedness, QC, guide integration, and multiple
#'   testing policy. Union pools assigned cells; singleton tests each guide;
#'   Bonferroni uses the existing minimum-p guide reduction before correction
#'   over the family. Failed-QC rows remain missing and are not discoveries.
#'   An all-tied eligible comparison returns p = 1. Empty groups are errors.
#'
#'   `fold_change` and `log_2_fold_change` describe the ratio of arithmetic mean
#'   CP10K responses, not a Wilcoxon location estimate or adjusted effect.
#'   `se_fold_change` is unavailable. No covariates are fitted, regardless of
#'   the
#'   object's formula. Library-size scaling does not remove batch, cell-state,
#'   donor, clonality, or assignment confounding. Cells must be independent and
#'   exchangeable under the null; this is not donor-level inference. Asymptotic
#'   p-values can be inaccurate for small or highly discrete samples. The usual
#'   assumptions of the chosen multiple-testing procedure still apply. A null
#'   control check is a diagnostic, not proof of FDR control.
#'
#'   P-value arithmetic delegates to the installed version of R's stats package.
#'   Numerical cancellation or underflow can return zero even when a tail
#'   probability remains representable. Raw zeros are retained. Log-scale
#'   plots apply a display floor only; they do not repair the probabilities.
#'   Extreme-tail magnitudes and extreme-tail calibration remain unqualified.
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
#'     set_analysis_parameters(
#'         discovery_pairs = pairs, control_group = "nt_cells"
#'     ) |>
#'     assign_grnas(method = "thresholding") |>
#'     run_qc()
#' alternative <- run_wilcoxon(object)
#' get_result(alternative, "run_discovery_analysis")
#' get_result_methods(alternative)

run_wilcoxon <- function(
    sceptre_object,
    analysis = "discovery_analysis",
    execution = "direct",
    parallel = FALSE,
    n_calibration_pairs = NULL,
    calibration_group_size = NULL
) {
    object <- validate_wilcoxon_support(sceptre_object, execution, parallel)
    analysis <- analysis_record_name(analysis)
    plan <- prepare_wilcoxon_plan(
        object,
        analysis,
        n_calibration_pairs,
        calibration_group_size
    )
    raw <- test_wilcoxon_family(object, plan)
    result <- finalize_wilcoxon_family(object, plan, raw)
    commit_wilcoxon_record(object, plan, result)
}

# Mode admission precedes storage access, pair sampling and execution.
validate_wilcoxon_support <- function(sceptre_object, execution, parallel) {
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
    object
}

# Reuse pair construction and its QC checks, with no unused resampling budget.
prepare_wilcoxon_plan <- function(object, analysis, n_pairs, group_size) {
    budgets <- c(B1 = 0L, B2 = 0L, B3 = 0L)
    if (identical(analysis, "run_calibration_check")) {
        prepared <- run_calibration_check_pt_1(
            object,
            n_calibration_pairs = n_pairs,
            calibration_group_size = group_size,
            budgets = budgets
        )
        pairs <- prepared$plan@pairs
        n_pairs <- prepared$plan@n_calibration_pairs
        group_size <- prepared$plan@calibration_group_size
    } else {
        if (!is.null(n_pairs) || !is.null(group_size)) {
            stop("Calibration sizes apply only to calibration_check.")
        }
        pairs <- analysis_pairs(object, analysis)
    }
    if (!nrow(pairs)) {
        stop("Wilcoxon requires a nonempty prepared family.")
    }
    make_run_plan(
        object,
        analysis,
        pairs,
        n_calibration_pairs = n_pairs,
        calibration_group_size = group_size,
        budgets = budgets
    )
}

prepare_wilcoxon_context <- function(object, plan) {
    eligible <- plan@pairs[plan@pairs$pass_qc, , drop = FALSE]
    raw <- eligible
    raw$p_value <- raw$fold_change <- raw$se_fold_change <-
        raw$log_2_fold_change <- rep(NA_real_, nrow(raw))
    index <- prepare_index_context(
        object@preparation@analysis_cells,
        ncol(get_response_matrix(object)),
        object@preparation@analysis_groups,
        eligible,
        rownames(get_response_matrix(object)),
        identical(plan@analysis, "run_calibration_check"),
        FALSE,
        TRUE
    )
    list(
        index = index,
        raw = raw,
        cells = object@preparation@analysis_cells,
        depth = object@data@covariate_data_frame$response_n_umis,
        alternative = c("less", "two.sided", "greater")[
            object@spec@side_code + 2L
        ]
    )
}

# Reuse worker error context while admitting only a single serial block.
test_wilcoxon_family <- function(object, plan) {
    context <- prepare_wilcoxon_context(object, plan)
    worker <- function(rows, proc_id = NULL) {
        test_wilcoxon_response_rows(object, context, rows)
    }
    execute_analysis_blocks(
        list(seq_len(nrow(context$raw))),
        worker,
        parallel = FALSE,
        analysis = plan@analysis,
        phase = "wilcoxon_tests"
    )[[1L]]
}

# Read each response once; guide views select from normalized QC-cell counts.
test_wilcoxon_response_rows <- function(object, context, rows) {
    raw <- context$raw
    index <- context$index
    for (response_idx in unique(index$pair_rows$response_idx[rows])) {
        y <- storage_read_response(
            get_response_matrix(object),
            response_idx,
            context$cells
        )
        normalized <- 10000 * y / context$depth[context$cells]
        request_rows <- rows[index$pair_rows$response_idx[rows] == response_idx]
        for (i in request_rows) {
            view <- make_fit_view(index, raw$grna_group[i])
            response <- normalized[view@fit_to_analysis]
            trt <- view@trt_fit_pos
            control <- setdiff(seq_along(response), trt)
            if (!length(trt) || !length(control)) {
                stop("Wilcoxon requires nonempty treatment and control groups.")
            }
            raw$p_value[i] <- wilcoxon_probability(
                log1p(response[trt]),
                log1p(response[control]),
                context$alternative
            )
            raw$fold_change[i] <- mean(response[trt]) / mean(response[control])
            raw$log_2_fold_change[i] <- log2(raw$fold_change[i])
        }
    }
    raw
}

finalize_wilcoxon_family <- function(object, plan, raw) {
    result <- assemble_pair_results(raw, plan@pairs)
    result <- switch(
        plan@analysis,
        run_power_check = apply_grouping_to_result(result, object),
        run_calibration_check = process_calibration_result(result, object),
        run_discovery_analysis = process_discovery_result(result, object)
    )
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
    result
}

commit_wilcoxon_record <- function(object, plan, result) {
    analysis <- plan@analysis
    object <- invalidate_other_inference_results(object, "wilcoxon_unadjusted")
    object <- commit_result_record(
        object,
        analysis,
        result,
        plan@pairs,
        1L,
        plan
    )
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
# preserves historical results and provenance without relabelling them.
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

alternative_inference_text <- function(object) {
    lines <- lapply(analysis_names(), function(analysis) {
        if (!result_record_is_current(object, analysis)) {
            return("")
        }
        method <- object@results@records[[analysis]]$provenance$inference_method
        if (!identical(method, "wilcoxon_unadjusted")) {
            return("")
        }
        paste0(
            "\n",
            analysis,
            ": experimental unadjusted Wilcoxon; ",
            "log1p(CP10K), full-library depth; no covariate adjustment.\n",
            "Effect: descriptive mean CP10K ratio; no effect SE.\n",
            "Configured formula/resampling settings are ignored.\n"
        )
    })
    paste0(unlist(lines, use.names = FALSE), collapse = "")
}

#' @rdname run_wilcoxon
#' @export
get_result_methods <- function(sceptre_object) {
    object <- migrate_sceptre_object(sceptre_object)
    current <- analysis_names()[vapply(
        analysis_names(),
        function(analysis) {
            result_record_is_current(object, analysis)
        },
        logical(1)
    )]
    stats::setNames(
        vapply(
            current,
            function(analysis) {
                method <- object@results@records[[
                    analysis
                ]]$provenance$inference_method
                if (is.null(method)) {
                    return("sceptre")
                }
                if (
                    !is.character(method) ||
                        length(method) != 1L ||
                        is.na(method) ||
                        !method %in% c("sceptre", "wilcoxon_unadjusted")
                ) {
                    stop("Invalid inference method in a current result record.")
                }
                method
            },
            character(1)
        ),
        current
    )
}

has_current_wilcoxon_result <- function(object) {
    any(get_result_methods(object) == "wilcoxon_unadjusted")
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
                    "Effect: descriptive mean CP10K ratio;",
                    "no covariate adjustment. Plot p floor: 1e-300.",
                    sep = "\n"
                )
            ) +
            ggplot2::theme(
                plot.caption = ggplot2::element_text(
                    size = 8,
                    hjust = 0,
                    colour = "black"
                ),
                plot.background = ggplot2::element_rect(
                    fill = "white",
                    colour = NA
                ),
                plot.caption.position = "plot"
            )
    }
    if (inherits(plot, "ggplot")) {
        return(label(plot))
    }
    lapply(plot, label)
}

require_supported_distributed_inference <- function(object) {
    if (has_current_wilcoxon_result(object)) {
        stop(
            "Experimental Wilcoxon is unsupported in distributed execution. ",
            "Use run_wilcoxon() with direct execution or a separate ",
            "prepared SCEPTRE object."
        )
    }
    invisible(NULL)
}

# Numerical cancellation or underflow can yield zero. Preserve raw results;
# alternative-only log-scale plots disclose their display floor.
alternative_plot_probabilities <- function(result) {
    if (
        !"inference_method" %in% names(result) ||
            !identical(unique(result$inference_method), "wilcoxon_unadjusted")
    ) {
        return(result)
    }
    result <- owned_result_value(result)
    result$p_value <- pmax(result$p_value, 1e-300)
    result
}

alternative_volcano_threshold <- function(threshold, result) {
    if (
        "inference_method" %in%
            names(result) &&
            identical(unique(result$inference_method), "wilcoxon_unadjusted")
    ) {
        return(pmax(threshold, 1e-12))
    }
    threshold
}
