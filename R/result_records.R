# Resolved scientific choices live in AnalysisSpec. Execution settings and the
# realized family belong to a transient plan and are frozen only with a result.
result_analysis_fields <- function() {
    stats::setNames(analysis_names(), analysis_names())
}

analysis_record_name <- function(analysis) {
    aliases <- c(
        calibration_check = "run_calibration_check",
        power_check = "run_power_check",
        discovery_analysis = "run_discovery_analysis"
    )
    if (!is.character(analysis) || length(analysis) != 1L || is.na(analysis)) {
        stop("Unknown analysis type.")
    }
    if (analysis %in% names(aliases)) {
        analysis <- aliases[[analysis]]
    }
    if (!analysis %in% analysis_names()) {
        stop("Unknown analysis type.")
    }
    analysis
}

analysis_pairs <- function(object, analysis) {
    switch(
        analysis_record_name(analysis),
        run_calibration_check = stop(
            "Calibration pairs must be supplied by a run plan."
        ),
        run_power_check = object@preparation@positive_control_pairs_with_info,
        run_discovery_analysis = object@preparation@discovery_pairs_with_info
    )
}

analysis_spec_snapshot <- function(object) {
    names <- methods::slotNames(object@spec)
    identity_value(stats::setNames(
        lapply(names, function(name) methods::slot(object@spec, name)),
        names
    ))
}

run_plan_identity <- function(object, plan) {
    options <- c(
        "side_code",
        "resampling_approximation",
        "control_group_complement",
        "run_permutations",
        "grna_integration_strategy",
        "multiple_testing_alpha",
        "multiple_testing_method"
    )
    list(
        data = object@data@identity,
        cells = object@preparation@analysis_cells,
        design = digest::digest(
            object@preparation@design_matrix,
            algo = "sha256"
        ),
        assignments = digest::digest(
            object@preparation@analysis_groups,
            algo = "sha256"
        ),
        choices = stats::setNames(
            lapply(options, function(name) methods::slot(object@spec, name)),
            options
        ),
        execution = digest::digest(run_plan_payload(plan), algo = "sha256")
    )
}

run_resampling_budget <- function(
    object,
    analysis,
    n_calibration_pairs = NULL
) {
    calibration <- identical(
        analysis_record_name(analysis),
        "run_calibration_check"
    )
    B1 <- 499L
    approximation <- identical(
        object@spec@resampling_approximation,
        "skew_normal"
    )
    B2 <- if (approximation) 4999L else 0L
    B3 <- if (approximation) {
        if (isTRUE(object@spec@run_permutations)) 24999L else 0L
    } else {
        # Preserve the established common discovery/power family policy.
        n <- if (calibration) {
            n_calibration_pairs
        } else {
            max(
                pair_count(object, "discovery", passing = TRUE),
                pair_count(object, "positive_control", passing = TRUE)
            )
        }
        compute_empirical_budget(
            n,
            object@spec@multiple_testing_alpha,
            object@spec@side_code,
            B1,
            B2
        )
    }
    c(B1 = B1, B2 = B2, B3 = B3)
}

# budgets is a private numerical-testing/extension seam. Public workflows always
# derive the established resampling policy from the resolved scientific choices.
make_run_plan <- function(
    object,
    analysis,
    pairs = NULL,
    n_calibration_pairs = NULL,
    calibration_group_size = NULL,
    budgets = NULL
) {
    analysis <- analysis_record_name(analysis)
    calibration <- identical(analysis, "run_calibration_check")
    if (is.null(pairs)) {
        pairs <- analysis_pairs(object, analysis)
    }
    if (calibration) {
        if (is.null(n_calibration_pairs)) {
            n_calibration_pairs <- nrow(pairs)
        }
        if (is.null(calibration_group_size)) {
            sizes <- lengths(strsplit(pairs$grna_group, "&", fixed = TRUE))
            if (!length(sizes) || length(unique(sizes)) != 1L) {
                stop("Calibration groups must have one common size.")
            }
            calibration_group_size <- sizes[[1L]]
        }
        check_native_scalar_count(
            n_calibration_pairs,
            "n_calibration_pairs",
            minimum = 1L
        )
        check_native_scalar_count(
            calibration_group_size,
            "calibration_group_size",
            minimum = 1L
        )
    } else if (
        !is.null(n_calibration_pairs) || !is.null(calibration_group_size)
    ) {
        stop(
            "Calibration-only RunPlan sizes must be absent for other analyses."
        )
    }
    if (is.null(budgets)) {
        budgets <- run_resampling_budget(object, analysis, n_calibration_pairs)
    }
    if (!all(c("B1", "B2", "B3") %in% names(budgets))) {
        stop("Run budgets require B1, B2 and B3.")
    }
    checked_resampling_total(budgets[["B1"]], budgets[["B2"]], budgets[["B3"]])
    plan <- list(
        analysis = analysis,
        pairs = copy_run_plan_pairs(pairs),
        B1 = as.integer(budgets[["B1"]]),
        B2 = as.integer(budgets[["B2"]]),
        B3 = as.integer(budgets[["B3"]]),
        n_calibration_pairs = if (calibration) {
            as.integer(n_calibration_pairs)
        } else {
            NULL
        },
        calibration_group_size = if (calibration) {
            as.integer(calibration_group_size)
        } else {
            NULL
        }
    )
    plan$identity <- run_plan_identity(object, plan)
    validate_run_plan(object, as_run_plan(plan))
}

validate_run_plan <- function(object, plan) {
    plan <- as_run_plan(plan)
    # Structural validity belongs to RunPlan. Membership and matching actual
    # screen/design/cell/assignment inputs require this contextual boundary.
    pairs <- plan@pairs
    response_ids <- rownames(object@data@response_matrix)
    if (!all(pairs$response_id %in% response_ids)) {
        stop("Run-plan responses must belong to the imported screen.")
    }
    if (identical(plan@analysis, "run_calibration_check")) {
        guides <- unlist(
            strsplit(pairs$grna_group, "&", fixed = TRUE),
            use.names = FALSE
        )
        if (
            !all(
                guides %in%
                    names(object@preparation@analysis_groups$indiv_nt_grna_idxs)
            )
        ) {
            stop(
                "Calibration run-plan groups must use prepared non-targeting guides."
            )
        }
    } else if (
        !all(
            pairs$grna_group[pairs$pass_qc] %in%
                names(object@preparation@analysis_groups$grna_group_idxs)
        )
    ) {
        stop("Eligible run-plan guide groups must have prepared assignments.")
    }
    if (!identical(plan@identity, run_plan_identity(object, plan))) {
        stop(
            "The run plan no longer matches its data, design, assignments, choices or execution scope."
        )
    }
    invisible(plan)
}

result_input_identity <- function(object, analysis, plan = NULL) {
    analysis <- analysis_record_name(analysis)
    if (is.null(plan)) {
        plan <- object@results@records[[analysis]]$provenance$plan
    }
    if (is.null(plan)) {
        return(NULL)
    }
    plan <- tryCatch(as_run_plan(plan), error = function(error) NULL)
    if (is.null(plan)) {
        return(NULL)
    }
    if (!identical(analysis, "run_calibration_check")) {
        plan@pairs <- analysis_pairs(object, analysis)
    }
    run_plan_identity(object, plan)
}

result_record_is_current <- function(object, analysis) {
    record <- object@results@records[[analysis]]
    stage_complete(object, "run_qc") &&
        is.list(record) &&
        is.data.frame(record$table) &&
        nrow(record$table) > 0L &&
        identical(record$scope, "complete_family") &&
        is.null(record$invalidated_by) &&
        !is.null(record$identity) &&
        identical(record$identity, result_input_identity(object, analysis))
}

invalidate_result_records <- function(object, reason) {
    for (analysis in names(result_analysis_fields())) {
        if (!is.null(object@results@records[[analysis]])) {
            object@results@records[[analysis]]$invalidated_by <- reason
        }
    }
    object
}

# Check the complete finalized hypothesis family, including failed-QC rows and
# legitimate repeated requests. This uses labels only, never statistical kernels.
validate_result_family <- function(object, result, plan) {
    expected <- data.frame(
        response_id = plan@pairs$response_id,
        grna_group = plan@pairs$grna_group,
        stringsAsFactors = FALSE
    )
    strategy <- object@spec@grna_integration_strategy
    if (identical(strategy, "union")) {
        names(expected)[names(expected) == "grna_group"] <- "grna_target"
    } else {
        names(expected)[names(expected) == "grna_group"] <- "grna_id"
        expected <- dplyr::left_join(
            expected,
            object@data@grna_target_data_frame[c("grna_id", "grna_target")],
            by = "grna_id"
        )
        if (
            identical(strategy, "bonferroni") &&
                !identical(plan@analysis, "run_calibration_check")
        ) {
            expected <- unique(as.data.frame(expected[c(
                "response_id",
                "grna_target"
            )]))
        }
    }
    keys <- names(expected)
    if (!all(keys %in% names(result))) {
        stop("The result lacks complete-family identity columns.")
    }
    canonical_keys <- function(table) {
        out <- as.data.frame(
            stats::setNames(
                lapply(keys, function(key) as.character(table[[key]])),
                keys
            ),
            stringsAsFactors = FALSE
        )
        out <- out[
            do.call(order, c(unname(out), list(na.last = TRUE))),
            ,
            drop = FALSE
        ]
        rownames(out) <- NULL
        out
    }
    if (!identical(canonical_keys(expected), canonical_keys(result))) {
        stop(
            "The result does not cover the complete run-plan hypothesis family."
        )
    }
    invisible(NULL)
}

commit_result_record <- function(
    object,
    analysis,
    result,
    pairs,
    output_amount,
    plan = NULL
) {
    analysis <- analysis_record_name(analysis)
    require_stage(object, "run_qc")
    if (
        !is.data.frame(result) || !nrow(result) || !"p_value" %in% names(result)
    ) {
        stop(
            "A completed analysis requires a nonempty result table with p-values."
        )
    }
    if (is.null(plan)) {
        plan <- make_run_plan(object, analysis, pairs)
    }
    plan <- validate_run_plan(object, plan)
    if (!identical(plan@analysis, analysis) || !identical(plan@pairs, pairs)) {
        stop("The completed result does not match its run plan.")
    }
    if (
        !identical(analysis, "run_calibration_check") &&
            !identical(plan@pairs, analysis_pairs(object, analysis))
    ) {
        stop("Only the complete prepared hypothesis family can be committed.")
    }
    validate_result_family(object, result, plan)
    p <- result$p_value
    if (
        !is.numeric(p) ||
            any(is.nan(p)) ||
            any(is.infinite(p)) ||
            any(p < 0 | p > 1, na.rm = TRUE)
    ) {
        stop(
            "Completed result p-values must be numeric probabilities or missing."
        )
    }
    if (analysis %in% c("run_calibration_check", "run_discovery_analysis")) {
        expected_significant <- stats::p.adjust(
            p,
            method = object@spec@multiple_testing_method
        ) <
            object@spec@multiple_testing_alpha
        if (
            !"significant" %in% names(result) ||
                !identical(
                    unname(result$significant),
                    unname(expected_significant)
                )
        ) {
            stop(
                "Discovery and calibration results must be adjusted across the complete hypothesis family."
            )
        }
    }
    previous <- object@results@records[[analysis]]
    history <- if (is.null(previous)) {
        list()
    } else {
        c(previous$history, list(previous[setdiff(names(previous), "history")]))
    }
    object@results@records[[analysis]] <- owned_result_value(list(
        table = result,
        scope = "complete_family",
        identity = plan@identity,
        provenance = list(
            schema_version = object@schema_version,
            package_version = as.character(utils::packageVersion("sceptre")),
            spec = analysis_spec_snapshot(object),
            plan = plan,
            output_amount = output_amount,
            response_fit_method = response_fit_method_id(),
            perturbation_fit_method = perturbation_fit_method_id()
        ),
        history = history
    ))
    object
}

current_result <- function(object, analysis, required = TRUE) {
    if (!result_record_is_current(object, analysis)) {
        if (!required) {
            return(data.frame())
        }
        stop(
            "No completed current result for `",
            analysis,
            "`. Rerun the analysis or use get_result(..., historical = TRUE) to inspect saved results."
        )
    }
    object@results@records[[analysis]]$table
}
