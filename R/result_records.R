# Result tables remain available after upstream edits, but only a completed run
# with matching inputs may be consumed as the current analysis.
result_analysis_fields <- function() c(
    run_calibration_check = "calibration_result",
    run_power_check = "power_result",
    run_discovery_analysis = "discovery_result"
)

result_input_identity <- function(object, analysis) {
    pair_field <- switch(analysis,
        run_calibration_check = "negative_control_pairs",
        run_power_check = "positive_control_pairs_with_info",
        run_discovery_analysis = "discovery_pairs_with_info")
    pairs <- if (analysis == "run_calibration_check") {
        object@spec@negative_control_pairs
    } else methods::slot(object@preparation, pair_field)
    option_names <- c("side_code", "resampling_approximation",
        "control_group_complement", "run_permutations", "B1", "B2",
        "grna_integration_strategy", "multiple_testing_alpha",
        "multiple_testing_method")
    # B3 is a working budget recomputed for each hypothesis family. A later
    # calibration may change it without changing a completed discovery run.
    list(data = object@data@identity,
        cells = object@preparation@cells_in_use,
        design = digest::digest(object@preparation@covariate_matrix, algo = "sha256"),
        assignments = digest::digest(object@preparation@grna_assignments, algo = "sha256"),
        pairs = digest::digest(pairs, algo = "sha256"),
        options = stats::setNames(lapply(option_names, function(name) {
            methods::slot(object@spec, name)
        }), option_names))
}

result_record_is_current <- function(object, analysis) {
    record <- object@results@records[[analysis]]
    isTRUE(object@functs_called[[analysis]]) &&
        is.list(record) && isTRUE(record$complete) &&
        identical(record$status, "current") &&
        identical(record$identity, result_input_identity(object, analysis))
}

invalidate_result_records <- function(object, reason) {
    for (analysis in names(result_analysis_fields())) {
        if (!is.null(object@results@records[[analysis]])) {
            object@results@records[[analysis]]$status <- "historical"
            object@results@records[[analysis]]$invalidation_reason <- reason
        }
        object@functs_called[analysis] <- FALSE
    }
    object
}

commit_result_record <- function(object, analysis, result, pairs, output_amount) {
    if (!is.data.frame(result) || !nrow(result) || !"p_value" %in% names(result)) {
        stop("A completed analysis requires a nonempty result table with p-values.")
    }
    field <- result_analysis_fields()[[analysis]]
    previous <- object@results@records[[analysis]]
    history <- if (is.null(previous)) list() else c(
        previous$history, list(c(previous[setdiff(names(previous), "history")],
            list(result = methods::slot(object@results, field)))))
    methods::slot(object@results, field, check = FALSE) <- result
    object@results@records[[analysis]] <- list(
        complete = TRUE, status = "current",
        scope = if (isTRUE(object@nf_pipeline)) "raw_shard" else "complete_family",
        identity = result_input_identity(object, analysis),
        provenance = list(schema_version = object@schema_version,
            package_version = as.character(utils::packageVersion("sceptre")),
            requested = object@spec@requested,
            effective = list(B1 = object@spec@B1, B2 = object@spec@B2,
                B3 = object@spec@B3,
                alpha = object@spec@multiple_testing_alpha,
                correction = object@spec@multiple_testing_method,
                output_amount = output_amount,
                assignment = list(method = object@spec@grna_assignment_method,
                    hyperparameters = object@spec@grna_assignment_hyperparameters),
                qc = list(cellwise = object@spec@cellwise_qc_thresholds,
                    treatment_nonzero = object@spec@n_nonzero_trt_thresh,
                    control_nonzero = object@spec@n_nonzero_cntrl_thresh),
                analysis = result_input_identity(object, analysis)$options,
                response_fit_method = "poisson-glm.fit-then-dispersion-v1"),
            pairs = pairs), history = history)
    object@functs_called[analysis] <- TRUE
    object@last_function_called <- analysis
    object
}

current_result <- function(object, analysis, required = TRUE) {
    if (!result_record_is_current(object, analysis)) {
        if (!required) return(data.frame())
        stop("No completed current result for `", analysis,
            "`. Rerun the analysis or use get_result(..., historical = TRUE) to inspect saved results.")
    }
    record <- object@results@records[[analysis]]
    if (!identical(record$scope, "complete_family")) {
        if (!required) return(data.frame())
        stop("This result contains raw shards; complete-family adjustment is required before reporting.")
    }
    methods::slot(object@results, result_analysis_fields()[[analysis]])
}
