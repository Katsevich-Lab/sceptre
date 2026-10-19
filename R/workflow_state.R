# A preparation stage exists only after successful validation. Results are
# completed records with their own scientific input identity. No call journal or
# duplicate workflow flags are saved in the public coordinator.
# Formula environments are execution contexts, not stable value identities. The
# actual compiled design and durable assignment values are fingerprinted below.
identity_value <- function(x) {
    if (inherits(x, "formula")) {
        return(list(formula = paste(deparse(x), collapse = "")))
    }
    if (is.list(x) && !is.data.frame(x)) {
        return(lapply(x, identity_value))
    }
    x
}

preparation_stage_names <- function() {
    c("design", "assignment", "cell_qc", "pair_qc")
}

stage_input_identity <- function(object, stage) {
    p <- object@preparation
    s <- object@spec
    inputs <- switch(
        stage,
        design = list(
            object@data@identity,
            s@formula_object,
            s@grna_integration_strategy,
            s@discovery_pairs,
            s@positive_control_pairs,
            p@design_matrix,
            p@guide_groups,
            p@discovery_pairs,
            p@positive_control_pairs
        ),
        assignment = list(
            object@data@identity,
            s@grna_assignment_method,
            s@grna_assignment_hyperparameters,
            p@guide_groups,
            p@guide_assignments,
            p@source_groups,
            p@assignment_exclusions
        ),
        cell_qc = list(
            object@data@identity,
            p@stages$assignment$identity,
            s@cellwise_qc_thresholds,
            s@additional_cells_to_remove,
            p@analysis_cells,
            p@cell_removal_metrics,
            p@analysis_groups
        ),
        pair_qc = list(
            object@data@identity,
            p@stages$cell_qc$identity,
            s@grna_integration_strategy,
            s@n_nonzero_trt_thresh,
            s@n_nonzero_cntrl_thresh,
            p@discovery_pairs,
            p@positive_control_pairs,
            p@discovery_pairs_with_info,
            p@positive_control_pairs_with_info,
            p@nt_nonzero_matrix,
            p@control_nonzero_totals
        ),
        stop("Unknown preparation stage: ", stage)
    )
    digest::digest(identity_value(inputs), algo = "sha256")
}

stage_complete <- function(object, stage) {
    if (stage == "import_data") {
        return(length(object@data@identity) == 1L)
    }
    stage <- switch(
        stage,
        set_analysis_parameters = "design",
        assign_grnas = "assignment",
        stage
    )
    if (stage == "run_qc") {
        return(
            stage_complete(object, "cell_qc") &&
                stage_complete(object, "pair_qc")
        )
    }
    if (!stage %in% preparation_stage_names()) {
        stop("Unknown preparation stage: ", stage)
    }
    record <- object@preparation@stages[[stage]]
    if (!is.list(record) || is.null(record$identity)) {
        return(FALSE)
    }
    dependency <- switch(
        stage,
        design = "import_data",
        assignment = "design",
        cell_qc = "assignment",
        pair_qc = "cell_qc"
    )
    stage_complete(object, dependency) &&
        identical(record$identity, stage_input_identity(object, stage))
}

require_stage <- function(object, stage) {
    if (!stage_complete(object, stage)) {
        public <- switch(
            stage,
            design = "set_analysis_parameters",
            assignment = "assign_grnas",
            cell_qc = "run_qc",
            pair_qc = "run_qc",
            stage
        )
        stop(
            "The function `",
            public,
            "()` must be called successfully before this stage."
        )
    }
    invisible(NULL)
}

begin_preparation_stage <- function(object, stage) {
    stage <- switch(
        stage,
        set_analysis_parameters = "design",
        assign_grnas = "assignment",
        run_qc = "cell_qc",
        stage
    )
    dependency <- switch(
        stage,
        design = "import_data",
        assignment = "design",
        cell_qc = "assignment",
        pair_qc = "cell_qc"
    )
    if (is.null(dependency)) {
        stop("Unknown preparation stage: ", stage)
    }
    require_stage(object, dependency)
    downstream <- switch(
        stage,
        design = preparation_stage_names(),
        assignment = c("assignment", "cell_qc", "pair_qc"),
        cell_qc = c("cell_qc", "pair_qc"),
        pair_qc = "pair_qc"
    )
    object@preparation@stages[downstream] <- NULL
    invalidate_result_records(object, paste(stage, "inputs changed"))
}

commit_preparation_stage <- function(object, stage) {
    validObject(object)
    p <- object@preparation
    if (stage == "design") {
        if (
            nrow(p@design_matrix) != ncol(get_response_matrix(object)) ||
                length(p@guide_groups) !=
                    nrow(object@data@grna_target_data_frame)
        ) {
            stop("Design and guide mapping must align with imported data.")
        }
    } else if (stage == "assignment") {
        if (
            !all(
                c("grna_group_idxs", "indiv_nt_grna_idxs") %in%
                    names(p@source_groups)
            )
        ) {
            stop(
                "Source guide groups must be prepared before assignment is committed."
            )
        }
    } else if (stage == "cell_qc") {
        if (
            !length(p@analysis_cells) ||
                !is.list(p@analysis_groups$grna_group_idxs)
        ) {
            stop(
                "Cell QC must produce retained cells and aligned guide groups."
            )
        }
    } else if (stage == "pair_qc") {
        for (type in c("discovery", "positive_control")) {
            pairs <- methods::slot(p, paste0(type, "_pairs_with_info"))
            if (
                nrow(pairs) != nrow(methods::slot(p, paste0(type, "_pairs"))) ||
                    (nrow(pairs) &&
                        (!all(
                            c(
                                "response_id",
                                "grna_group",
                                "n_nonzero_trt",
                                "n_nonzero_cntrl",
                                "pass_qc"
                            ) %in%
                                names(pairs)
                        ) ||
                            anyNA(pairs$pass_qc)))
            ) {
                stop(
                    "Pair QC must cover the complete prepared hypothesis family."
                )
            }
        }
    } else {
        stop("Unknown preparation stage: ", stage)
    }
    object@preparation@stages[[stage]] <- list(
        identity = stage_input_identity(object, stage)
    )
    object
}

workflow_status <- function(object) {
    preparation <- c(
        import_data = stage_complete(object, "import_data"),
        set_analysis_parameters = stage_complete(object, "design"),
        assign_grnas = stage_complete(object, "assignment"),
        run_qc = stage_complete(object, "run_qc")
    )
    analyses <- names(result_analysis_fields())
    c(
        preparation,
        stats::setNames(
            vapply(
                analyses,
                function(name) result_record_is_current(object, name),
                logical(1)
            ),
            analyses
        )
    )
}
