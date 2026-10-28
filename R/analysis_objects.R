# First draft of the Screen and Analysis objects (sceptre development plan,
# slide 5). The Screen is a SingleCellExperiment holding the imported data; an
# Analysis holds one analysis of a Screen: its specification, the records of
# its stages, its result tables and its reusable model fits.
#
# For now an Analysis is converted to and from the existing sceptre_object, so
# every analysis function keeps working unchanged. analysis_layout says where
# each sceptre_object slot lives.

analysis_layout <- list(
    screen = c(
        "response_matrix",
        "grna_matrix",
        "covariate_data_frame",
        "grna_target_data_frame",
        "low_moi",
        "import_grna_assignment_info"
    ),
    spec = c(
        "discovery_pairs",
        "positive_control_pairs",
        "formula_object",
        "side_code",
        "resampling_approximation",
        "control_group_complement",
        "run_permutations",
        "n_nonzero_trt_thresh",
        "n_nonzero_cntrl_thresh",
        "B1",
        "B2",
        "B3",
        "grna_integration_strategy",
        "grna_assignment_method",
        "grna_assignment_hyperparameters",
        "multiple_testing_alpha",
        "multiple_testing_method",
        "cellwise_qc_thresholds"
    ),
    stages = list(
        assignment = c(
            "initial_grna_assignment_list",
            "grna_assignments_raw",
            "grna_assignments",
            "grnas_per_cell",
            "cells_w_zero_or_twoplus_grnas",
            "mean_cells_per_grna"
        ),
        qc = c("cell_removal_metrics", "cells_in_use"),
        design = c("covariate_names", "covariate_matrix"),
        pairs = c(
            "M_matrix",
            "n_nonzero_tot_vector",
            "discovery_pairs_with_info",
            "positive_control_pairs_with_info",
            "n_discovery_pairs",
            "n_ok_discovery_pairs",
            "n_positive_control_pairs",
            "n_ok_positive_control_pairs"
        ),
        calibration = c(
            "calibration_group_size",
            "negative_control_pairs",
            "n_calibration_pairs"
        )
    ),
    results = c(
        calibration = "calibration_result",
        power = "power_result",
        discovery = "discovery_result"
    ),
    fits = "response_precomputations",
    status = c("last_function_called", "functs_called"),
    # unused, or used only by the on-disk (ondisc) pipeline
    legacy = c(
        "mitochondrial_gene",
        "nuclear",
        "nf_pipeline",
        "integer_id",
        "elements_to_analyze"
    )
)

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
NULL

setClass("AnalysisSpec", slots = c(settings = "list"))

# A stage's outputs and a hash of them, so a later stage can record exactly
# which version of its inputs it used.
setClass("StageRecord", slots = c(contents = "list", hash = "character"))

setClass(
    "Analysis",
    slots = c(
        screen = "SingleCellExperiment",
        spec = "AnalysisSpec",
        stages = "list",
        results = "list",
        fits = "list",
        status = "list",
        legacy = "list"
    )
)

stage_record <- function(contents) {
    methods::new(
        "StageRecord",
        contents = contents,
        hash = rlang::hash(contents)
    )
}

slot_values <- function(object, names) {
    stats::setNames(lapply(names, methods::slot, object = object), names)
}

#' Convert between a sceptre_object and an Analysis
#'
#' `as_analysis()` splits a `sceptre_object` into the components of the new
#' design: a `SingleCellExperiment` holding the imported data (the Screen), the
#' analysis specification, one record per analysis stage, the result tables and
#' the cached model fits. `as_sceptre_object()` converts back. The conversion
#' loses nothing, so the existing analysis functions can be used on either.
#'
#' @param sceptre_object a `sceptre_object`
#' @param analysis an `Analysis`
#' @return an `Analysis` or a `sceptre_object`
#' @keywords internal
as_analysis <- function(sceptre_object) {
    screen_slots <- slot_values(sceptre_object, analysis_layout$screen)
    screen <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = screen_slots$response_matrix[[1]]),
        colData = S4Vectors::DataFrame(
            screen_slots$covariate_data_frame,
            check.names = FALSE
        ),
        altExps = list(
            gRNA = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = screen_slots$grna_matrix[[1]]),
                # each cell's top gRNA and its share of the cell's gRNA UMIs
                colData = S4Vectors::DataFrame(
                    screen_slots$import_grna_assignment_info,
                    check.names = FALSE
                )
            )
        ),
        metadata = list(
            low_moi = screen_slots$low_moi,
            grna_target_data_frame = screen_slots$grna_target_data_frame,
            covariate_row_names = attr(
                screen_slots$covariate_data_frame,
                "row.names"
            )
        )
    )
    do.call(
        methods::new,
        c(list("Analysis", screen = screen), analysis_parts(sceptre_object))
    )
}

# Everything in an Analysis except the Screen.
analysis_parts <- function(sceptre_object) {
    list(
        spec = methods::new(
            "AnalysisSpec",
            settings = slot_values(sceptre_object, analysis_layout$spec)
        ),
        stages = lapply(analysis_layout$stages, function(names) {
            stage_record(slot_values(sceptre_object, names))
        }),
        results = lapply(
            analysis_layout$results,
            methods::slot,
            object = sceptre_object
        ),
        fits = methods::slot(sceptre_object, analysis_layout$fits),
        status = slot_values(sceptre_object, analysis_layout$status),
        legacy = slot_values(sceptre_object, analysis_layout$legacy)
    )
}

#' @rdname as_analysis
#' @keywords internal
as_sceptre_object <- function(analysis) {
    screen <- analysis@screen
    meta <- S4Vectors::metadata(screen)
    covariates <- as.data.frame(
        SummarizedExperiment::colData(screen),
        optional = TRUE
    )
    attr(covariates, "row.names") <- meta$covariate_row_names
    grna <- SingleCellExperiment::altExp(screen, "gRNA")
    values <- c(
        list(
            response_matrix = list(SummarizedExperiment::assay(
                screen,
                "counts"
            )),
            grna_matrix = list(SummarizedExperiment::assay(grna, "counts")),
            covariate_data_frame = covariates,
            grna_target_data_frame = meta$grna_target_data_frame,
            low_moi = meta$low_moi,
            import_grna_assignment_info = as.list(
                SummarizedExperiment::colData(grna)
            )
        ),
        analysis@spec@settings,
        unlist(
            lapply(unname(analysis@stages), methods::slot, name = "contents"),
            recursive = FALSE
        ),
        stats::setNames(analysis@results, analysis_layout$results),
        stats::setNames(list(analysis@fits), analysis_layout$fits),
        analysis@status,
        analysis@legacy
    )
    object <- methods::new("sceptre_object")
    for (name in names(values)) {
        methods::slot(object, name, check = FALSE) <- values[[name]]
    }
    object
}
