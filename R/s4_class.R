# Domain state only. RunPlan and CellView are temporary typed values.
#' @import Matrix
#' @importFrom methods new slot validObject
setClassUnion("CountHandle", c("externalptr", "NULL"))
setClass(
    "BackedCountStore",
    slots = c(descriptor = "list", handle = "CountHandle"),
    prototype = list(handle = NULL)
)
setClassUnion("CountStore", c("dgRMatrix", "BackedCountStore", "NULL"))

setClass(
    "ScreenData",
    slots = c(
        response_names = "character",
        response_matrix = "CountStore",
        grna_matrix = "CountStore",
        covariate_data_frame = "data.frame",
        grna_target_data_frame = "data.frame",
        low_moi = "logical",
        import_grna_assignment_info = "list",
        identity = "character"
    )
)
setClass(
    "AnalysisSpec",
    slots = c(
        discovery_pairs = "data.frame",
        positive_control_pairs = "data.frame",
        formula_object = "formula",
        side_code = "integer",
        resampling_approximation = "character",
        control_group_complement = "logical",
        run_permutations = "logical",
        n_nonzero_trt_thresh = "integer",
        n_nonzero_cntrl_thresh = "integer",
        grna_integration_strategy = "character",
        grna_assignment_method = "character",
        grna_assignment_hyperparameters = "list",
        multiple_testing_alpha = "numeric",
        multiple_testing_method = "character",
        cellwise_qc_thresholds = "list",
        additional_cells_to_remove = "integer"
    )
)
setClass(
    "PreparationRecord",
    slots = c(
        design_matrix = "matrix",
        nt_nonzero_matrix = "matrix",
        control_nonzero_totals = "integer",
        discovery_pairs_with_info = "data.frame",
        positive_control_pairs_with_info = "data.frame",
        guide_assignments = "list",
        source_groups = "list",
        analysis_groups = "list",
        guide_counts_per_cell = "integer",
        assignment_exclusions = "integer",
        analysis_cells = "integer",
        cell_removal_metrics = "integer",
        stages = "list",
        guide_groups = "character",
        discovery_pairs = "data.frame",
        positive_control_pairs = "data.frame"
    )
)
setClass(
    "FitCache",
    slots = c(
        response_fits = "list",
        response_identity = "list",
        perturbation_fits = "list",
        perturbation_identity = "list"
    )
)
setClass("ResultRecord", slots = c(records = "list", import_metadata = "list"))

setClass(
    "RunPlan",
    slots = c(
        analysis = "character",
        pairs = "data.frame",
        B1 = "integer",
        B2 = "integer",
        B3 = "integer",
        n_calibration_pairs = "integer",
        calibration_group_size = "integer",
        identity = "list"
    )
)

setClass(
    "CellView",
    slots = c(
        fit_to_analysis = "integer",
        source_cells = "integer",
        trt_fit_pos = "integer",
        n_trt = "integer"
    )
)
# The coordinator owns workflow transitions and references, never duplicate data.
setClass(
    "sceptre_object",
    slots = c(
        data = "ScreenData",
        spec = "AnalysisSpec",
        preparation = "PreparationRecord",
        fits = "FitCache",
        results = "ResultRecord",
        schema_version = "integer"
    ),
    prototype = list(
        data = new("ScreenData"),
        spec = new("AnalysisSpec"),
        preparation = new("PreparationRecord"),
        fits = new("FitCache"),
        results = new("ResultRecord"),
        schema_version = 3L
    )
)
