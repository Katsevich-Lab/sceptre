# class union
#' @import Matrix
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
# setClassUnion("grna_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix",
# "dgTMatrix", "lgCMatrix", "lgRMatrix", "lgTMatrix"))
# setClassUnion("response_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix",
# "dgTMatrix"))

# sceptre object class
setClass(
    "sceptre_object",
    slots = list(
        # raw data: response counts in assay "counts", gRNA counts in assay
        # "counts" of altExp "gRNA", cell covariates in colData
        data = "SingleCellExperiment",
        covariate_matrix = "matrix",
        grna_target_data_frame = "data.frame",
        low_moi = "logical",
        covariate_names = "character",

        # analysis parameters
        discovery_pairs = "data.frame",
        positive_control_pairs = "data.frame",
        formula_object = "formula",
        side_code = "integer",
        resampling_approximation = "character",
        control_group_complement = "logical",
        run_permutations = "logical",
        n_nonzero_trt_thresh = "integer",
        n_nonzero_cntrl_thresh = "integer",
        B1 = "integer",
        B2 = "integer",
        B3 = "integer",
        grna_integration_strategy = "character",
        grna_assignment_method = "character",
        grna_assignment_hyperparameters = "list",
        multiple_testing_alpha = "numeric",
        multiple_testing_method = "character",
        seed = "integer",
        rng_protocol = "character",
        cell_removal_metrics = "integer",
        cellwise_qc_thresholds = "list",

        # computed objects
        mitochondrial_gene = "logical",
        M_matrix = "matrix",
        n_nonzero_tot_vector = "integer",
        discovery_pairs_with_info = "data.frame",
        positive_control_pairs_with_info = "data.frame",
        negative_control_pairs = "data.frame",
        initial_grna_assignment_list = "list",
        grna_assignments_raw = "list",
        grna_assignments = "list",
        grnas_per_cell = "integer",
        cells_w_zero_or_twoplus_grnas = "integer",
        cells_in_use = "integer",
        n_discovery_pairs = "integer",
        n_ok_discovery_pairs = "integer",
        n_positive_control_pairs = "integer",
        n_ok_positive_control_pairs = "integer",
        calibration_group_size = "integer",
        n_calibration_pairs = "integer",
        import_grna_assignment_info = "list",
        mean_cells_per_grna = "numeric",

        # cached objects
        response_precomputations = "list",

        # flags
        last_function_called = "character",
        functs_called = "logical",

        # results
        calibration_result = "data.frame",
        power_result = "data.frame",
        discovery_result = "data.frame",

        # odm-related entries
        nuclear = "logical",
        nf_pipeline = "logical",
        integer_id = "integer",
        elements_to_analyze = "character"
    )
)


# Builds the SingleCellExperiment held in the `data` slot: the response counts
# in assay "counts", the gRNA counts in assay "counts" of altExp "gRNA" and the
# cell covariates in colData. The matrices are stored as given.
new_screen <- function(
    response_matrix,
    grna_matrix,
    covariate_data_frame = NULL
) {
    screen <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = response_matrix)
    )
    if (!is.null(covariate_data_frame)) {
        SummarizedExperiment::colData(screen) <- S4Vectors::DataFrame(
            covariate_data_frame,
            check.names = FALSE
        )
    }
    SingleCellExperiment::altExp(screen, "gRNA", withDimnames = FALSE) <-
        SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = grna_matrix)
        )
    return(screen)
}


#' Update a sceptre_object
#'
#' `updateObject()` updates a `sceptre_object` created by an earlier version of
#' `sceptre` (0.99.1 or earlier, which held the response matrix, gRNA matrix
#' and cell covariates in separate slots) to the current class definition,
#' which holds them in a `SingleCellExperiment`. The analysis state and results
#' are carried over unchanged. Current objects are returned unchanged.
#'
#' @param object a `sceptre_object`
#' @param ... additional arguments (ignored)
#' @param verbose ignored
#'
#' @return a `sceptre_object` of the current class definition
#' @importFrom BiocGenerics updateObject
#' @exportMethod updateObject
#' @examples
#' data(lowmoi_example_data)
#' sceptre_object <- import_data(
#'   response_matrix = lowmoi_example_data$response_matrix,
#'   grna_matrix = lowmoi_example_data$grna_matrix,
#'   grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'   moi = "low"
#' )
#' sceptre_object <- updateObject(sceptre_object)
#' @name updateObject,sceptre_object-method
#' @aliases updateObject,sceptre_object-method
setMethod(
    "updateObject",
    signature = signature("sceptre_object"),
    function(object, ..., verbose = FALSE) {
        update_sceptre_object(object)
    }
)


update_sceptre_object <- function(object) {
    old_slots <- attributes(object)
    if ("data" %in% names(old_slots)) {
        return(object)
    }
    # 0.99.1 and earlier: the matrices in one-element lists, the covariates in
    # their own slot
    response_matrix <- old_slots$response_matrix
    grna_matrix <- old_slots$grna_matrix
    if (
        !is.list(response_matrix) ||
            length(response_matrix) != 1L ||
            length(dim(response_matrix[[1]])) != 2L ||
            !is.list(grna_matrix) ||
            length(grna_matrix) != 1L ||
            length(dim(grna_matrix[[1]])) != 2L
    ) {
        stop(
            "This `sceptre_object` does not contain its response and gRNA ",
            "matrices (for example, it was saved by ",
            "`write_ondisc_backed_sceptre_object()`) and cannot be updated. ",
            "Re-import the data with the current version of sceptre.",
            call. = FALSE
        )
    }
    updated <- methods::new("sceptre_object")
    for (slot_name in intersect(
        names(old_slots),
        methods::slotNames("sceptre_object")
    )) {
        methods::slot(updated, slot_name, check = FALSE) <- old_slots[[
            slot_name
        ]]
    }
    updated@data <- new_screen(
        response_matrix = response_matrix[[1]],
        grna_matrix = grna_matrix[[1]],
        covariate_data_frame = old_slots$covariate_data_frame
    )
    return(updated)
}
