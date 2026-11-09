# class union
#' @import Matrix
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
NULL
# setClassUnion("grna_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix",
# "dgTMatrix", "lgCMatrix", "lgRMatrix", "lgTMatrix"))
# setClassUnion("response_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix",
# "dgTMatrix"))

# package_version is an S3 class; register it so that it can be a slot class
methods::setOldClass(c("package_version", "numeric_version"))

#' The sceptre_object class
#'
#' A `sceptre_object` holds the data of a single-cell CRISPR screen together
#' with the parameters, intermediate values and results of its analysis. It is
#' created by [import_data()], [import_data_from_cellranger()],
#' [import_data_from_parse()] or [import_data_from_sce()], or read by
#' [read_sceptre_object()], and updated by the pipeline functions
#' [set_analysis_parameters()], [assign_grnas()], [run_qc()],
#' [run_calibration_check()], [run_power_check()] and
#' [run_discovery_analysis()]. Read a `sceptre_object` through the accessors
#' below rather than through its slots, which are internal and may change
#' between versions; [updateObject()] converts an object saved by an earlier
#' version.
#'
#' @section Accessors:
#' - [get_response_matrix()], [get_grna_matrix()] and [get_cell_covariates()]
#'   return the response UMI counts, the gRNA UMI counts and the cell
#'   covariates.
#' - [get_grna_assignments()] returns the gRNA-to-cell assignments.
#' - [get_result()] returns the results of the calibration check, the power
#'   check or the discovery analysis.
#' - `print()` summarizes the data and the analysis; `plot()` plots the step
#'   run last; [write_outputs_to_directory()] writes the results, plots and a
#'   summary to a directory.
#'
#' @slot version the version of `sceptre` that created the object, a
#' `package_version`; [updateObject()] sets it to `0.99.1` on an object
#' converted from an earlier version. An object created by a later version of
#' `sceptre` than the one installed is refused with a message to upgrade.
#' @slot data a `SingleCellExperiment` holding the response UMI counts (assay
#' `"counts"`), the gRNA UMI counts (assay `"counts"` of the alternative
#' experiment `"grna"`) and the cell covariates (`colData`)
#' @slot covariate_matrix the design matrix of the cell covariates given by
#' `formula_object`
#' @slot grna_target_data_frame the data frame mapping each gRNA to its target
#' (columns `grna_id` and `grna_target`, optionally `chr`, `start` and `end`)
#' @slot low_moi `TRUE` for a low-MOI screen, `FALSE` for a high-MOI screen
#' @slot covariate_names the names of the cell covariates
#' @slot discovery_pairs the discovery pairs set in [set_analysis_parameters()]
#' @slot positive_control_pairs the positive control pairs set in
#' [set_analysis_parameters()]
#' @slot formula_object the formula of the cell covariates adjusted for
#' @slot side_code the sidedness of the test (`-1L` left, `0L` both, `1L`
#' right)
#' @slot resampling_approximation the approximation to the null distribution
#' @slot control_group_complement `TRUE` if the control group is the complement
#' set, `FALSE` if it is the non-targeting cells
#' @slot run_permutations `TRUE` if the resampling mechanism is permutations,
#' `FALSE` if it is conditional resampling
#' @slot B1 the number of resamples of the first round
#' @slot B2 the number of resamples of the second round
#' @slot B3 the number of resamples of the third round
#' @slot grna_integration_strategy the gRNA integration strategy
#' @slot multiple_testing_alpha the level of the multiple testing correction
#' @slot multiple_testing_method the multiple testing correction method
#' @slot seed the seed of the random draws of the analysis
#' @slot n_nonzero_trt_thresh the minimum number of treatment cells with
#' nonzero expression of a pair passing pairwise quality control
#' @slot n_nonzero_cntrl_thresh the minimum number of control cells with
#' nonzero expression of a pair passing pairwise quality control
#' @slot cellwise_qc_thresholds the cellwise quality-control thresholds set in
#' [run_qc()]
#' @slot grna_assignment_method the gRNA assignment method set in
#' [assign_grnas()]
#' @slot grna_assignment_hyperparameters the parameters of the gRNA assignment
#' method
#' @slot initial_grna_assignment_list for each gRNA, the indices of the cells
#' to which it is assigned
#' @slot grna_assignments_raw for each gRNA group, the indices of the cells to
#' which it is assigned
#' @slot grna_assignments `grna_assignments_raw` after cellwise quality control
#' @slot grnas_per_cell the number of gRNAs assigned to each cell
#' @slot cells_w_zero_or_twoplus_grnas the cells assigned zero or several gRNAs
#' @slot mean_cells_per_grna the mean number of cells to which a gRNA is
#' assigned
#' @slot import_grna_assignment_info for each cell, the gRNA with the most UMIs
#' and its share of the cell's gRNA UMIs, computed at import
#' @slot cell_removal_metrics the numbers of cells removed by each cellwise
#' quality-control filter
#' @slot cells_in_use the indices of the cells that pass cellwise quality
#' control
#' @slot mitochondrial_gene unused
#' @slot M_matrix the numbers of cells with nonzero expression of each response
#' among the cells of each non-targeting gRNA, computed in pairwise quality
#' control
#' @slot n_nonzero_tot_vector the number of cells with nonzero expression of
#' each response, computed in pairwise quality control
#' @slot discovery_pairs_with_info the discovery pairs with their pairwise
#' quality-control information
#' @slot positive_control_pairs_with_info the positive control pairs with their
#' pairwise quality-control information
#' @slot n_discovery_pairs the number of discovery pairs
#' @slot n_ok_discovery_pairs the number of discovery pairs that pass pairwise
#' quality control
#' @slot n_positive_control_pairs the number of positive control pairs
#' @slot n_ok_positive_control_pairs the number of positive control pairs that
#' pass pairwise quality control
#' @slot negative_control_pairs the negative control pairs of the calibration
#' check
#' @slot calibration_group_size the number of non-targeting gRNAs in each
#' negative control gRNA group
#' @slot n_calibration_pairs the number of negative control pairs
#' @slot response_precomputations for each response, the fitted coefficients
#' and dispersion of its regression on the cell covariates, reused across
#' analyses
#' @slot last_function_called the pipeline function called last
#' @slot functs_called for each pipeline function, whether it has been called
#' @slot calibration_result the results of the calibration check
#' @slot power_result the results of the power check
#' @slot discovery_result the results of the discovery analysis
#'
#' @exportClass sceptre_object
#' @name sceptre_object-class
#' @rdname sceptre_object-class
setClass(
    "sceptre_object",
    slots = list(
        # the version of sceptre that created the object, set where objects
        # are created (init_sceptre_object(), update_sceptre_object()); a
        # blank object has an empty version
        version = "package_version",

        # raw data: response counts in assay "counts", gRNA counts in assay
        # "counts" of altExp "grna", cell covariates in colData
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
        discovery_result = "data.frame"
    ),
    prototype = methods::prototype(version = package_version(character()))
)


setValidity("sceptre_object", function(object) {
    problems <- character()
    if (length(object@version) > 1L) {
        problems <- c(problems, "`version` must be a single package version.")
    }
    screen <- object@data
    if (length(SummarizedExperiment::assayNames(screen)) >= 1L) {
        if (!("counts" %in% SummarizedExperiment::assayNames(screen))) {
            problems <- c(
                problems,
                "`data` must hold the response counts in assay \"counts\"."
            )
        }
        if ("grna" %in% SingleCellExperiment::altExpNames(screen)) {
            grna_experiment <- SingleCellExperiment::altExp(screen, "grna")
            if (
                !("counts" %in%
                    SummarizedExperiment::assayNames(grna_experiment))
            ) {
                problems <- c(
                    problems,
                    paste0(
                        "`data` must hold the gRNA counts in assay \"counts\" ",
                        "of the alternative experiment \"grna\"."
                    )
                )
            }
        }
    }
    if (length(object@low_moi) > 1L || anyNA(object@low_moi)) {
        problems <- c(problems, "`low_moi` must be a single TRUE or FALSE.")
    }
    grna_target_data_frame <- object@grna_target_data_frame
    if (
        ncol(grna_target_data_frame) >= 1L &&
            !all(
                c("grna_id", "grna_target") %in%
                    colnames(grna_target_data_frame)
            )
    ) {
        problems <- c(
            problems,
            "`grna_target_data_frame` must contain columns `grna_id` and `grna_target`."
        )
    }
    cells_in_use <- object@cells_in_use
    if (
        length(cells_in_use) >= 1L &&
            (anyNA(cells_in_use) ||
                min(cells_in_use) < 1L ||
                max(cells_in_use) > ncol(screen))
    ) {
        problems <- c(
            problems,
            "`cells_in_use` must contain indices of cells in `data`."
        )
    }
    if (
        length(object@side_code) > 1L ||
            (length(object@side_code) == 1L &&
                !(object@side_code %in% c(-1L, 0L, 1L)))
    ) {
        problems <- c(problems, "`side_code` must be -1, 0 or 1.")
    }
    for (slot_name in c("B1", "B2", "B3")) {
        value <- methods::slot(object, slot_name)
        if (length(value) > 1L || anyNA(value) || any(value < 0L)) {
            problems <- c(
                problems,
                paste0("`", slot_name, "` must be a single nonnegative integer.")
            )
        }
    }
    alpha <- object@multiple_testing_alpha
    if (
        length(alpha) > 1L ||
            anyNA(alpha) ||
            any(alpha <= 0 | alpha >= 1)
    ) {
        problems <- c(
            problems,
            "`multiple_testing_alpha` must be a single number between 0 and 1."
        )
    }
    if (length(problems) == 0L) TRUE else problems
})


# Builds the SingleCellExperiment held in the `data` slot: the response counts
# in assay "counts", the gRNA counts in assay "counts" of altExp "grna" and the
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
        covariates <- S4Vectors::DataFrame(
            covariate_data_frame,
            check.names = FALSE
        )
        # DataFrame() drops row names stored as integers (for example, those
        # of a subset of a data frame); keep them, as character row names
        if (.row_names_info(covariate_data_frame) > 0L) {
            rownames(covariates) <- rownames(covariate_data_frame)
        }
        SummarizedExperiment::colData(screen) <- covariates
    }
    SingleCellExperiment::altExp(screen, "grna", withDimnames = FALSE) <-
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
#' are carried over unchanged; in sceptre 0.99.2, an analysis continued on the
#' updated object draws its random numbers as 0.99.1 and earlier did (the
#' fixed seed 4). The `version` slot of the updated object is `0.99.1`.
#' Current objects are returned unchanged.
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


# Stops, naming updateObject(), when `object` was saved by a version of sceptre
# whose class definition had no `data` slot (0.99.1 and earlier); such an
# object cannot be used until it is updated. Stops, asking to upgrade sceptre,
# when `object` was created by a later version of sceptre than the one
# installed.
check_object_version <- function(object) {
    if (!methods::.hasSlot(object, "data")) {
        stop(
            "This `sceptre_object` was saved by an earlier version of sceptre; ",
            "convert it with `updateObject()` before using it.",
            call. = FALSE
        )
    }
    installed <- utils::packageVersion("sceptre")
    if (
        methods::.hasSlot(object, "version") &&
            length(object@version) == 1L &&
            object@version > installed
    ) {
        stop(
            "This `sceptre_object` was created by sceptre ",
            as.character(object@version),
            ", later than the installed sceptre ",
            as.character(installed),
            "; upgrade sceptre to use it.",
            call. = FALSE
        )
    }
    invisible(object)
}


# The data frame slots of a sceptre_object may hold data.tables (pairs and
# results) or tibbles (for example, a gRNA target data frame read with readr, or
# that of the example data). S4 recognizes these as data frames only once
# data.table and tibble, which register their classes with setOldClass(), are
# loaded; tibble is loaded with dplyr. Loads both, so that an object read from
# a file in a fresh R session passes validObject().
load_data_frame_classes <- function() {
    loadNamespace("data.table")
    loadNamespace("dplyr")
    invisible(NULL)
}


update_sceptre_object <- function(object) {
    # BiocGenerics::updateObject() validates the object this returns
    load_data_frame_classes()
    old_slots <- attributes(object)
    if ("data" %in% names(old_slots)) {
        if (!("version" %in% names(old_slots))) {
            # an object with `data` but no `version` was created by a
            # pre-release 0.99.2; it is stamped 0.99.1 like the objects of
            # the released versions that lack the slot
            attr(object, "version") <- package_version("0.99.1")
        }
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
    # 0.99.1 and earlier drew their random numbers with the fixed seed 4
    if (isTRUE(updated@functs_called["set_analysis_parameters"])) {
        updated@seed <- 4L
    }
    # 0.99.1 and earlier recorded no version
    updated@version <- package_version("0.99.1")
    return(updated)
}
