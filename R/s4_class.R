# class union
# setClassUnion("grna_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix", "dgTMatrix", "lgCMatrix", "lgRMatrix", "lgTMatrix"))
# setClassUnion("response_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix", "dgTMatrix"))


#' The sceptre_object S4 class
#'
#' @import Matrix
#'
#' @description
#' An object of class \code{sceptre_object} holds the inputs, parameters,
#' intermediate computations, and results of a single-cell CRISPR screen
#' analysis. The object is built up incrementally across the pipeline
#' (\code{\link{import_data}} \eqn{\rightarrow}
#' \code{\link{set_analysis_parameters}} \eqn{\rightarrow}
#' \code{\link{assign_grnas}} \eqn{\rightarrow} \code{\link{run_qc}}
#' \eqn{\rightarrow} \code{\link{run_calibration_check}} /
#' \code{\link{run_power_check}} /
#' \code{\link{run_discovery_analysis}}). Slots not yet populated by the
#' current pipeline stage are intentionally left at their length-0 defaults.
#'
#' Users should construct and modify \code{sceptre_object}s through the
#' exported pipeline functions rather than by calling
#' \code{new("sceptre_object")} directly.
#'
#' @slot response_matrix list holding the response (e.g. gene-by-cell
#'   expression) matrix.
#' @slot grna_matrix list holding the gRNA-by-cell count matrix.
#' @slot covariate_data_frame data frame of per-cell covariates.
#' @slot covariate_matrix numeric model matrix derived from
#'   \code{covariate_data_frame} and \code{formula_object}.
#' @slot grna_target_data_frame data frame mapping each gRNA (\code{grna_id})
#'   to its target (\code{grna_target}); both columns are character.
#' @slot low_moi logical flag: \code{TRUE} for low multiplicity-of-infection
#'   data, \code{FALSE} for high MOI.
#' @slot covariate_names character vector of covariate column names
#'   (sorted) used for analysis.
#' @slot discovery_pairs data frame of (\code{grna_target},
#'   \code{response_id}) pairs to test for discovery.
#' @slot positive_control_pairs data frame of positive-control pairs.
#' @slot formula_object the regression \code{formula} used in the analysis.
#' @slot side_code integer encoding test sidedness: \code{-1L} (left),
#'   \code{0L} (both), or \code{1L} (right).
#' @slot resampling_approximation character; one of \code{"skew_normal"} or
#'   \code{"no_approximation"}.
#' @slot control_group_complement logical; if \code{TRUE}, the control group
#'   is the complement of the treatment group.
#' @slot run_permutations logical; \code{TRUE} for permutation resampling,
#'   \code{FALSE} for the conditional randomization test.
#' @slot n_nonzero_trt_thresh integer minimum number of non-zero treatment
#'   cells required to test a pair.
#' @slot n_nonzero_cntrl_thresh integer minimum number of non-zero control
#'   cells required to test a pair.
#' @slot B1,B2,B3 integer numbers of resamples used at each stage of the
#'   resampling scheme.
#' @slot grna_integration_strategy character strategy for integrating gRNAs
#'   targeting the same site (e.g. \code{"union"}, \code{"bonferroni"},
#'   \code{"singleton"}).
#' @slot grna_assignment_method character method for assigning gRNAs to
#'   cells (e.g. \code{"mixture"}, \code{"maximum"}, \code{"thresholding"}).
#' @slot grna_assignment_hyperparameters list of hyperparameters for
#'   \code{grna_assignment_method}.
#' @slot multiple_testing_alpha numeric significance level in \eqn{(0, 1)}.
#' @slot multiple_testing_method character multiple-testing correction
#'   method (e.g. \code{"BH"}).
#' @slot cell_removal_metrics integer vector summarising cells removed at
#'   QC.
#' @slot cellwise_qc_thresholds list of per-cell QC thresholds applied in
#'   \code{\link{run_qc}}.
#' @slot mitochondrial_gene logical vector flagging mitochondrial responses.
#' @slot M_matrix numeric matrix of design-matrix projections used in
#'   resampling.
#' @slot n_nonzero_tot_vector integer per-response count of non-zero cells.
#' @slot discovery_pairs_with_info data frame of discovery pairs annotated
#'   with per-pair diagnostics.
#' @slot positive_control_pairs_with_info data frame of positive-control
#'   pairs annotated with per-pair diagnostics.
#' @slot negative_control_pairs data frame of negative-control pairs
#'   constructed from non-targeting gRNAs.
#' @slot initial_grna_assignment_list list of initial (pre-QC) gRNA
#'   assignments per gRNA.
#' @slot grna_assignments_raw list of raw (post-assignment, pre-integration)
#'   gRNA assignments.
#' @slot grna_assignments list of gRNA assignments used for analysis (after
#'   integration).
#' @slot grnas_per_cell integer vector of the number of gRNAs assigned to
#'   each cell.
#' @slot cells_w_zero_or_twoplus_grnas integer indices of cells with zero or
#'   two-plus gRNA assignments (low-MOI only).
#' @slot cells_in_use integer indices of cells retained after QC.
#' @slot n_discovery_pairs,n_ok_discovery_pairs integer counts of total and
#'   QC-passing discovery pairs.
#' @slot n_positive_control_pairs,n_ok_positive_control_pairs integer
#'   counts of total and QC-passing positive-control pairs.
#' @slot calibration_group_size integer size of each calibration group.
#' @slot n_calibration_pairs integer number of calibration pairs.
#' @slot import_grna_assignment_info list of provisional gRNA-assignment
#'   information carried over from \code{import_data}.
#' @slot mean_cells_per_grna numeric mean number of cells per gRNA.
#' @slot response_precomputations list of cached per-response fits used to
#'   accelerate resampling.
#' @slot last_function_called character; the most recently completed
#'   pipeline step. One of \code{import_data},
#'   \code{set_analysis_parameters}, \code{assign_grnas}, \code{run_qc},
#'   \code{run_calibration_check}, \code{run_power_check},
#'   \code{run_discovery_analysis}.
#' @slot functs_called named logical vector recording which pipeline
#'   functions have run, with names matching the steps listed under
#'   \code{last_function_called}.
#' @slot calibration_result data frame of calibration-check results.
#' @slot power_result data frame of power-check results.
#' @slot discovery_result data frame of discovery-analysis results.
#' @slot nuclear logical flag enabling a stripped-down (\dQuote{nuclear})
#'   analysis path.
#' @slot nf_pipeline logical flag: \code{TRUE} when running inside the
#'   Nextflow pipeline.
#' @slot integer_id integer identifier used by the on-disk
#'   (\code{ondisc}) backend.
#' @slot elements_to_analyze character vector naming which pipeline
#'   elements to analyze in the on-disk / Nextflow setting.
#'
#' @name sceptre_object-class
#' @aliases sceptre_object
#' @docType class
#' @keywords classes
setClass("sceptre_object",
  slots = list(
    # raw data
    response_matrix = "list",
    grna_matrix = "list",
    covariate_data_frame = "data.frame",
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
    B1 = "integer", B2 = "integer", B3 = "integer",
    grna_integration_strategy = "character",
    grna_assignment_method = "character",
    grna_assignment_hyperparameters = "list",
    multiple_testing_alpha = "numeric",
    multiple_testing_method = "character",
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
