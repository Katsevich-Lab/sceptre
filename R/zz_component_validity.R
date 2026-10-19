# CellView is transient: no counts are copied into it. The existing numerical
# drivers use dollar access; all fields still have typed S4 ownership.
#' Internal analysis cell selection
#'
#' Small temporary mappings between source, analysis and fitting cell positions.
#' Dollar access exposes typed fields to the internal numerical drivers.
#' @name CellView-class
#' @aliases $,CellView-method
#' @keywords internal
NULL

setMethod("$", "CellView", function(x, name) slot(x, name))
setValidity("CellView", function(object) {
    if (length(object@source_cells) != length(object@fit_to_analysis) ||
        anyNA(c(object@source_cells, object@fit_to_analysis, object@trt_fit_pos)) ||
        any(object@source_cells < 1L) || any(object@fit_to_analysis < 1L) ||
        anyDuplicated(object@source_cells) || anyDuplicated(object@fit_to_analysis) ||
        any(object@trt_fit_pos < 1L | object@trt_fit_pos > length(object@source_cells)) ||
        !identical(object@n_trt, as.integer(length(object@trt_fit_pos)))) {
        return("CellView has inconsistent ordered source, analysis or treatment positions")
    }
    TRUE
})

setValidity("ScreenData", function(object) {
    if (!length(object@response_matrix) || !length(object@grna_matrix)) return(TRUE)
    response <- object@response_matrix[[1L]]
    guide <- object@grna_matrix[[1L]]
    if (is.null(dim(response)) || is.null(dim(guide))) return(TRUE) # closed backing
    if (ncol(response) != ncol(guide)) return("Modalities must have aligned cells")
    if (!is.null(colnames(response)) && !is.null(colnames(guide)) &&
        !identical(colnames(response), colnames(guide))) return("Modality cell order differs")
    if (nrow(object@covariate_data_frame) && nrow(object@covariate_data_frame) != ncol(response)) {
        return("Cell covariates must match the source cells")
    }
    TRUE
})

setValidity("sceptre_object", function(object) {
    if (!identical(object@schema_version, 2L)) return("Unsupported sceptre schema version")
    valid <- validObject(object@data, test = TRUE)
    if (!isTRUE(valid)) return(valid)
    if (length(object@data@response_matrix) &&
        !is.null(dim(object@data@response_matrix[[1L]]))) {
        n <- ncol(object@data@response_matrix[[1L]])
        cells <- object@preparation@cells_in_use
        if (anyNA(cells) || anyDuplicated(cells) || any(cells < 1L | cells > n)) {
            return("Prepared cells must be unique source-cell positions")
        }
        design <- object@preparation@covariate_matrix
        if (length(design) && nrow(design) != n) return("Design must match source cells")
        assignments <- object@preparation@initial_grna_assignment_list
        if (length(assignments) && any(vapply(assignments, function(x) {
            !is.numeric(x) || anyNA(x) || any(x != as.integer(x)) ||
                any(x < 1L | x > n) || anyDuplicated(x) > 0L
        }, logical(1)))) return("Assignments must use canonical source-cell coordinates")
    }
    TRUE
})
