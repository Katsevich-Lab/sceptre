#' Data getter functions
#'
#' The data getter functions (i.e., `get_response_matrix()`,
#' `get_grna_matrix()`, `get_cell_covariates()`) return of a specified data
#' field from a `sceptre_object`.
#'
#' @param sceptre_object a `sceptre_object`
#'
#' @return
#' `get_response_matrix()` returns the response matrix contained within a
#' `sceptre_object`.
#'
#' `get_grna_matrix()` returns the gRNA matrix contained within a
#' `sceptre_object`.
#'
#' `get_cell_covariates()` returns the cell covariate data frame contained
#' within a `sceptre_object`.
#' @export
#' @examples
#' # 1. create a sceptre_object
#'
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#'
#' # 2. extract data fields from the sceptre_object
#' response_matrix <- get_response_matrix(sceptre_object)
#' grna_matrix <- get_grna_matrix(sceptre_object)
#' cell_covariates <- get_cell_covariates(sceptre_object)
get_response_matrix <- function(sceptre_object) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    store <- sceptre_object@data@response_matrix
    return(count_store_value(store))
}

#' Get gRNA matrix
#'
#' @rdname get_response_matrix
#' @export
get_grna_matrix <- function(sceptre_object) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    store <- sceptre_object@data@grna_matrix
    return(count_store_value(store))
}

#' Get cell covariates
#'
#' @rdname get_response_matrix
#' @export
get_cell_covariates <- function(sceptre_object) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    return(sceptre_object@data@covariate_data_frame)
}


#' Get gRNA assignments
#'
#' `get_grna_assignments()` returns the gRNA-to-cell assignments contained
#' within a `sceptre_object`. The output is a sparse logical matrix, with gRNAs
#' in the rows and cells in the columns. A given entry of the matrix is set to
#' `TRUE` if the given gRNA is assigned to the given cell (and `FALSE`
#' otherwise).
#'
#' When using the "maximum" assignment strategy, exactly one gRNA is assigned to
#' a given cell. In other words, each column of the gRNA-to-cell assignment
#' matrix contains exactly one TRUE entry.
#'
#' @param sceptre_object a `sceptre_object` that has had `assign_grnas()` called
#' on it
#' @param apply_cellwise_qc a logical value (i.e., `TRUE` or `FALSE`) indicating
#' whether to return the gRNA-to-cell assignment matrix after cellwise QC has
#' been applied (default `FALSE`)
#'
#' @return a sparse logical matrix containing the gRNA-to-cell assignments.
#'   With `apply_cellwise_qc = TRUE`, the existing Matrix subsetting behavior
#'   may return a logical vector when one retained dimension has length one.
#' @export
#' @examples
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' discovery_pairs <- construct_cis_pairs(sceptre_object)
#' sceptre_object <- sceptre_object |>
#'   set_analysis_parameters(
#'     discovery_pairs = discovery_pairs,
#'     side = "left"
#'   ) |>
#'   assign_grnas(method = "mixture") |>
#'   run_qc()
#'
#' grna_assignment_matrix <- get_grna_assignments(
#'   sceptre_object = sceptre_object
#' )
#' grna_assignment_matrix_with_qc <- get_grna_assignments(
#'   sceptre_object = sceptre_object,
#'   apply_cellwise_qc = TRUE
#' )
get_grna_assignments <- function(sceptre_object, apply_cellwise_qc = FALSE) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    if (!stage_complete(sceptre_object, "assignment")) {
        stop(
            "`assign_grnas()` has not yet been called on the `sceptre_object`."
        )
    }
    initial_grna_assignment_list <- sceptre_object@preparation@guide_assignments
    grna_ids <- names(initial_grna_assignment_list)
    names(initial_grna_assignment_list) <- NULL
    # Saved assignments are original cell columns starting at 1. Subtracting
    # one creates a new vector for Matrix's zero-based column numbers, leaving
    # the saved assignments unchanged.
    j <- unlist(initial_grna_assignment_list) - 1L
    l <- vapply(initial_grna_assignment_list, length, FUN.VALUE = integer(1))
    p <- c(0L, cumsum(l))
    mat <- Matrix::sparseMatrix(j = 1, p = c(0, 1), repr = "R")
    mat@j <- j
    mat@p <- p
    mat@Dim <- c(
        length(initial_grna_assignment_list),
        get_grna_matrix(sceptre_object) |> ncol()
    )
    rownames(mat) <- grna_ids
    if (apply_cellwise_qc) {
        if (!stage_complete(sceptre_object, "run_qc")) {
            stop("QC has not yet been called on this sceptre_object.")
        }
        # Keep Matrix's usual subsetting: a single remaining row/column may
        # become a vector rather than a matrix.
        mat <- mat[, sceptre_object@preparation@analysis_cells]
    }
    return(mat)
}

#' Replace the response matrix stored in the object
#'
#' Normalize a replacement count store and refresh its imported summaries.
#' Return the updated object after dependent preparation and fits are cleared.
#' @noRd
set_response_matrix <- function(sceptre_object, response_matrix) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    previous <- sceptre_object@data@response_matrix
    if (
        length(sceptre_object@data@identity) &&
            !inherits(response_matrix, "odm") &&
            !methods::is(response_matrix, "BackedCountStore") &&
            !is.null(dim(response_matrix))
    ) {
        response_matrix <- set_matrix_accessibility(
            response_matrix,
            make_row_accessible = TRUE
        )
    }
    if (
        length(previous) &&
            !identical(dimnames(previous), dimnames(response_matrix))
    ) {
        stop(
            "Replacement counts must preserve aligned cells and feature identities; use import_data() for a new screen."
        )
    }
    sceptre_object@data@response_matrix <- as_count_store(response_matrix)
    validObject(sceptre_object@data)
    if (
        length(sceptre_object@data@identity) &&
            !identical(previous, sceptre_object@data@response_matrix)
    ) {
        sceptre_object <- refresh_imported_summaries(sceptre_object, "response")
        sceptre_object <- invalidate_screen_dependents(sceptre_object)
    }
    return(sceptre_object)
}

#' Replace the guide matrix stored in the object
#'
#' Replace normalized guide counts and refresh the assignment summaries.
#' @noRd
set_grna_matrix <- function(sceptre_object, grna_matrix) {
    sceptre_object <- migrate_sceptre_object(sceptre_object)
    previous <- sceptre_object@data@grna_matrix
    if (
        length(sceptre_object@data@identity) &&
            !inherits(grna_matrix, "odm") &&
            !methods::is(grna_matrix, "BackedCountStore") &&
            !is.null(dim(grna_matrix))
    ) {
        grna_matrix <- set_matrix_accessibility(
            grna_matrix,
            make_row_accessible = TRUE
        )
    }
    if (
        length(previous) &&
            !identical(dimnames(previous), dimnames(grna_matrix))
    ) {
        stop(
            "Replacement counts must preserve aligned cells and feature identities; use import_data() for a new screen."
        )
    }
    sceptre_object@data@grna_matrix <- as_count_store(grna_matrix)
    validObject(sceptre_object@data)
    if (
        length(sceptre_object@data@identity) &&
            !identical(previous, sceptre_object@data@grna_matrix)
    ) {
        sceptre_object <- refresh_imported_summaries(sceptre_object, "grna")
        sceptre_object <- invalidate_screen_dependents(sceptre_object)
    }
    return(sceptre_object)
}

# Content replacement is atomic with its imported summaries. Backing relocation
# uses verified descriptors instead; replacing backed content requires reimport.
refresh_imported_summaries <- function(object, modality) {
    store <- methods::slot(object@data, paste0(modality, "_matrix"))
    if (methods::is(store, "BackedCountStore") || is.null(dim(store))) {
        stop(
            "Use import_data() to replace backed screen content and its import summaries."
        )
    }
    labels <- object@data@response_names
    if (
        modality == "response" &&
            !length(labels) &&
            "response_p_mito" %in% names(object@data@covariate_data_frame)
    ) {
        stop(
            "Original mitochondrial feature labels are unknown; use import_data() to replace response counts."
        )
    }
    summaries <- compute_cell_covariates(
        store,
        feature_names = if (modality == "response") labels else character(),
        compute_p_mito = modality == "response",
        compute_max_feature = modality == "grna"
    )
    if (modality == "grna") {
        object@data@import_grna_assignment_info <- list(
            max_grna = summaries$feature_w_max_expression,
            max_grna_frac_umis = summaries$frac_umis_max_feature
        )
        summaries$feature_w_max_expression <- summaries$frac_umis_max_feature <- NULL
    }
    for (field in names(summaries)) {
        object@data@covariate_data_frame[[paste0(
            modality,
            "_",
            field
        )]] <- summaries[[field]]
    }
    object
}
