#' Import complete raw counts from a Seurat object
#'
#' Use one unsplit counts layer each for the response and guide assays.
#' Assay5 and legacy Assay containers are supported through SeuratObject 5.
#' All assay cells must equal the complete object's cells in source order.
#' No normalization, joining, feature selection or automatic alignment occurs.
#' Supply the original full feature library; upstream filtering cannot be
#' undone.
#'
#' @param seurat_object A Seurat object with in-memory raw count assays.
#' @param grna_target_data_frame Complete guide metadata in assay row order.
#' @param moi Either `"low"` or `"high"`.
#' @param response_assay Name of the response assay, default `"RNA"`.
#' @param grna_assay Name of the guide assay, default `"gRNA"`.
#' @param covariates Names of additional per-cell metadata columns to include.
#' @param response_name_column Optional response feature metadata column of
#'   symbols used to identify mitochondrial genes. Otherwise use feature IDs.
#' @return A new sceptre_object, with original counts and full-library
#'   covariates.
#' @export
#' @examples
#' if (requireNamespace("SeuratObject", quietly = TRUE) &&
#'     utils::packageVersion("SeuratObject") >= "5.0.0") {
#'     data("lowmoi_example_data")
#'     x <- lowmoi_example_data
#'     ids <- paste0("guide", seq_len(nrow(x$grna_matrix)))
#'     rownames(x$grna_matrix) <- x$grna_target_data_frame$grna_id <- ids
#'     s <- SeuratObject::CreateSeuratObject(
#'         Matrix::Matrix(x$response_matrix, sparse = TRUE),
#'         min.cells = 0, min.features = 0)
#'     s[["gRNA"]] <- SeuratObject::CreateAssay5Object(
#'         Matrix::Matrix(x$grna_matrix, sparse = TRUE),
#'         min.cells = 0, min.features = 0)
#'     object <- import_from_seurat(s, x$grna_target_data_frame, "low")
#' }
import_from_seurat <- function(
    seurat_object,
    grna_target_data_frame,
    moi,
    response_assay = "RNA",
    grna_assay = "gRNA",
    covariates = character(),
    response_name_column = NULL
) {
    require_seurat_object(seurat_object)
    response <- seurat_count_layer(seurat_object, response_assay)
    grna <- seurat_count_layer(seurat_object, grna_assay)
    if (
        !is.data.frame(grna_target_data_frame) ||
            !identical(
                as.character(grna_target_data_frame$grna_id),
                rownames(grna)
            )
    ) {
        stop(
            "Guide metadata grna_id must match all guide rows exactly in order."
        )
    }
    cd <- seurat_object[[]]
    if (!identical(rownames(cd), colnames(response))) {
        stop("Cell metadata must match counts exactly in source order.")
    }
    import_data(
        response,
        grna,
        grna_target_data_frame,
        moi,
        extra_covariates = interop_extra_covariates(
            cd,
            covariates,
            colnames(response)
        ),
        response_names = seurat_response_names(
            seurat_object,
            response_assay,
            response_name_column
        )
    )
}

seurat_response_names <- function(object, assay_name, column) {
    if (is.null(column)) {
        return(NA_character_)
    }
    rd <- object[[assay_name]][[]]
    if (
        !is.character(column) ||
            length(column) != 1L ||
            is.na(column) ||
            !column %in% colnames(rd) ||
            !identical(rownames(rd), rownames(object[[assay_name]]))
    ) {
        stop(
            "response_name_column must name aligned response feature metadata."
        )
    }
    labels <- as.character(rd[[column]])
    if (anyNA(labels) || any(!nzchar(labels))) {
        stop("Response names must be nonmissing and nonempty.")
    }
    labels
}

#' Return current SCEPTRE pair results to a Seurat object
#'
#' The result is a copy with owned pair tables and provenance in
#' SeuratObject::Misc(result, slot = "sceptre"). Counts, cell metadata,
#' reductions and existing unrelated miscellaneous information are unchanged.
#' Existing sceptre miscellaneous data is rejected to avoid overwriting results.
#'
#' @param sceptre_object A sceptre_object with completed current results.
#' @param seurat_object A Seurat object with exactly the imported response
#'   counts.
#' @param analyses Analysis names accepted by get_result(), or short aliases.
#' @param response_assay Name of the original response assay.
#' @return A copied Seurat object with owned pair results and provenance.
#' @export
#' @examples
#' if (requireNamespace("SeuratObject", quietly = TRUE) &&
#'     utils::packageVersion("SeuratObject") >= "5.0.0") {
#'     data("lowmoi_example_data")
#'     x <- lowmoi_example_data
#'     ids <- paste0("guide", seq_len(nrow(x$grna_matrix)))
#'     rownames(x$grna_matrix) <- x$grna_target_data_frame$grna_id <- ids
#'     s <- SeuratObject::CreateSeuratObject(
#'         Matrix::Matrix(x$response_matrix, sparse = TRUE),
#'         min.cells = 0, min.features = 0)
#'     s[["gRNA"]] <- SeuratObject::CreateAssay5Object(
#'         Matrix::Matrix(x$grna_matrix, sparse = TRUE),
#'         min.cells = 0, min.features = 0)
#'     o <- import_from_seurat(s, x$grna_target_data_frame, "low")
#'     pairs <- data.frame(response_id = rownames(x$response_matrix)[1L],
#'         grna_target = x$grna_target_data_frame$grna_target[1L])
#'     o <- set_analysis_parameters(o, discovery_pairs = pairs,
#'         formula_object = ~ log(response_n_umis),
#'         control_group = "complement", side = "left",
#'         resampling_approximation = "no_approximation")
#'     o <- set_grna_assignments(o, x$grna_matrix >= 5, "example caller")
#'     o <- run_qc(o, response_n_umis_range = c(0, 1),
#'         response_n_nonzero_range = c(0, 1), p_mito_threshold = 1,
#'         n_nonzero_trt_thresh = 3L, n_nonzero_cntrl_thresh = 3L)
#'     o <- run_discovery_analysis(o, print_progress = FALSE)
#'     returned <- export_results_to_seurat(o, s)
#' }
export_results_to_seurat <- function(
    sceptre_object,
    seurat_object,
    analyses = "discovery_analysis",
    response_assay = "RNA"
) {
    require_seurat_object(seurat_object)
    object <- migrate_sceptre_object(sceptre_object)
    counts <- seurat_count_layer(seurat_object, response_assay)
    if (
        !identical(
            set_matrix_accessibility(counts),
            get_response_matrix(object)
        )
    ) {
        stop("Seurat response counts and IDs must match the imported library.")
    }
    if (!is.character(analyses) || !length(analyses) || anyNA(analyses)) {
        stop("analyses must name distinct completed analyses.")
    }
    canonical <- vapply(analyses, analysis_record_name, character(1))
    if (anyDuplicated(canonical)) {
        stop("analyses must name distinct analyses.")
    }
    records <- stats::setNames(
        lapply(canonical, function(a) get_result(object, a)),
        canonical
    )
    if ("sceptre" %in% names(SeuratObject::Misc(seurat_object))) {
        stop("Seurat object already contains sceptre miscellaneous data.")
    }
    SeuratObject::Misc(seurat_object, slot = "sceptre") <-
        owned_result_value(interop_result_metadata(object, counts, records))
    seurat_object
}

require_seurat_object <- function(object) {
    if (
        !requireNamespace("SeuratObject", quietly = TRUE) ||
            utils::packageVersion("SeuratObject") < "5.0.0"
    ) {
        stop("Install optional SeuratObject >= 5.0.0 for this interface.")
    }
    if (!methods::is(object, "Seurat")) {
        stop("Expected a Seurat object.")
    }
    methods::validObject(object)
    invisible(NULL)
}

seurat_count_layer <- function(object, assay_name) {
    if (
        !is.character(assay_name) ||
            length(assay_name) != 1L ||
            is.na(assay_name) ||
            !assay_name %in% SeuratObject::Assays(object)
    ) {
        stop("Missing count assay: ", assay_name)
    }
    assay <- object[[assay_name]]
    if (!methods::is(assay, "Assay5") && !methods::is(assay, "Assay")) {
        stop("Only in-memory Assay5 and Assay containers are supported.")
    }
    layers <- SeuratObject::Layers(assay, search = NA)
    count_layers <- grep("^counts($|[.])", layers, value = TRUE)
    if (!identical(count_layers, "counts")) {
        stop(
            "Require one unsplit counts layer; explicitly join split layers ",
            "and verify the complete library before import."
        )
    }
    x <- SeuratObject::LayerData(assay, layer = "counts")
    check_interop_counts(x, assay_name)
    check_interop_ids(rownames(x), nrow(x), "feature")
    check_interop_ids(colnames(x), ncol(x), "cell")
    if (
        !identical(rownames(x), rownames(assay)) ||
            !identical(colnames(x), colnames(assay)) ||
            !identical(colnames(x), colnames(object))
    ) {
        stop(
            "Count layer must contain all assay features and object cells ",
            "in exact source order."
        )
    }
    seurat_check_physical_names(assay, x)
    x
}

seurat_check_physical_names <- function(assay, counts) {
    raw <- if (methods::is(assay, "Assay5")) {
        SeuratObject::LayerData(assay, layer = "counts", fast = TRUE)
    } else {
        SeuratObject::LayerData(assay, layer = "counts")
    }
    for (i in seq_len(2L)) {
        stored <- dimnames(raw)[[i]]
        if (!is.null(stored) && !identical(stored, dimnames(counts)[[i]])) {
            stop("Stored count layer IDs conflict with Seurat assay maps.")
        }
    }
    invisible(NULL)
}
