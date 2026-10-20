#' Import raw counts from a SingleCellExperiment
#'
#' Import the complete response library and one guide alternative experiment.
#' No normalization, feature selection, or automatic name matching is performed.
#' Cell covariates are computed from all supplied counts before any SCEPTRE QC.
#' Supply the original, unfiltered feature library: counts omitted upstream
#' cannot be recovered. Delayed or disk-backed assays must be converted
#' explicitly by the caller after checking their memory requirements.
#'
#' @param sce A SingleCellExperiment containing raw response counts.
#' @param grna_target_data_frame Complete guide metadata in guide assay row
#'   order.
#'   With vector_id, assignments subsequently refer to the collapsed vector IDs.
#' @param moi Either `"low"` or `"high"`.
#' @param response_assay Raw response assay name; default `"counts"`.
#' @param grna_assay Raw guide assay name; default `"counts"`.
#' @param grna_alt_exp Name of the guide alternative experiment.
#' @param covariates Names of colData columns to include, in source cell order.
#' @param response_name_column Optional rowData column of response names used
#'   to identify mitochondrial genes. Otherwise use response row names.
#' @return A new sceptre_object. The input SCE is unchanged.
#' @export
#' @examples
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'     data("lowmoi_example_data")
#'     x <- lowmoi_example_data
#'     sce <- SingleCellExperiment::SingleCellExperiment(
#'         assays = list(counts = x$response_matrix))
#'     SingleCellExperiment::altExp(sce, "gRNA") <-
#'         SummarizedExperiment::SummarizedExperiment(
#'             assays = list(counts = x$grna_matrix))
#'     object <- import_from_sce(sce, x$grna_target_data_frame, "low")
#' }
import_from_sce <- function(
    sce,
    grna_target_data_frame,
    moi,
    response_assay = "counts",
    grna_alt_exp = "gRNA",
    grna_assay = "counts",
    covariates = character(),
    response_name_column = NULL
) {
    require_sce(sce)
    response <- sce_count_assay(sce, response_assay, "response")
    if (!grna_alt_exp %in% SingleCellExperiment::altExpNames(sce)) {
        stop("Guide alternative experiment is missing: ", grna_alt_exp)
    }
    guides <- SingleCellExperiment::altExp(
        sce,
        grna_alt_exp,
        withDimnames = FALSE
    )
    grna <- sce_count_assay(guides, grna_assay, "guide")
    if (!identical(colnames(response), colnames(grna))) {
        stop("Response and guide cell IDs must be identical and in order.")
    }
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
    import_data(
        response,
        grna,
        grna_target_data_frame,
        moi,
        extra_covariates = sce_extra_covariates(sce, covariates),
        response_names = sce_response_names(sce, response_name_column)
    )
}

sce_extra_covariates <- function(sce, covariates) {
    cd <- SummarizedExperiment::colData(sce)
    reserved <- c(
        "response_n_umis",
        "response_n_nonzero",
        "response_p_mito",
        "grna_n_umis",
        "grna_n_nonzero",
        "grna_frac_umis_max_feature",
        "grna_feature_w_max_expression"
    )
    if (
        !is.character(covariates) ||
            anyNA(covariates) ||
            anyDuplicated(covariates) ||
            !all(covariates %in% colnames(cd))
    ) {
        stop("covariates must name distinct existing colData columns.")
    }
    if (any(covariates %in% reserved)) {
        stop(
            "Selected colData covariates conflict with computed count ",
            "summaries."
        )
    }
    extra <- as.data.frame(cd[, covariates, drop = FALSE])
    rownames(extra) <- colnames(sce)
    extra
}

sce_response_names <- function(sce, response_name_column) {
    response_names <- NA_character_
    if (!is.null(response_name_column)) {
        rd <- SummarizedExperiment::rowData(sce)
        if (
            length(response_name_column) != 1L ||
                !response_name_column %in% colnames(rd)
        ) {
            stop("response_name_column must name an existing rowData column.")
        }
        response_names <- as.character(rd[[response_name_column]])
        if (anyNA(response_names) || any(!nzchar(response_names))) {
            stop("Response names must be nonmissing and nonempty.")
        }
    }
    response_names
}

#' Supply externally computed guide assignments
#'
#' Use after set_analysis_parameters(), instead of assign_grnas(). Assignments
#' refer to every original cell and every canonical guide (vector IDs after
#' vector collapsing). All-zero columns explicitly represent unassigned cells.
#' Low-MOI zero/multiple-guide exclusions and target grouping follow the usual
#' SCEPTRE rules. Changing assignments invalidates QC, results and affected
#' fits.
#'
#' @param sceptre_object A configured sceptre_object.
#' @param assignments Named binary numeric/logical guide-by-cell matrix.
#'   Row IDs must equal canonical guide metadata order; column IDs must equal
#'   original count matrix cell order. No reordering or positional fallback.
#' @param provenance Nonempty description of the assignment method/version.
#' @return An updated sceptre_object; caller-owned inputs remain unchanged.
#' @export
#' @examples
#' data("lowmoi_example_data")
#' x <- lowmoi_example_data
#' object <- import_data(x$response_matrix, x$grna_matrix,
#'                       x$grna_target_data_frame, "low")
#' object <- set_analysis_parameters(object)
#' object <- set_grna_assignments(object, x$grna_matrix >= 5,
#'                                "demonstration caller: UMI >=5")
set_grna_assignments <- function(sceptre_object, assignments, provenance) {
    object <- migrate_sceptre_object(sceptre_object)
    require_stage(object, "design")
    guides <- object@data@grna_target_data_frame$grna_id
    cells <- colnames(get_grna_matrix(object))
    check_interop_ids(guides, length(guides), "canonical guide")
    check_interop_ids(cells, ncol(get_grna_matrix(object)), "source cell")
    check_binary_assignments(assignments, guides, cells)
    if (
        !is.character(provenance) ||
            length(provenance) != 1L ||
            is.na(provenance) ||
            !nzchar(trimws(provenance))
    ) {
        stop("provenance must describe the external assignment method/version.")
    }
    params <- list(
        provenance = provenance,
        assignments_sha256 = digest::digest(assignments, algo = "sha256")
    )
    object <- begin_preparation_stage(object, "assignment")
    object <- invalidate_fit_caches(
        object,
        "assignments",
        response_cache_assignments_changed(object, "external", params)
    )
    object@spec@grna_assignment_method <- "external"
    object@spec@grna_assignment_hyperparameters <- params
    object@preparation@assignment_exclusions <- integer()
    object@preparation@guide_assignments <- stats::setNames(
        lapply(seq_along(guides), function(i) {
            as.integer(which(assignments[i, ] != 0))
        }),
        guides
    )
    object <- process_initial_assignment_list(object)
    if (object@data@low_moi) {
        object@preparation@assignment_exclusions <- as.integer(which(
            object@preparation@guide_counts_per_cell != 1L
        ))
    }
    commit_preparation_stage(object, "assignment")
}

#' Export SCEPTRE results into an original SingleCellExperiment
#'
#' Return a copy with metadata(sce)$sceptre containing completed pair-level
#' result records, original-cell covariates, retained cell IDs, and assignment
#' provenance. Pair tables remain tables: many targets can map to one response.
#' Counts, rowData, colData and alternative experiments are left unchanged.
#' Existing sceptre metadata is rejected to avoid accidental overwrite.
#'
#' @param sceptre_object A sceptre_object with completed requested analyses.
#' @param sce The original SCE with exactly the imported response counts/IDs.
#' @param analyses Analysis names accepted by get_result(), or short aliases
#'   such as `"discovery_analysis"`. Default discovery.
#' @param response_assay The original raw response assay name.
#' @return A new SingleCellExperiment with owned result records in metadata.
#' @export
#' @examples
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'     data("lowmoi_example_data")
#'     x <- lowmoi_example_data
#'     sce <- SingleCellExperiment::SingleCellExperiment(
#'         assays = list(counts = x$response_matrix))
#'     SingleCellExperiment::altExp(sce, "gRNA") <-
#'         SummarizedExperiment::SummarizedExperiment(
#'             assays = list(counts = x$grna_matrix))
#'     object <- import_from_sce(sce, x$grna_target_data_frame, "low")
#'     pairs <- data.frame(
#'         response_id = rownames(sce)[1L],
#'         grna_target = x$grna_target_data_frame$grna_target[1L])
#'     object <- set_analysis_parameters(object, discovery_pairs = pairs,
#'         formula_object = ~ log(response_n_umis),
#'         control_group = "complement",
#'         side = "left", resampling_approximation = "no_approximation")
#'     object <- set_grna_assignments(object, x$grna_matrix >= 5,
#'                                    "demonstration caller: UMI >=5")
#'     object <- run_qc(object, response_n_umis_range = c(0, 1),
#'         response_n_nonzero_range = c(0, 1), p_mito_threshold = 1,
#'         n_nonzero_trt_thresh = 3L, n_nonzero_cntrl_thresh = 3L)
#'     object <- run_discovery_analysis(object, print_progress = FALSE)
#'     returned <- export_results_to_sce(object, sce)
#' }
export_results_to_sce <- function(
    sceptre_object,
    sce,
    analyses = "discovery_analysis",
    response_assay = "counts"
) {
    require_sce(sce)
    object <- migrate_sceptre_object(sceptre_object)
    counts <- sce_count_assay(sce, response_assay, "response")
    if (
        !identical(
            set_matrix_accessibility(counts, make_row_accessible = TRUE),
            get_response_matrix(object)
        )
    ) {
        stop(
            "SCE response counts and IDs must exactly match the imported ",
            "library."
        )
    }
    if (
        !is.character(analyses) ||
            !length(analyses) ||
            anyNA(analyses) ||
            anyDuplicated(analyses)
    ) {
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
    md <- S4Vectors::metadata(sce)
    if ("sceptre" %in% names(md)) {
        stop("SCE already contains sceptre metadata.")
    }
    md$sceptre <- owned_result_value(sce_result_metadata(
        object,
        counts,
        records
    ))
    S4Vectors::metadata(sce) <- md
    sce
}

sce_result_metadata <- function(object, counts, records) {
    list(
        package_version = as.character(utils::packageVersion("sceptre")),
        schema_version = object@schema_version,
        import_identity = object@data@identity,
        response_ids = rownames(counts),
        cell_ids = colnames(counts),
        retained_cell_ids = colnames(counts)[object@preparation@analysis_cells],
        cell_covariates = get_cell_covariates(object),
        assignment_method = object@spec@grna_assignment_method,
        assignment_provenance = object@spec@grna_assignment_hyperparameters,
        result_provenance = stats::setNames(
            lapply(names(records), function(a) {
                object@results@records[[a]]$provenance
            }),
            names(records)
        ),
        results = records
    )
}

# Internal boundary checks never realize delayed assays or infer missing names.
require_sce <- function(sce) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop(
            "Install the optional SingleCellExperiment package to use ",
            "this interface."
        )
    }
    if (!methods::is(sce, "SingleCellExperiment")) {
        stop("sce must be a SingleCellExperiment.")
    }
    methods::validObject(sce)
    invisible(NULL)
}

check_interop_ids <- function(ids, n, label) {
    if (
        is.null(ids) ||
            length(ids) != n ||
            !n ||
            anyNA(ids) ||
            any(!nzchar(ids)) ||
            anyDuplicated(ids)
    ) {
        stop(label, " IDs must be present, unique, nonmissing and nonempty.")
    }
}

sce_count_assay <- function(sce, assay_name, label) {
    if (
        !is.character(assay_name) ||
            length(assay_name) != 1L ||
            is.na(assay_name) ||
            !assay_name %in% SummarizedExperiment::assayNames(sce)
    ) {
        stop(label, " raw count assay is missing.")
    }
    sce_check_assay_dimnames(sce, assay_name)
    x <- SummarizedExperiment::assay(sce, assay_name)
    if (
        !(is.matrix(x) && is.numeric(x)) &&
            !any(vapply(
                c("dgCMatrix", "dgRMatrix", "dgTMatrix"),
                function(cl) methods::is(x, cl),
                logical(1)
            ))
    ) {
        stop(
            label,
            " counts must be an in-memory numeric matrix or dgC/dgR/dgTMatrix."
        )
    }
    # Canonicalization sums duplicate sparse coordinates before validation.
    values <- if (is.matrix(x)) x else methods::as(x, "RsparseMatrix")@x
    if (any(!is.finite(values)) || any(values < 0 | values != floor(values))) {
        stop(
            label,
            " counts must be finite nonnegative integers; supply raw counts."
        )
    }
    x
}

sce_check_assay_dimnames <- function(sce, assay_name) {
    raw <- SummarizedExperiment::assay(sce, assay_name, withDimnames = FALSE)
    labels <- c("feature", "cell")
    for (i in 1:2) {
        check_interop_ids(dimnames(sce)[[i]], dim(sce)[i], labels[i])
        stored <- dimnames(raw)[[i]]
        if (!is.null(stored) && !identical(stored, dimnames(sce)[[i]])) {
            stop("Stored assay ", labels[i], " IDs conflict with SCE IDs.")
        }
    }
    invisible(NULL)
}

check_binary_assignments <- function(x, guides, cells) {
    if (!is.matrix(x) && !methods::is(x, "sparseMatrix")) {
        stop("assignments must be a binary matrix.")
    }
    if (!identical(rownames(x), guides) || !identical(colnames(x), cells)) {
        stop("Assignment guide and cell IDs must exactly match source order.")
    }
    # Also interprets sparse pattern matrices as TRUE membership.
    z <- if (is.matrix(x)) {
        x
    } else {
        methods::as(methods::as(x, "dMatrix"), "RsparseMatrix")@x
    }
    if (!(is.numeric(z) || is.logical(z)) || anyNA(z) || any(!z %in% c(0, 1))) {
        stop(
            "assignments must contain only nonmissing binary 0/1 or ",
            "FALSE/TRUE values."
        )
    }
}
