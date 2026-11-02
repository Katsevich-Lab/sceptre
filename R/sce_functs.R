#' Import data from a SingleCellExperiment
#'
#' `import_data_from_sce()` imports data from a `SingleCellExperiment` holding
#' the response UMI counts in its `counts` assay and the gRNA UMI counts in the
#' `counts` assay of an alternative experiment (`altExp`), the Bioconductor
#' convention for a second feature type and the layout that
#' [as_singlecellexperiment()] returns. It is otherwise the same as
#' [import_data()].
#'
#' Unless `grna_target_data_frame` is supplied, the gRNA targets are read from
#' the alternative experiment: the gRNA IDs from its row names, the targets from
#' the `grna_target` column of its `rowData` (gRNAs whose target is `NA` are not
#' analyzed), the vectors, if any, from the `vector_id` column, and the gRNA
#' positions, if any, from its `rowRanges` (one range per gRNA, none for gRNAs
#' without a position, such as non-targeting gRNAs). Other columns of the
#' `rowData` are kept in the gRNA target data frame.
#'
#' @param sce a `SingleCellExperiment`, with responses in rows and cells in
#' columns
#' @param moi a string indicating the MOI of the dataset, either "low" or
#' "high".
#' @param grna_target_data_frame (optional) the gRNA target data frame, as in
#' [import_data()]; by default, read from the alternative experiment
#' `grna_altexp` as described above
#' @param grna_altexp the name of the alternative experiment holding the gRNA
#' counts (default `"gRNA"`)
#' @param extra_covariates (optional) the names of columns of `colData(sce)` to
#' use as extra covariates (e.g., batch); by default none. The covariates that
#' `sceptre` computes are always computed anew.
#' @param response_names (optional) a vector of human-readable response names,
#' for example `SummarizedExperiment::rowData(sce)$Symbol`; names with the
#' prefix "MT-" are taken as mitochondrial genes and are used to compute the
#' covariate `response_p_mito`. By default, the row names of `sce`.
#' @inheritParams import_data
#'
#' @return an initialized `sceptre_object`
#' @seealso [as_singlecellexperiment()] for the reverse direction
#' @export
#' @examples
#' data(lowmoi_example_data)
#' grna_experiment <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = lowmoi_example_data$grna_matrix)
#' )
#' grna_target_data_frame <- lowmoi_example_data$grna_target_data_frame
#' SummarizedExperiment::rowData(grna_experiment)$grna_target <-
#'   grna_target_data_frame$grna_target[
#'     match(rownames(grna_experiment), grna_target_data_frame$grna_id)
#'   ]
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = lowmoi_example_data$response_matrix),
#'   colData = lowmoi_example_data$extra_covariates,
#'   altExps = list(gRNA = grna_experiment)
#' )
#' sceptre_object <- import_data_from_sce(
#'   sce,
#'   moi = "low",
#'   extra_covariates = "batch"
#' )
import_data_from_sce <- function(
    sce,
    moi,
    grna_target_data_frame = NULL,
    grna_altexp = "gRNA",
    extra_covariates = character(),
    response_names = NA_character_,
    backed = FALSE,
    directory_to_write = NULL
) {
    if (!methods::is(sce, "SingleCellExperiment")) {
        stop("`sce` must be a SingleCellExperiment.")
    }
    altexp_names <- SingleCellExperiment::altExpNames(sce)
    if (
        !(is.character(grna_altexp) && length(grna_altexp) == 1L) ||
            !(grna_altexp %in% altexp_names)
    ) {
        stop(
            "`grna_altexp` must name an alternative experiment of `sce`; ",
            "`sce` has ",
            if (length(altexp_names) == 0L) {
                "none"
            } else {
                paste0("\"", altexp_names, "\"", collapse = ", ")
            },
            "."
        )
    }
    grna_experiment <- SingleCellExperiment::altExp(sce, grna_altexp)
    for (experiment in list(sce, grna_experiment)) {
        if (!("counts" %in% SummarizedExperiment::assayNames(experiment))) {
            stop(
                "`sce` and its alternative experiment `grna_altexp` must each ",
                "have a `counts` assay."
            )
        }
    }
    covariates <- SummarizedExperiment::colData(sce)
    if (
        !is.character(extra_covariates) ||
            !all(extra_covariates %in% colnames(covariates))
    ) {
        stop(
            "`extra_covariates` must name columns of `colData(sce)`; ",
            "`colData(sce)` has no column ",
            paste0(
                "`",
                setdiff(extra_covariates, colnames(covariates)),
                "`",
                collapse = ", "
            ),
            "."
        )
    }
    if (is.null(grna_target_data_frame)) {
        grna_target_data_frame <- grna_targets_from_experiment(grna_experiment)
    }
    import_data(
        response_matrix = sce_counts(sce),
        grna_matrix = sce_counts(grna_experiment),
        grna_target_data_frame = grna_target_data_frame,
        moi = moi,
        extra_covariates = if (length(extra_covariates) == 0L) {
            data.frame()
        } else {
            as.data.frame(
                covariates[, extra_covariates, drop = FALSE],
                optional = TRUE
            )
        },
        response_names = response_names,
        backed = backed,
        directory_to_write = directory_to_write
    )
}


#' Export the data of a sceptre_object as a SingleCellExperiment
#'
#' `as_singlecellexperiment()` returns the data of a `sceptre_object` as a
#' `SingleCellExperiment`, for use with other Bioconductor packages:
#' * the response UMI counts in the `counts` assay;
#' * the cell covariates (those computed by `sceptre` and the extra
#'   covariates) in `colData`;
#' * the gRNA UMI counts in the `counts` assay of the alternative experiment
#'   `"gRNA"`, with the gRNA IDs and targets in the `grna_id` and
#'   `grna_target` columns of its `rowData` (`grna_target` is `NA` for gRNAs
#'   absent from the gRNA target data frame), the other columns of the gRNA
#'   target data frame (such as `vector_id`) in its other `rowData` columns,
#'   and the gRNA positions, if any, in its `rowRanges` (one range per gRNA,
#'   none for gRNAs without a position). This `rowRanges` is a
#'   `GRangesList` that [import_data()] takes as `grna_target_data_frame`.
#'
#' Sparse counts are returned in column-compressed (`dgCMatrix`) format.
#' [import_data_from_sce()] imports the result.
#'
#' @param sceptre_object a `sceptre_object`
#'
#' @return a `SingleCellExperiment`
#' @seealso [import_data_from_sce()] for the reverse direction
#' @export
#' @examples
#' data(lowmoi_example_data)
#' sceptre_object <- import_data(
#'   response_matrix = lowmoi_example_data$response_matrix,
#'   grna_matrix = lowmoi_example_data$grna_matrix,
#'   grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'   extra_covariates = lowmoi_example_data$extra_covariates,
#'   moi = "low"
#' )
#' sce <- as_singlecellexperiment(sceptre_object)
#' # the gRNA targets
#' SummarizedExperiment::rowData(SingleCellExperiment::altExp(sce, "gRNA"))
as_singlecellexperiment <- function(sceptre_object) {
    if (!methods::is(sceptre_object, "sceptre_object")) {
        stop("`sceptre_object` must be a `sceptre_object`.")
    }
    grna_matrix <- as_column_sparse(get_grna_matrix(sceptre_object))
    grna_experiment <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = grna_matrix),
        rowRanges = grna_ranges_from_target_data_frame(
            sceptre_object@grna_target_data_frame,
            rownames(grna_matrix)
        )
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            counts = as_column_sparse(get_response_matrix(sceptre_object))
        ),
        colData = SummarizedExperiment::colData(sceptre_object@data)
    )
    SingleCellExperiment::altExp(sce, "gRNA", withDimnames = FALSE) <-
        grna_experiment
    return(sce)
}


# the counts assay of a SummarizedExperiment, a delayed (for example,
# HDF5-backed) matrix read into a sparse matrix
sce_counts <- function(experiment) {
    counts <- SummarizedExperiment::assay(experiment, "counts")
    if (methods::is(counts, "DelayedArray")) {
        counts <- methods::as(counts, "dgCMatrix")
    }
    return(counts)
}


as_column_sparse <- function(mat) {
    if (methods::is(mat, "sparseMatrix")) {
        mat <- methods::as(mat, "CsparseMatrix")
    }
    return(mat)
}


# The gRNA targets held by a gRNA experiment (see import_data_from_sce()), as a
# GRangesList (one element per gRNA, empty for gRNAs without a position) when
# the experiment holds positions, and as a data frame otherwise
grna_targets_from_experiment <- function(grna_experiment) {
    targets <- SummarizedExperiment::rowData(grna_experiment)
    if (!("grna_target" %in% colnames(targets))) {
        stop(
            "The `rowData` of the gRNA alternative experiment must have a ",
            "column `grna_target`, or `grna_target_data_frame` must be supplied."
        )
    }
    targets <- cbind(
        S4Vectors::DataFrame(grna_id = rownames(grna_experiment)),
        targets[, setdiff(colnames(targets), "grna_id"), drop = FALSE]
    )
    keep <- !is.na(targets$grna_target)
    positions <- if (
        methods::is(grna_experiment, "RangedSummarizedExperiment")
    ) {
        SummarizedExperiment::rowRanges(grna_experiment)
    } else {
        NULL
    }
    if (methods::is(positions, "GRanges")) {
        positions <- methods::as(positions, "GRangesList")
    }
    if (is.null(positions) || all(lengths(positions) == 0L)) {
        return(as.data.frame(targets[keep, , drop = FALSE], optional = TRUE))
    }
    positions <- positions[keep]
    S4Vectors::mcols(positions) <- targets[keep, , drop = FALSE]
    names(positions) <- NULL
    return(positions)
}


# The positions of the gRNAs `grna_ids` given by the columns `chr`, `start` and
# `end` of the gRNA target data frame, as a GRangesList (one element per gRNA,
# named by gRNA ID, empty for gRNAs without a position), with the gRNA IDs and
# the other columns as metadata columns
grna_ranges_from_target_data_frame <- function(
    grna_target_data_frame,
    grna_ids
) {
    if (anyDuplicated(grna_target_data_frame$grna_id)) {
        stop(
            "The gRNA target data frame lists some gRNAs more than once; ",
            "the gRNA targets cannot be stored one row per gRNA."
        )
    }
    targets <- grna_target_data_frame[
        match(grna_ids, grna_target_data_frame$grna_id),
        setdiff(
            colnames(grna_target_data_frame),
            c("grna_id", "chr", "start", "end")
        ),
        drop = FALSE
    ]
    rownames(targets) <- NULL
    targets <- cbind(data.frame(grna_id = grna_ids), targets)
    chr <- character(length(grna_ids))
    start <- end <- integer(length(grna_ids))
    has_position <- rep(FALSE, length(grna_ids))
    if (all(c("chr", "start", "end") %in% colnames(grna_target_data_frame))) {
        idx <- match(grna_ids, grna_target_data_frame$grna_id)
        chr <- as.character(grna_target_data_frame$chr[idx])
        start <- grna_target_data_frame$start[idx]
        end <- grna_target_data_frame$end[idx]
        has_position <- !is.na(chr) & !is.na(start) & !is.na(end)
    }
    ranges <- GenomicRanges::GRanges(
        seqnames = chr[has_position],
        ranges = IRanges::IRanges(
            start = start[has_position],
            end = end[has_position]
        )
    )
    positions <- S4Vectors::splitAsList(
        ranges,
        factor(which(has_position), levels = seq_along(grna_ids))
    )
    names(positions) <- grna_ids
    S4Vectors::mcols(positions) <- S4Vectors::DataFrame(
        targets,
        check.names = FALSE
    )
    return(positions)
}
