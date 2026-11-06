#' Import data from a SingleCellExperiment
#'
#' `import_data_from_sce()` imports data from a `SingleCellExperiment` holding
#' the response UMI counts in its `counts` assay and the gRNA UMI counts in the
#' `counts` assay of an alternative experiment (`altExp`), the Bioconductor
#' convention for a second feature type. It is otherwise the same as
#' [import_data()].
#'
#' Unless `grna_target_data_frame` is supplied, the gRNA targets are read from
#' the alternative experiment: the gRNA IDs from its row names, the targets from
#' the `grna_target` column of its `rowData` (gRNAs whose target is `NA` are not
#' analyzed), the vectors, if any, from the `vector_id` column, and the
#' chromosome, start and end of the region that each gRNA targets, if any, from
#' the `target_chr`, `target_start` and `target_end` columns (`NA` for gRNAs
#' without a target region, such as non-targeting gRNAs), which become the
#' columns `chr`, `start` and `end` of the gRNA target data frame. Other columns
#' of the `rowData` are kept in the gRNA target data frame. The `rowRanges` of
#' the alternative experiment, which would give the positions of the gRNAs
#' themselves, are not used.
#'
#' The `rowRanges` of `sce`, the positions of the responses (for example, of the
#' genes), are not used; [construct_cis_pairs()] takes the position of each
#' response, its transcription start site, as a data frame
#' (`response_position_data_frame`).
#'
#' A `SingleCellExperiment` made by another package may name the alternative
#' experiment of the gRNA counts differently; for example, one made from Cell
#' Ranger output with `DropletUtils::read10xCounts()` and
#' `SingleCellExperiment::splitAltExps()` names it after the 10x feature type,
#' `"CRISPR Guide Capture"`. Set `grna_altexp` to that name.
#'
#' The counts of `sce` are read into memory, also with `backed = TRUE`, in
#' which case they are then written to `screen.h5mu` in `directory_to_write`;
#' a `SingleCellExperiment` too large for memory cannot be imported.
#'
#' @param sce a `SingleCellExperiment`, with responses in rows and cells in
#' columns
#' @param moi a string indicating the MOI of the dataset, either "low" or
#' "high".
#' @param grna_target_data_frame (optional) the gRNA target data frame, as in
#' [import_data()]; by default, read from the alternative experiment
#' `grna_altexp` as described above
#' @param grna_altexp the name of the alternative experiment holding the gRNA
#' counts (default `"grna"`)
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
#'   altExps = list(grna = grna_experiment)
#' )
#' sceptre_object <- import_data_from_sce(
#'   sce,
#'   moi = "low",
#'   extra_covariates = "batch"
#' )
#'
#' # the gRNA counts in an alternative experiment named after the 10x feature
#' # type, as in a SingleCellExperiment made from Cell Ranger output
#' SingleCellExperiment::altExpNames(sce) <- "CRISPR Guide Capture"
#' sceptre_object <- import_data_from_sce(
#'   sce,
#'   moi = "low",
#'   grna_altexp = "CRISPR Guide Capture",
#'   extra_covariates = "batch"
#' )
import_data_from_sce <- function(
    sce,
    moi,
    grna_target_data_frame = NULL,
    grna_altexp = "grna",
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
    sceptre_object <- import_data(
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
    return(sceptre_object)
}


# the counts assay of a SummarizedExperiment; counts in any other class than
# a base or Matrix matrix (for example, a SparseArray or a delayed,
# HDF5-backed matrix) are read into a dgCMatrix
sce_counts <- function(experiment) {
    counts <- SummarizedExperiment::assay(experiment, "counts")
    if (!(is.matrix(counts) || methods::is(counts, "Matrix"))) {
        counts <- methods::as(counts, "dgCMatrix")
    }
    return(counts)
}


# The columns of the gRNA target data frame that give the region each gRNA
# targets, named as in the rowData of a gRNA experiment (see
# import_data_from_sce())
target_region_columns <- c(
    chr = "target_chr",
    start = "target_start",
    end = "target_end"
)


# The gRNA target data frame held by a gRNA experiment (see
# import_data_from_sce()): its rowData, with the gRNA IDs from its row names and
# the target region columns renamed
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
    present <- target_region_columns %in% colnames(targets)
    if (any(present)) {
        if (any(names(target_region_columns) %in% colnames(targets))) {
            stop(
                "The `rowData` of the gRNA alternative experiment must give ",
                "the target regions in the columns `target_chr`, ",
                "`target_start` and `target_end` or in `chr`, `start` and ",
                "`end`, not both."
            )
        }
        idx <- match(target_region_columns[present], colnames(targets))
        colnames(targets)[idx] <- names(target_region_columns)[present]
    }
    keep <- !is.na(targets$grna_target)
    return(as.data.frame(targets[keep, , drop = FALSE], optional = TRUE))
}
