#' Write or read a sceptre_object
#'
#' `write_sceptre_object()` writes a `sceptre_object` to a directory, and
#' `read_sceptre_object()` reads it back. The directory holds two files:
#' `screen.h5mu`, a MuData file holding the data (the response and gRNA UMI
#' counts and the cell covariates), which Bioconductor's `MuData` and Python's
#' `mudata` packages can read, and `sceptre_object.rds`, which holds the
#' parameters, the intermediate values and the results of the analysis. The
#' directory has the same layout for a standard and a disk-backed
#' `sceptre_object`; a disk-backed `sceptre_object` whose response counts are
#' already in `screen.h5mu` in the directory (for example, one imported with
#' `backed = TRUE` and `directory_to_write` set to the directory) is written
#' without copying them.
#'
#' `read_sceptre_object()` reads the response counts into memory
#' (`backed = FALSE`) or leaves them on disk (`backed = TRUE`), in which case
#' the directory must stay in place while the `sceptre_object` is used. Either
#' way, a `sceptre_object` written by `write_sceptre_object()` can be read as
#' a standard or a disk-backed `sceptre_object`.
#'
#' @param sceptre_object a `sceptre_object`
#' @param directory a file path to a directory
#' @param backed (default `FALSE`) a logical value indicating whether to leave
#' the response counts on disk (`TRUE`) or to read them into memory (`FALSE`)
#'
#' @return `write_sceptre_object()` returns `NULL` invisibly, and
#' `read_sceptre_object()` returns a `sceptre_object`.
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
#' directory <- file.path(tempdir(), "sceptre_object")
#' write_sceptre_object(sceptre_object, directory)
#' list.files(directory)
#' sceptre_object <- read_sceptre_object(directory)
#' sceptre_object_backed <- read_sceptre_object(directory, backed = TRUE)
write_sceptre_object <- function(sceptre_object, directory) {
    if (!methods::is(sceptre_object, "sceptre_object")) {
        stop("`sceptre_object` must be a `sceptre_object`.")
    }
    directory <- prepare_directory_to_write(directory)
    file <- file.path(directory, "screen.h5mu")
    response_matrix <- get_response_matrix(sceptre_object)
    if (
        !(normalizePath(file, winslash = "/", mustWork = FALSE) %in%
            backing_files(response_matrix))
    ) {
        write_screen(sceptre_object@data, file)
    }
    # the counts are in screen.h5mu; keep only their dimensions and names
    sceptre_object <- set_response_matrix(
        sceptre_object,
        empty_counts(response_matrix)
    )
    sceptre_object <- set_grna_matrix(
        sceptre_object,
        empty_counts(get_grna_matrix(sceptre_object))
    )
    saveRDS(sceptre_object, file.path(directory, "sceptre_object.rds"))
    invisible(NULL)
}


#' @rdname write_sceptre_object
#' @export
read_sceptre_object <- function(directory, backed = FALSE) {
    if (!(isTRUE(backed) || isFALSE(backed))) {
        stop("`backed` must be TRUE or FALSE.")
    }
    directory <- prepare_file_path(directory, "`directory`")
    files <- file.path(directory, c("sceptre_object.rds", "screen.h5mu"))
    if (!all(file.exists(files))) {
        stop(
            "The directory ",
            directory,
            " must contain the files sceptre_object.rds and screen.h5mu ",
            "written by `write_sceptre_object()`."
        )
    }
    sceptre_object <- readRDS(files[1L])
    if (!methods::is(sceptre_object, "sceptre_object")) {
        stop(files[1L], " does not contain a `sceptre_object`.")
    }
    screen <- read_screen(files[2L], backed = backed)
    for (modality in c("response", "grna")) {
        stored <- if (modality == "response") {
            get_response_matrix(sceptre_object)
        } else {
            get_grna_matrix(sceptre_object)
        }
        counts <- screen[[modality]]
        # screen.h5mu names the cells by their indices when the sceptre_object
        # has no cell names (see write_screen())
        cell_ids <- colnames(stored)
        if (is.null(cell_ids)) {
            cell_ids <- as.character(seq_len(ncol(stored)))
        }
        if (
            !identical(dim(counts), dim(stored)) ||
                !identical(rownames(counts), rownames(stored)) ||
                !identical(as.character(colnames(counts)), cell_ids)
        ) {
            stop(
                files[2L],
                " does not hold the data of the `sceptre_object` in ",
                files[1L],
                "."
            )
        }
        dimnames(counts) <- dimnames(stored)
        sceptre_object <- if (modality == "response") {
            set_response_matrix(sceptre_object, counts)
        } else {
            set_grna_matrix(sceptre_object, counts)
        }
    }
    methods::validObject(sceptre_object)
    return(sceptre_object)
}


# --- screen.h5mu ---------------------------------------------------------------

# The format stamp of screen.h5mu; a reader refuses a stamp of a later format.
sceptre_format_stamp <- "{\"format\":3}"

# Writes the data of a sceptre_object (its `data` slot) to `file`, a MuData
# file with two modalities, "response" and "grna", each holding its counts as
# X and the feature names as var names; the cell covariates are the obs of
# the response modality. The counts are written feature-contiguous (AnnData
# csc_matrix; see backed_counts.R): MuData writes the gRNA counts in that
# layout from a dgRMatrix, and the response counts, in memory or on disk, are
# copied into the file block by block (write_backed_counts()) after MuData has
# written everything else around an empty placeholder; for in-memory counts
# this is several times faster than MuData's writer.
#
# MuData 1.16 reads NA in a text column back as an empty string (a
# categorical column keeps it), so text covariates are written as categoricals
# and listed in the metadata entry `covariate_text_columns`; h5mu stores a
# missing number as NaN, the AnnData convention. MultiAssayExperiment needs
# cell names, so a screen without them is written with the cell indices as
# names and the metadata entry `cell_names_absent`. The file is written next
# to its destination and moved there once complete. `response_source`, if
# given, holds the response counts to write in place of the counts of
# `screen`, which are then an empty matrix with their dimensions and names
# (see import_data_backed()).
write_screen <- function(screen, file, response_source = NULL) {
    counts <- SummarizedExperiment::assay(screen, "counts", withDimnames = FALSE)
    if (
        is_backed(counts) &&
            normalizePath(file, winslash = "/", mustWork = FALSE) %in%
                backing_files(counts)
    ) {
        stop(
            "The counts are read from ",
            file,
            ", which cannot be overwritten."
        )
    }
    meta <- list(sceptre_format = sceptre_format_stamp)
    cell_ids <- colnames(screen)
    if (is.null(cell_ids)) {
        meta$cell_names_absent <- TRUE
        cell_ids <- as.character(seq_len(ncol(screen)))
    } else if (anyDuplicated(cell_ids)) {
        stop(
            "The cell IDs must be unique to write screen.h5mu; ",
            "some are duplicated, for example ",
            cell_ids[anyDuplicated(cell_ids)],
            "."
        )
    }
    covariates <- as.data.frame(
        SummarizedExperiment::colData(screen),
        optional = TRUE
    )
    text <- names(covariates)[vapply(covariates, is.character, logical(1))]
    covariates[text] <- lapply(covariates[text], factor)
    meta$covariate_text_columns <- text
    rownames(covariates) <- cell_ids
    response <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = empty_counts(counts)),
        colData = S4Vectors::DataFrame(covariates, check.names = FALSE)
    )
    grna <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            counts = as_row_counts(SummarizedExperiment::assay(
                SingleCellExperiment::altExp(
                    screen,
                    "gRNA",
                    withDimnames = FALSE
                ),
                "counts",
                withDimnames = FALSE
            ))
        )
    )
    colnames(response) <- colnames(grna) <- cell_ids
    temporary_file <- tempfile(
        pattern = "screen",
        tmpdir = dirname(file),
        fileext = ".h5mu"
    )
    on.exit(unlink(temporary_file), add = TRUE)
    MuData::writeH5MU(
        MultiAssayExperiment::MultiAssayExperiment(
            list(response = response, grna = grna),
            metadata = meta
        ),
        temporary_file
    )
    write_h5mu_mapping_encodings(temporary_file)
    write_nullable_boolean_arrays(temporary_file)
    rhdf5::h5delete(temporary_file, response_counts_group)
    write_backed_counts(
        if (is.null(response_source)) counts else response_source,
        temporary_file,
        response_counts_group
    )
    if (!file.rename(temporary_file, file)) {
        stop("Could not write ", file, ".")
    }
    invisible(file)
}


# Reads the counts in screen.h5mu: list(response, grna), features x cells.
# With `backed = TRUE` the response counts stay on disk, as a DelayedMatrix
# read through HDF5Array; the gRNA counts are read into memory either way.
read_screen <- function(file, backed = FALSE) {
    # MultiAssayExperiment warns that it rebuilds the sample map read from the
    # file; that is expected
    mae <- withCallingHandlers(
        MuData::readH5MU(file, backed = backed),
        warning = function(w) {
            if (grepl("sampleMap", conditionMessage(w), fixed = TRUE)) {
                invokeRestart("muffleWarning")
            }
        }
    )
    check_screen_format(S4Vectors::metadata(mae)$sceptre_format, file)
    response <- SummarizedExperiment::assay(mae[["response"]], 1L)
    if (!backed) {
        response <- as_row_counts(response)
    }
    grna <- as_row_counts(SummarizedExperiment::assay(mae[["grna"]], 1L))
    list(response = response, grna = grna)
}


check_screen_format <- function(stamp, file) {
    stamp <- as.vector(stamp)
    pattern <- "^\\{\"format\":([0-9]+)[,}].*$"
    if (
        !is.character(stamp) || length(stamp) != 1L || !grepl(pattern, stamp)
    ) {
        stop(file, " is not a screen.h5mu file written by sceptre.")
    }
    format <- as.integer(sub(pattern, "\\1", stamp))
    if (format > 3L) {
        stop(
            file,
            " was written by a newer version of sceptre; upgrade sceptre to ",
            "read it."
        )
    }
    invisible(NULL)
}
