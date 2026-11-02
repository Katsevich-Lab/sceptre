############
# R MATRICES
############
#' Import data
#'
#' `import_data()` imports data from a collection of R objects to create a
#' `sceptre_object`. Users can create either a standard `sceptre` object or a
#' disk-backed `sceptre` object, whose response counts stay on disk; the latter
#' is more appropriate for large-scale data. See
#' \href{https://timothy-barry.github.io/sceptre-book/import-data.html#import-data-from-a-collection-of-r-objects}{Chapter 1 of the manual}
#' for more detailed information about this function.
#'
#' @param response_matrix a matrix of response UMI counts, with responses in
#' rows and cells in columns. The matrix should be of type `"matrix"`,
#' `"dgCMatrix"`, `"dgRMatrix"`, or `"dgTMatrix"`. The row names of the matrix
#' should give the response IDs. The column names, when supplied, should give
#' the cell IDs. The entries must be integer counts; a matrix with non-integer
#' entries (for example, normalized expression values) is an error.
#' @param grna_matrix a matrix of gRNA UMI counts, with gRNAs in rows and cells
#' in columns. The matrix should be of type `"matrix"`, `"dgCMatrix"`,
#' `"dgRMatrix"`, or `"dgTMatrix"`. The row names of the matrix should give the
#' gRNA IDs. The column names, when supplied, should give the cell IDs. The
#' entries must be integer counts; a matrix with non-integer entries is an
#' error.
#' @param grna_target_data_frame a data frame containing columns `grna_id` and
#' `grna_target` mapping each individual gRNA to its target. Non-targeting gRNAs
#' should be assigned a label of "non-targeting". Optionally,
#' `grna_target_data_frame` can contain columns `chr`, `start`, and `end`,
#' giving the chromosome, start coordinate, and end coordiante, respectively, of
#' each gRNA. Additionally, `grna_target_data_frame` can contain the column
#' `vector_id` specifying the vector to which a given gRNA belongs.
#' `grna_target_data_frame` can instead be a `GRanges` object, one range per
#' gRNA, or a `GRangesList` object, one element per gRNA, empty for gRNAs
#' without a position (such as non-targeting gRNAs), as in the `rowRanges` of
#' a `SummarizedExperiment`; the ranges give the positions (`chr`, `start`, and
#' `end`), and the metadata columns give `grna_id`, `grna_target`, and
#' optionally `vector_id`.
#' @param moi a string indicating the MOI of the dataset, either "low" or
#' "high".
#' @param extra_covariates (optional) a data frame containing extra covariates
#' (e.g., batch, biological replicate) beyond those that `sceptre` can compute.
#' The rows should follow the same cell ordering as the columns of
#' `response_matrix` and `grna_matrix`. Cell IDs may be supplied as row names.
#' When cell IDs are supplied for at least two of these three objects, they must
#' be identical.
#' @param response_names (optional) a vector of human-readable response names;
#' names with the prefix "MT-" are taken as mitochondrial genes and are used to
#' compute the covariate `response_p_mito`.
#' @param backed (default `FALSE`) a logical value (i.e., `TRUE` or `FALSE`)
#' indicating whether to create a disk-backed `sceptre_object` (`TRUE`), whose
#' response counts stay on disk in the file `screen.h5mu` in
#' `directory_to_write` and are read as they are needed, or a standard
#' `sceptre_object` (`FALSE`), which holds them in memory. The gRNA counts and
#' the cell covariates are held in memory either way. The analysis of a
#' disk-backed `sceptre_object` uses permutations (`resampling_mechanism =
#' "permutations"` in [set_analysis_parameters()]).
#' @param directory_to_write (optional) a file path to a directory in which to
#' write `screen.h5mu`, a MuData file holding the response and gRNA counts and
#' the cell covariates. Must be supplied if `backed` is set to `TRUE`. Save
#' the `sceptre_object` in the same directory with [write_sceptre_object()].
#' @param use_ondisc deprecated; use `backed`.
#'
#' @return an initialized `sceptre_object`
#' @export
#' @examples
#' data("lowmoi_example_data")
#' # 1. initialize a standard sceptre_object from R objects
#' sceptre_object <- import_data(
#'   response_matrix = lowmoi_example_data$response_matrix,
#'   grna_matrix = lowmoi_example_data$grna_matrix,
#'   grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'   extra_covariates = lowmoi_example_data$extra_covariates,
#'   moi = "low"
#' )
#'
#' # 2. initialize a disk-backed sceptre_object from R objects
#' sceptre_object <- import_data(
#'   response_matrix = lowmoi_example_data$response_matrix,
#'   grna_matrix = lowmoi_example_data$grna_matrix,
#'   grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'   extra_covariates = lowmoi_example_data$extra_covariates,
#'   moi = "low",
#'   backed = TRUE,
#'   directory_to_write = tempdir()
#' )
import_data <- function(
    response_matrix,
    grna_matrix,
    grna_target_data_frame,
    moi,
    extra_covariates = data.frame(),
    response_names = NA_character_,
    backed = FALSE,
    directory_to_write = NULL,
    use_ondisc = NULL
) {
    backed <- resolve_backed(backed, use_ondisc)
    if (backed) {
        directory_to_write <- prepare_directory_to_write(directory_to_write)
    }
    if (
        methods::is(grna_target_data_frame, "GRanges") ||
            methods::is(grna_target_data_frame, "GRangesList")
    ) {
        grna_target_data_frame <- grna_target_data_frame_from_granges(
            grna_target_data_frame
        )
    }
    sceptre_object <- import_data_memory(
        response_matrix,
        grna_matrix,
        grna_target_data_frame,
        moi,
        extra_covariates,
        response_names
    )
    if (backed) {
        sceptre_object <- back_response_counts(
            sceptre_object,
            directory_to_write
        )
    }
    return(sceptre_object)
}


# `use_ondisc`, the argument of 0.99.1 and earlier, is deprecated in favor of
# `backed`
resolve_backed <- function(backed, use_ondisc) {
    if (!is.null(use_ondisc)) {
        warning(
            "The argument `use_ondisc` is deprecated; use `backed` instead.",
            call. = FALSE
        )
        backed <- use_ondisc
    }
    if (!(isTRUE(backed) || isFALSE(backed))) {
        stop("`backed` must be TRUE or FALSE.")
    }
    return(backed)
}


# Writes the data of an in-memory sceptre_object to `directory/screen.h5mu`
# and keeps its response counts there, on disk.
back_response_counts <- function(sceptre_object, directory) {
    file <- file.path(directory, "screen.h5mu")
    write_screen(sceptre_object@data, file)
    response_matrix <- get_response_matrix(sceptre_object)
    backed_counts <- read_screen(file, backed = TRUE)$response
    dimnames(backed_counts) <- dimnames(response_matrix)
    set_response_matrix(sceptre_object, backed_counts)
}

# convert a GRanges or GRangesList of gRNAs, whose metadata columns hold
# `grna_id`, `grna_target` (and optionally `vector_id`), to a gRNA target data
# frame; an empty GRangesList element is a gRNA without a position
grna_target_data_frame_from_granges <- function(grna_ranges) {
    grna_target_data_frame <- GenomicRanges::as.data.frame(
        GenomicRanges::mcols(grna_ranges)
    )
    if (
        !all(c("grna_id", "grna_target") %in% colnames(grna_target_data_frame))
    ) {
        stop(
            "A `GRanges` or `GRangesList` `grna_target_data_frame` must have metadata columns `grna_id` and `grna_target`."
        )
    }
    positioned <- rep(TRUE, length(grna_ranges))
    if (methods::is(grna_ranges, "GRangesList")) {
        n_ranges <- lengths(grna_ranges)
        if (any(n_ranges >= 2L)) {
            stop(
                "Each element of a `GRangesList` `grna_target_data_frame` must contain at most one range (none for a gRNA without a position)."
            )
        }
        positioned <- n_ranges == 1L
        grna_ranges <- unlist(grna_ranges, use.names = FALSE)
    }
    grna_target_data_frame$chr <- NA_character_
    grna_target_data_frame$start <- NA_integer_
    grna_target_data_frame$end <- NA_integer_
    grna_target_data_frame$chr[positioned] <- as.character(
        GenomicRanges::seqnames(grna_ranges)
    )
    grna_target_data_frame$start[positioned] <- GenomicRanges::start(
        grna_ranges
    )
    grna_target_data_frame$end[positioned] <- GenomicRanges::end(grna_ranges)
    return(grna_target_data_frame)
}

# memory
import_data_memory <- function(
    response_matrix,
    grna_matrix,
    grna_target_data_frame,
    moi,
    extra_covariates,
    response_names
) {
    #  perform initial check
    check_import_data_inputs(
        response_matrix,
        grna_matrix,
        grna_target_data_frame,
        moi,
        extra_covariates
    ) |>
        invisible()

    # collapse gRNA matrix (if vector_id supplied)
    if ("vector_id" %in% colnames(grna_target_data_frame)) {
        grna_matrix <- collapse_grna_matrix(
            grna_matrix = grna_matrix,
            grna_target_data_frame = grna_target_data_frame
        )
        grna_target_data_frame <- collapse_grna_target_data_frame(
            grna_target_data_frame
        )
    }

    # process covariates and matrices
    processed_inputs <- process_covariates_and_matrices(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        extra_covariates = extra_covariates,
        response_names = response_names
    )
    # initialize the sceptre_object
    sceptre_object <- init_sceptre_object(
        response_matrix = processed_inputs$response_matrix,
        grna_matrix = processed_inputs$grna_matrix,
        covariate_data_frame = processed_inputs$covariate_data_frame,
        moi = moi,
        grna_target_data_frame = grna_target_data_frame
    )
}

############
# CELLRANGER
############
#' Import data from Cell Ranger
#'
#' `import_data_from_cellranger()` imports data from the output of one or more
#' calls to Cell Ranger count. Each directory supplied as an input to this
#' function should be in feature-barcode format, containing the files
#' `features.tsv.gz` and `matrix.mtx.gz` (and optionally `barcodes.tsv.gz`).
#' When every supplied directory contains a barcode file, the cell barcodes are retained as the column names of the
#' response and gRNA matrices and the row names of the cell covariate data frame.
#' When multiple directories are supplied, the retained cell barcodes are
#' prefixed with their batch identifiers (for example, `b1_` and `b2_`) to keep
#' the cell IDs unique and stable across imports.
#' Supplying barcode files for only some directories is an error.
#' When multiple directories are supplied, the directory from which each cell
#' was loaded is recorded in the cell covariate `cellranger_directory`, a
#' factor that enters the default formula as a covariate. Columns of
#' `extra_covariates`, including one named `batch`, are kept as supplied; with
#' multiple directories, a column of `extra_covariates` named
#' `cellranger_directory` is an error.
#' A malformed matrix file (for example, one with an entry outside the dimensions declared in its header), a matrix file
#' that lists the same (feature, cell) entry more than once, or a matrix file
#' whose number of rows differs from the number of rows of `features.tsv`, is an
#' error.
#' Users can create either a standard `sceptre` object or a disk-backed
#' `sceptre` object, whose response counts stay on disk; the latter is more
#' appropriate for large-scale data. A disk-backed import holds the matrix of
#' one directory in memory at a time and writes the response counts to
#' `screen.h5mu` through temporary files in `directory_to_write` (a
#' subdirectory named `sceptre_import_...`, removed when the import ends). See
#' \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-whole_game_import_data}{the introductory chapter}
#' or
#' \href{https://timothy-barry.github.io/sceptre-book/import-data.html}{Chapter 1}
#' of the manual for more information about this function.
#'
#' @param directories a character vector of file paths to directories containing
#' the output of one or more calls to Cell Ranger count. Each directory should
#' contain the files `matrix.mtx.gz` and `features.tsv.gz` (and optionally
#' `barcodes.tsv.gz`).
#' @param moi a string indicating the MOI of the dataset, either "low" or
#' "high".
#' @param grna_target_data_frame a data frame containing columns `grna_id` and
#' `grna_target` mapping each individual gRNA to its target. Non-targeting gRNAs
#' should be assigned a label of "non-targeting". Optionally,
#' `grna_target_data_frame` can contain columns `chr`, `start`, and `end`,
#' giving the chromosome, start coordinate, and end coordiante, respectively, of
#' each gRNA. Additionally, `grna_target_data_frame` can contain the column
#' `vector_id` specifying the vector to which a given gRNA belongs.
#' @param extra_covariates (optional) a data frame containing extra covariates
#' (e.g., batch, biological replicate) beyond those that `sceptre` can compute.
#' If row names are supplied and cell barcodes are retained, the row names must
#' match the retained (and, for multiple directories, batch-prefixed) cell IDs.
#' @inheritParams import_data
#'
#' @return an initialized `sceptre_object`
#' @export
#' @examples
#' data(grna_target_data_frame_highmoi)
#' directories <- paste0(
#'   system.file("extdata", package = "sceptre"),
#'   "/highmoi_example/gem_group_", c(1, 2)
#' )
#'
#' # 1. create a standard sceptre_object from Cell Ranger output
#' sceptre_object <- import_data_from_cellranger(
#'   directories = directories,
#'   moi = "high",
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#' )
#'
#' # 2. create a disk-backed sceptre_object from Cell Ranger output
#' sceptre_object <- import_data_from_cellranger(
#'   directories = directories,
#'   moi = "high",
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   backed = TRUE,
#'   directory_to_write = tempdir()
#' )
import_data_from_cellranger <- function(
    directories,
    moi,
    grna_target_data_frame,
    extra_covariates = data.frame(),
    backed = FALSE,
    directory_to_write = NULL,
    use_ondisc = NULL
) {
    backed <- resolve_backed(backed, use_ondisc)
    if (backed) {
        directory_to_write <- prepare_directory_to_write(directory_to_write)
    }
    # the directory of each cell is recorded in `cellranger_directory`; do not
    # overwrite a user's column of that name
    if (
        length(directories) >= 2L &&
            "cellranger_directory" %in% colnames(extra_covariates)
    ) {
        stop(
            "`extra_covariates` contains a column named `cellranger_directory`, which `import_data_from_cellranger()` uses to record the directory from which each cell was loaded. Rename that column."
        )
    }
    import_data_from_cellranger_files(
        directories,
        moi,
        grna_target_data_frame,
        extra_covariates,
        directory_to_write = if (backed) directory_to_write else NULL
    )
}

# Reads the matrices into memory; or, given `directory_to_write`, holds one
# directory's matrix in memory at a time: its response counts are written to a
# temporary HDF5 file in `directory_to_write`, and the temporary files are then
# copied, a block of genes at a time, into `directory_to_write/screen.h5mu`,
# from which the returned disk-backed sceptre_object reads them. Both give the
# same sceptre_object.
import_data_from_cellranger_files <- function(
    directories,
    moi,
    grna_target_data_frame,
    extra_covariates,
    directory_to_write = NULL
) {
    backed <- !is.null(directory_to_write)
    # 1. check that directories exist
    for (curr_directory in directories) {
        if (!dir.exists(curr_directory)) {
            stop("The directory ", curr_directory, " does not exist.")
        }
    }

    # 2. create the vector of matrix, features, and barcode files
    input_files <- lapply(X = directories, function(curr_directory) {
        fs <- list.files(curr_directory)
        grep_strs <- c("*features.tsv($|.gz)", "*matrix.mtx($|.gz)")
        out <- vapply(
            grep_strs,
            function(grep_str) {
                file_names <- grep(pattern = grep_str, x = fs, value = TRUE)
                if (length(file_names) >= 2L) {
                    stop(
                        "There are multiple ",
                        grep_str,
                        " files within the directory ",
                        curr_directory,
                        "."
                    )
                }
                if (length(file_names) == 0L) {
                    stop(
                        "The directory ",
                        curr_directory,
                        " contains zero ",
                        grep_str,
                        " files."
                    )
                }
                return(paste0(curr_directory, "/", file_names))
            },
            FUN.VALUE = character(1)
        ) |>
            stats::setNames(c("features", "matrix"))
        barcode_file_names <- grep(
            pattern = "^barcodes\\.tsv(\\.gz)?$",
            x = fs,
            value = TRUE
        )
        if (length(barcode_file_names) >= 2L) {
            stop(
                "There are multiple barcode files within the directory ",
                curr_directory,
                "."
            )
        }
        barcode_fp <- if (length(barcode_file_names) == 1L) {
            file.path(curr_directory, barcode_file_names)
        } else {
            NA_character_
        }
        c(out, barcodes = barcode_fp)
    }) |>
        stats::setNames(NULL)

    # 3. obtain the barcodes, features, and matrix fp vectors
    feature_fps <- vapply(
        X = input_files,
        FUN = function(i) i[["features"]],
        FUN.VALUE = character(1)
    )
    matrix_fps <- vapply(
        X = input_files,
        FUN = function(i) i[["matrix"]],
        FUN.VALUE = character(1)
    )
    barcode_fps <- vapply(
        X = input_files,
        FUN = function(i) i[["barcodes"]],
        FUN.VALUE = character(1)
    )
    barcode_files_present <- !is.na(barcode_fps)
    if (any(barcode_files_present) && !all(barcode_files_present)) {
        missing_directories <- directories[!barcode_files_present]
        stop(
            "Barcode files must either be present in every input directory ",
            "or absent from every input directory. No barcode file was found ",
            "in: ",
            paste(missing_directories, collapse = ", "),
            "."
        )
    }
    retain_cell_barcodes <- all(barcode_files_present)
    cell_barcodes <- if (retain_cell_barcodes) {
        vector(mode = "list", length = length(directories))
    } else {
        NULL
    }

    # 4. obtain the feature data frame; determine the ids of each feature
    feature_df <- data.table::fread(
        file = feature_fps[1],
        colClasses = c("character", "character", "character"),
        col.names = c("feature_id", "feature_name", "modality"),
        header = FALSE
    )
    if (!("CRISPR Guide Capture" %in% feature_df$modality)) {
        stop(
            "The features.tsv file should contain the modality `CRISPR Guide Capture`."
        )
    }
    modalities <- unique(feature_df$modality)
    modalities <- modalities |> stats::setNames(modalities)
    modality_idxs <- lapply(modalities, function(modality) {
        which(feature_df$modality == modality)
    })
    if (!("Gene Expression" %in% modalities)) {
        stop(
            "The features.tsv file should contain the modality `Gene Expression`."
        )
    }
    gene_names <- feature_df$feature_name[modality_idxs[["Gene Expression"]]]

    # 5. loop over directories, loading the matrix corresponding to each; when
    # backed, each matrix's response counts go to a temporary file in
    # `directory_to_write`, with their per-cell covariates computed on the way
    out_list <- vector(mode = "list", length = length(directories))
    n_cells_per_matrix <- vector(mode = "integer", length = length(directories))
    matrix_fps <- path.expand(matrix_fps)
    if (backed) {
        scratch <- tempfile("sceptre_import_", tmpdir = directory_to_write)
        dir.create(scratch)
        on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
        response_files <- file.path(
            scratch,
            paste0("response_", seq_along(directories), ".h5")
        )
        response_covariates <- vector(mode = "list", length(directories))
    }
    for (i in seq_along(directories)) {
        message("Processing directory ", i, ".")
        # check that the features df matches
        features_fp <- feature_fps[i]
        curr_feature_df <- data.table::fread(
            file = features_fp,
            colClasses = c("character", "character", "character"),
            col.names = c("feature_id", "feature_name", "modality"),
            header = FALSE
        )
        if (!isTRUE(all.equal(feature_df, curr_feature_df))) {
            stop("The features.tsv files must match across directories.")
        }

        # process the matrix
        mat <- read_mtx(matrix_fps[i], transpose = FALSE)
        mat_dim <- dim(mat)
        if (mat_dim[1L] != nrow(feature_df)) {
            stop(
                "The matrix file in directory ",
                directories[i],
                " has ",
                mat_dim[1L],
                " rows, but its features.tsv file has ",
                nrow(feature_df),
                "."
            )
        }
        n_cells_per_matrix[i] <- mat_dim[2L]
        if (retain_cell_barcodes) {
            curr_barcodes <- data.table::fread(
                file = barcode_fps[i],
                header = FALSE,
                colClasses = "character",
                showProgress = FALSE
            )
            if (ncol(curr_barcodes) != 1L) {
                stop(
                    "The barcode file in directory ",
                    directories[i],
                    " must contain exactly one column."
                )
            }
            curr_barcodes <- curr_barcodes[[1L]]
            if (length(curr_barcodes) != mat_dim[2L]) {
                barcode_label <- if (length(curr_barcodes) == 1L) {
                    "barcode"
                } else {
                    "barcodes"
                }
                stop(
                    "The barcode file in directory ",
                    directories[i],
                    " contains ",
                    length(curr_barcodes),
                    " ",
                    barcode_label,
                    ", but its matrix contains ",
                    mat_dim[2L],
                    " cells."
                )
            }
            if (anyDuplicated(curr_barcodes)) {
                stop(
                    "The barcode file in directory ",
                    directories[i],
                    " contains duplicate cell barcodes."
                )
            }
            cell_barcodes[[i]] <- curr_barcodes
        }
        # split the rows by feature type; the features of a type need not be
        # contiguous in features.tsv
        subsetted_mats <- lapply(modality_idxs, function(idxs) {
            mat[idxs, , drop = FALSE]
        })
        rm(mat)
        if (backed) {
            # the response counts of this directory are checked and summarized
            # per cell as the in-memory import checks and summarizes the
            # combined counts, then written to their temporary file
            response <- subsetted_mats[["Gene Expression"]]
            check_integer_counts(response, "response_matrix")
            response_covariates[[i]] <- compute_cell_covariates(
                matrix_in = response,
                feature_names = gene_names,
                compute_p_mito = TRUE,
                compute_max_feature = FALSE
            )
            write_backed_counts(response, response_files[i], "counts", level = 0L)
            subsetted_mats[["Gene Expression"]] <- NULL
            rm(response)
        }
        gc() |> invisible()
        out_list[[i]] <- subsetted_mats
        message(crayon::green(" \u2713"))
    }

    # 6. combine the matrices via do.call cbind; when backed, the response
    # counts are read from the temporary files, side by side
    message("Combining matrices across directories.")
    out_mats <- lapply(modalities, function(modality) {
        l <- lapply(out_list, function(mat) mat[[modality]])
        if (is.null(l[[1L]])) {
            # the response counts, when backed: an empty matrix of their size
            return(empty_row_counts(c(
                length(modality_idxs[[modality]]),
                sum(n_cells_per_matrix)
            )))
        }
        do.call(what = "cbind", args = l)
    })
    message(crayon::green(" \u2713"))
    rm(out_list)
    gc() |> invisible()

    # 7. collect metadata pieces
    cell_ids <- if (retain_cell_barcodes) {
        unlist(cell_barcodes, use.names = FALSE)
    } else {
        NULL
    }
    if (!is.null(cell_ids) && length(directories) >= 2L) {
        cell_ids <- unlist(
            Map(
                function(barcodes, batch_id) {
                    paste0("b", batch_id, "_", barcodes)
                },
                cell_barcodes,
                seq_along(cell_barcodes)
            ),
            use.names = FALSE
        )
    }
    for (modality in names(out_mats)) {
        rownames(out_mats[[modality]]) <- feature_df$feature_id[modality_idxs[[
            modality
        ]]]
        if (retain_cell_barcodes) {
            colnames(out_mats[[modality]]) <- cell_ids
        }
    }

    # 8. set up the extra covariates; with multiple directories, record the
    # directory of each cell in `cellranger_directory`, keeping the user's
    # columns (including any named `batch`)
    if (length(directories) >= 2L) {
        cellranger_directory <- rep(
            paste0("b", seq_along(n_cells_per_matrix)),
            n_cells_per_matrix
        ) |>
            factor()
        if (nrow(extra_covariates) == 0L) {
            extra_covariates <- data.frame(
                cellranger_directory = cellranger_directory
            )
        } else {
            extra_covariates$cellranger_directory <- cellranger_directory
        }
    }

    # 9. initialize the sceptre object
    message("Creating the sceptre object.")
    sceptre_object <- if (backed) {
        response_source <- do.call(
            BiocGenerics::cbind,
            lapply(response_files, function(file) {
                DelayedArray::DelayedArray(
                    HDF5Array::H5SparseMatrixSeed(file, "counts")
                )
            })
        )
        dimnames(response_source) <- dimnames(out_mats[["Gene Expression"]])
        import_data_backed(
            response_matrix = out_mats[["Gene Expression"]],
            response_source = response_source,
            response_covariates = do.call(rbind, response_covariates),
            grna_matrix = out_mats[["CRISPR Guide Capture"]],
            grna_target_data_frame = grna_target_data_frame,
            moi = moi,
            extra_covariates = extra_covariates,
            response_names = gene_names,
            directory_to_write = directory_to_write
        )
    } else {
        import_data_memory(
            response_matrix = out_mats[["Gene Expression"]],
            grna_matrix = out_mats[["CRISPR Guide Capture"]],
            grna_target_data_frame = grna_target_data_frame,
            moi = moi,
            extra_covariates = extra_covariates,
            response_names = gene_names
        )
    }
    gc() |> invisible()
    message(crayon::green(" \u2713"))
    return(sceptre_object)
}


# As import_data_memory(), for response counts held in temporary files
# (`response_source`, a DelayedMatrix reading them; `response_matrix` is an
# empty matrix with their dimensions and names) whose per-cell covariates
# (`response_covariates`) were computed as the files were written: the counts
# are copied from the temporary files to `directory_to_write/screen.h5mu`,
# from which the returned sceptre_object reads them.
import_data_backed <- function(
    response_matrix,
    response_source,
    response_covariates,
    grna_matrix,
    grna_target_data_frame,
    moi,
    extra_covariates,
    response_names,
    directory_to_write
) {
    check_import_data_inputs(
        response_matrix,
        grna_matrix,
        grna_target_data_frame,
        moi,
        extra_covariates
    ) |>
        invisible()
    if ("vector_id" %in% colnames(grna_target_data_frame)) {
        grna_matrix <- collapse_grna_matrix(
            grna_matrix = grna_matrix,
            grna_target_data_frame = grna_target_data_frame
        )
        grna_target_data_frame <- collapse_grna_target_data_frame(
            grna_target_data_frame
        )
    }
    if (identical(response_names, NA_character_)) {
        response_names <- rownames(response_matrix)
    }
    covariate_data_frame <- auto_compute_cell_covariates(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        extra_covariates = extra_covariates,
        response_names = response_names,
        response_covariates = response_covariates
    )
    sceptre_object <- init_sceptre_object(
        response_matrix = response_matrix,
        grna_matrix = set_matrix_accessibility(
            grna_matrix,
            make_row_accessible = TRUE
        ),
        covariate_data_frame = covariate_data_frame,
        moi = moi,
        grna_target_data_frame = grna_target_data_frame
    )
    file <- file.path(directory_to_write, "screen.h5mu")
    write_screen(sceptre_object@data, file, response_source = response_source)
    backed_counts <- read_screen(file, backed = TRUE)$response
    dimnames(backed_counts) <- dimnames(response_matrix)
    set_response_matrix(sceptre_object, backed_counts)
}


#######
# PARSE
#######
#' Import data from Parse (experimental)
#'
#' `import_data_from_parse()` imports data from the output of the Parse count
#' matrix generation program. See
#' \href{https://timothy-barry.github.io/sceptre-book/import-data.html#import-from-the-parse-program-experimental}{Chapter 1 of the manual}
#' for more information about this function.
#'
#' `import_data_from_parse()` is experimental, and the API of this function is
#' subject to change. We expect the API to solidify as we learn more about the
#' Parse platform and the structure of the Parse count matrix generation program
#' output.
#'
#' A malformed matrix file (for example, one with an entry outside the
#' dimensions declared in its header), or one that lists the same
#' (cell, feature) entry more than once, is an error.
#'
#' The example data shipped with the package (used in the examples below) are a small subset of Parse Biosciences' "Technical Performance of CRISPR Detect in Cell Lines" demo dataset, used under the Creative Commons Attribution 4.0 International (CC BY 4.0) license, with attribution to Parse Biosciences. See \url{https://www.parsebiosciences.com/datasets/technical-performance-of-crispr-detect-in-cell-lines/}.
#'
#' @param gene_mat_fp file path to the gene `count_matrix.mtx` file.
#' @param grna_mat_fp file path to the gRNA `count_matrix.mtx` file.
#' @param all_genes_fp file path to the `all_genes.csv` file containing the gene
#' IDs.
#' @param all_grnas_fp file path to the `all_guides.csv` file containing the
#' gRNA IDs. The gRNA IDs are assumed to be in the second column (i.e., the
#' "gene_name" column) of this file.
#' @param moi a string indicating the MOI of the dataset, either "low" or
#' "high".
#' @param grna_target_data_frame a data frame containing columns `grna_id` and
#' `grna_target` mapping each individual gRNA to its target. Non-targeting gRNAs
#' should be assigned a label of "non-targeting". Optionally,
#' `grna_target_data_frame` can contain columns `chr`, `start`, and `end`,
#' giving the chromosome, start coordinate, and end coordiante, respectively, of
#' each gRNA. Additionally, `grna_target_data_frame` can contain the column
#' `vector_id` specifying the vector to which a given gRNA belongs.
#' @param extra_covariates (optional) a data frame containing extra covariates
#' (e.g., batch, biological replicate) beyond those that `sceptre` can compute.
#'
#' @return an initialized `sceptre_object`
#' @export
#'
#' @examples
#' directory <- paste0(
#'   system.file("extdata", package = "sceptre"),
#'   "/parse_example/"
#' )
#' gene_mat_fp <- paste0(directory, "gene_mat.mtx")
#' grna_mat_fp <- paste0(directory, "grna_mat.mtx")
#' all_genes_fp <- paste0(directory, "all_genes.csv")
#' all_grnas_fp <- paste0(directory, "all_grnas.csv")
#' grna_target_data_frame <- data.frame(
#'   grna_id = c("guide_A", "guide_B", "guide_C"),
#'   grna_target = c("target-A", "target-B", "non-targeting")
#' )
#' sceptre_object <- import_data_from_parse(
#'   gene_mat_fp = gene_mat_fp,
#'   grna_mat_fp = grna_mat_fp,
#'   all_genes_fp = all_genes_fp,
#'   all_grnas_fp = all_grnas_fp,
#'   moi = "low",
#'   grna_target_data_frame = grna_target_data_frame
#' )
import_data_from_parse <- function(
    gene_mat_fp,
    grna_mat_fp,
    all_genes_fp,
    all_grnas_fp,
    moi,
    grna_target_data_frame,
    extra_covariates = data.frame()
) {
    # 1. check that files exist
    fps <- c(
        gene_mat_fp = gene_mat_fp,
        grna_mat_fp = grna_mat_fp,
        all_genes_fp = all_genes_fp,
        all_grnas_fp = all_grnas_fp
    )
    for (fp_name in names(fps)) {
        fp <- fps[[fp_name]]
        if (!file.exists(fp)) {
            stop(
                "File ",
                fp,
                " does not exist. Set `",
                fp_name,
                "` to a different value."
            )
        }
    }

    # 2. load the gene.mtx and grna.mtx data; Parse writes cells in rows
    gene_mat <- read_mtx(path.expand(gene_mat_fp), transpose = TRUE)
    grna_mat <- read_mtx(path.expand(grna_mat_fp), transpose = TRUE)

    # 3. load the gene ids; update rownames of gene mat
    gene_feature_df <- data.table::fread(
        file = all_genes_fp,
        colClasses = c("character", "character", "character"),
        col.names = c("feature_id", "feature_name", "genome"),
        header = TRUE
    )
    rownames(gene_mat) <- gene_feature_df$feature_id[seq(1L, nrow(gene_mat))]
    gene_names <- gene_feature_df$feature_name[seq(1L, nrow(gene_mat))]

    # 4. load the grna ids; update rownames of grna mat
    grna_feature_df <- data.table::fread(
        file = all_grnas_fp,
        colClasses = c("character", "character", "character"),
        col.names = c("feature_id", "feature_name", "genome"),
        header = TRUE
    )
    rownames(grna_mat) <- grna_feature_df$feature_name

    # 5. create the sceptre_object
    sceptre_object <- import_data_memory(
        response_matrix = gene_mat,
        grna_matrix = grna_mat,
        grna_target_data_frame = grna_target_data_frame,
        moi = moi,
        extra_covariates = extra_covariates,
        response_names = gene_names
    )
    return(sceptre_object)
}


##################
# HELPER FUNCTIONS
##################
init_sceptre_object <- function(
    response_matrix,
    grna_matrix,
    covariate_data_frame,
    moi,
    grna_target_data_frame
) {
    # create sceptre_object
    sceptre_object <- methods::new("sceptre_object")

    # process the covariate data frame
    sceptre_object@import_grna_assignment_info <- list(
        max_grna = covariate_data_frame$grna_feature_w_max_expression,
        max_grna_frac_umis = covariate_data_frame$grna_frac_umis_max_feature
    )
    covariate_data_frame$grna_feature_w_max_expression <- covariate_data_frame$grna_frac_umis_max_feature <- NULL

    # set data fields
    sceptre_object@data <- new_screen(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        covariate_data_frame = covariate_data_frame
    )
    sceptre_object@low_moi <- (moi == "low")
    sceptre_object@covariate_names <- sort(colnames(covariate_data_frame))
    sceptre_object@grna_target_data_frame <- grna_target_data_frame |>
        dplyr::mutate(
            grna_id = as.character(grna_id),
            grna_target = as.character(grna_target)
        )
    # 5. initialize flags
    sceptre_object@last_function_called <- "import_data"
    sceptre_object@functs_called <- c(
        import_data = TRUE,
        set_analysis_parameters = FALSE,
        assign_grnas = FALSE,
        run_qc = FALSE,
        run_calibration_check = FALSE,
        run_power_check = FALSE,
        run_discovery_analysis = FALSE
    )
    return(sceptre_object)
}


prepare_directory_to_write <- function(directory_to_write) {
    if (is.null(directory_to_write)) {
        stop("`directory_to_write` must be supplied.")
    }
    if (
        !(is.character(directory_to_write) &&
            length(directory_to_write) == 1L &&
            !is.na(directory_to_write))
    ) {
        stop("`directory_to_write` must be a single file path.")
    }
    directory_to_write <- path.expand(directory_to_write)
    if (!dir.exists(directory_to_write)) {
        dir.create(directory_to_write, recursive = TRUE)
    }
    directory_to_write <- normalizePath(
        directory_to_write,
        winslash = "/",
        mustWork = TRUE
    )
    check_path_has_no_tilde(directory_to_write, "`directory_to_write`")
    return(directory_to_write)
}


prepare_file_path <- function(file_path, arg_name) {
    if (is.null(file_path)) {
        stop(arg_name, " must be supplied.")
    }
    if (
        !(is.character(file_path) &&
            length(file_path) == 1L &&
            !is.na(file_path))
    ) {
        stop(arg_name, " must be a single file path.")
    }
    file_path <- normalizePath(
        path.expand(file_path),
        winslash = "/",
        mustWork = TRUE
    )
    check_path_has_no_tilde(file_path, arg_name)
    return(file_path)
}


check_path_has_no_tilde <- function(file_path, arg_name) {
    if (grepl("~", file_path, fixed = TRUE)) {
        stop(
            arg_name,
            " resolves to '",
            file_path,
            "', which contains a '~' character ",
            "inside a directory or file name. HDF5 misinterprets embedded tildes as ",
            "user-home expansion and produces a corrupted path. Please choose a path ",
            "where no directory or file name contains '~' (a leading '~/' is fine - ",
            "it is expanded before this check).",
            call. = FALSE
        )
    }
    return(NULL)
}


process_covariates_and_matrices <- function(
    response_matrix,
    grna_matrix,
    extra_covariates,
    response_names
) {
    covariate_data_frame <- auto_compute_cell_covariates(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        extra_covariates = extra_covariates,
        response_names = if (identical(response_names, NA_character_)) {
            rownames(response_matrix)
        } else {
            response_names
        }
    )
    response_matrix <- set_matrix_accessibility(
        response_matrix,
        make_row_accessible = TRUE
    )
    grna_matrix <- set_matrix_accessibility(
        grna_matrix,
        make_row_accessible = TRUE
    )
    return(list(
        covariate_data_frame = covariate_data_frame,
        response_matrix = response_matrix,
        grna_matrix = grna_matrix
    ))
}


collapse_grna_target_data_frame <- function(grna_target_data_frame) {
    vector_nested_within_target <- (grna_target_data_frame |>
        dplyr::group_by(vector_id) |>
        dplyr::summarize(n_distinct_targets = length(unique(grna_target))) |>
        dplyr::pull(n_distinct_targets) ==
        1) |>
        all()
    if (!vector_nested_within_target) {
        stop("`vector_id` must be nested within `grna_target`.")
    }

    grna_target_data_frame |>
        dplyr::select(-grna_id) |>
        dplyr::rename(grna_id = vector_id) |>
        dplyr::distinct()
}


collapse_grna_matrix <- function(grna_matrix, grna_target_data_frame) {
    # temporary function; we may rewrite this function to improve its
    # speed and memory efficiency if multiguide vector data become more
    # prevalent.
    vector_ids <- grna_target_data_frame$vector_id |> unique()
    grna_matrix <- as.matrix(grna_matrix)
    n_cells <- ncol(grna_matrix)
    grna_matrix_collapsed <- vapply(
        vector_ids,
        function(curr_vector_id) {
            curr_grna_ids <- grna_target_data_frame |>
                dplyr::filter(vector_id == curr_vector_id) |>
                dplyr::pull(grna_id)
            Matrix::colSums(grna_matrix[curr_grna_ids, , drop = FALSE])
        },
        FUN.VALUE = numeric(n_cells)
    ) |>
        t()
    return(grna_matrix_collapsed)
}
