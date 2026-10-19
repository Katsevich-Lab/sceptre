############
# R MATRICES
############
#' Import data
#'
#' `import_data()` imports data from a collection of R objects to create a
#' `sceptre_object`. Users can create either a standard `sceptre` object or an
#' `ondisc`-backed `sceptre` object; the latter is more appropriate for
#' large-scale data. See
#' \href{https://timothy-barry.github.io/sceptre-book/import-data.html#import-data-from-a-collection-of-r-objects}{Chapter 1 of the manual}
#' for more detailed information about this function.
#'
#' @param response_matrix a matrix of response UMI counts, with responses in
#' rows and cells in columns. The matrix should be of type `"matrix"`,
#' `"dgCMatrix"`, `"dgRMatrix"`, or `"dgTMatrix"`. The row names of the matrix
#' should give the response IDs. The column names, when supplied, should give
#' the cell IDs.
#' @param grna_matrix a matrix of gRNA UMI counts, with gRNAs in rows and cells
#' in columns. The matrix should be of type `"matrix"`, `"dgCMatrix"`,
#' `"dgRMatrix"`, or `"dgTMatrix"`. The row names of the matrix should give the
#' gRNA IDs. The column names, when supplied, should give the cell IDs.
#' @param grna_target_data_frame a data frame containing columns `grna_id` and
#' `grna_target` mapping each individual gRNA to its target. Non-targeting gRNAs
#' should be assigned a label of "non-targeting". Optionally,
#' `grna_target_data_frame` can contain columns `chr`, `start`, and `end`,
#' giving the chromosome, start coordinate, and end coordiante, respectively, of
#' each gRNA. Additionally, `grna_target_data_frame` can contain the column
#' `vector_id` specifying the vector to which a given gRNA belongs.
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
#' @param use_ondisc (default `FALSE`) a logical value (i.e., `TRUE` or `FALSE`)
#' indicating whether to create an `ondisc`-backed `sceptre_object` (`TRUE`) or
#' a standard `sceptre_object` (`FALSE`).
#' @param directory_to_write (optional) a file path to a directory in which to
#' write the backing `.odm` files for the response and gRNA expression matrices.
#' Must be supplied if `use_ondisc` is set to `TRUE`.
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
#' # 2. initialize an ondisc-backed sceptre_object from R objects
#' if (requireNamespace("ondisc", quietly = TRUE)) {
#'   sceptre_object <- import_data(
#'     response_matrix = lowmoi_example_data$response_matrix,
#'     grna_matrix = lowmoi_example_data$grna_matrix,
#'     grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'     extra_covariates = lowmoi_example_data$extra_covariates,
#'     moi = "low",
#'     use_ondisc = TRUE,
#'     directory_to_write = tempdir()
#'   )
#' }
import_data <- function(
    response_matrix,
    grna_matrix,
    grna_target_data_frame,
    moi,
    extra_covariates = data.frame(),
    response_names = NA_character_,
    use_ondisc = FALSE,
    directory_to_write = NULL
) {
    if (use_ondisc) {
        # use ondisc
        import_data_use_ondisc(
            response_matrix,
            grna_matrix,
            grna_target_data_frame,
            moi,
            extra_covariates,
            response_names,
            directory_to_write
        )
    } else {
        # do not use ondisc
        import_data_memory(
            response_matrix,
            grna_matrix,
            grna_target_data_frame,
            moi,
            extra_covariates,
            response_names
        )
    }
}

# memory
#' Import in-memory count matrices and create a sceptre object
#' @param response_matrix,grna_matrix Feature-by-cell matrices with unique feature
#'   names and cells in the same order. Supplied cell names must also match.
#' @param grna_target_data_frame Guide-to-target metadata. If vector_id is present,
#'   sum guides within each vector before calculating cell covariates.
#' @param moi `"low"` or `"high"`.
#' @param extra_covariates Covariates in the same cell order, or an empty data frame.
#' @param response_names Feature labels matching response rows; NA_character_ uses
#'   row names when identifying mitochondrial genes.
#' @return Initialized object with row-accessible sparse counts and per-cell
#'   covariates, retaining original cell order. Conversion may allocate new matrices.
#' @noRd
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
        response_names = processed_inputs$response_names,
        moi = moi,
        grna_target_data_frame = grna_target_data_frame
    )
}

#' Import R count matrices by writing them to ondisc files
#' @param response_matrix,grna_matrix,grna_target_data_frame,moi,extra_covariates,response_names
#'   Inputs as for `import_data_memory()`. In this path, guide rows are written
#'   as supplied; they are not first summed by vector_id.
#' @param directory_to_write Output directory, created if necessary.
#' @return Initialized object and backing files response.odm and grna.odm, sharing
#'   a newly sampled integer_id. Sampling this ID advances R's random state;
#'   matrix construction follows the existing response-then-guide order.
#' @noRd
import_data_use_ondisc <- function(
    response_matrix,
    grna_matrix,
    grna_target_data_frame,
    moi,
    extra_covariates,
    response_names,
    directory_to_write
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
    # check directory_to_write
    directory_to_write <- prepare_directory_to_write(directory_to_write)

    # process covariates and matrices
    processed_inputs <- process_covariates_and_matrices(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        extra_covariates = extra_covariates,
        response_names = response_names
    )
    response_matrix <- processed_inputs$response_matrix
    grna_matrix <- processed_inputs$grna_matrix
    covariate_data_frame <- processed_inputs$covariate_data_frame

    # write matrices to disk
    integer_id <- sample(x = seq(0L, .Machine$integer.max), size = 1L)
    out <- write_matrices_to_disk(
        response_matrix = response_matrix,
        grna_matrix = grna_matrix,
        directory_to_write = directory_to_write,
        integer_id = integer_id
    )

    # initialize the sceptre_object
    sceptre_object <- init_sceptre_object(
        response_matrix = out$response_matrix,
        grna_matrix = out$grna_matrix,
        covariate_data_frame = covariate_data_frame,
        response_names = processed_inputs$response_names,
        moi = moi,
        grna_target_data_frame = grna_target_data_frame
    )

    # add integer id
    return(sceptre_object)
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
#' For standard in-memory objects, when every supplied directory contains a
#' barcode file, the cell barcodes are retained as the column names of the
#' response and gRNA matrices and the row names of the cell covariate data
#' frame.
#' When multiple directories are supplied, the retained cell barcodes are
#' prefixed with their batch identifiers (for example, `b1_` and `b2_`) to keep
#' the cell IDs unique and stable across imports.
#' Supplying barcode files for only some directories is an error.
#' Users can create either a standard `sceptre` object or an `ondisc`-backed
#' `sceptre` object; the latter is more appropriate for large-scale data. See
#' \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-whole_game_import_data}{the introductory chapter}
#' or
#' \href{https://timothy-barry.github.io/sceptre-book/import-data.html}{Chapter 1}
#' of the manual for more information about this function.
#'
#' For in-memory import, count-file entries may be unsorted. Each cell-feature
#' coordinate must occur only once; repeated coordinates cause an import error
#' rather than being summed. This does not reject repeated guide-metadata rows.
#' Backed Cell Ranger import delegates file parsing to ondisc, which has its own
#' coordinate handling. See the README for the qualified ondisc installation.
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
#' @param use_ondisc (optional; default `FALSE`) a logical indicating whether to
#' store the expression data in a disk-backed `ondisc` matrix (`TRUE`) or an
#' in-memory sparse matrix (`FALSE`).
#' @param directory_to_write (optional) a string indicating the directory in
#' which to write the backing `.odm` files (must be specified if `use_ondisc` is
#' set to `TRUE`).
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
#' # 2. create an ondisc-backed sceptre_object from Cell Ranger output
#' if (requireNamespace("ondisc", quietly = TRUE)) {
#'   sceptre_object <- import_data_from_cellranger(
#'     directories = directories,
#'     moi = "high",
#'     grna_target_data_frame = grna_target_data_frame_highmoi,
#'     use_ondisc = TRUE,
#'     directory_to_write = tempdir()
#'   )
#' }
import_data_from_cellranger <- function(
    directories,
    moi,
    grna_target_data_frame,
    extra_covariates = data.frame(),
    use_ondisc = FALSE,
    directory_to_write = NULL
) {
    # take cases on use_ondisc
    if (!use_ondisc) {
        import_data_from_cellranger_memory(
            directories,
            moi,
            grna_target_data_frame,
            extra_covariates
        )
    } else {
        import_data_from_cellranger_use_ondisc(
            directories,
            moi,
            grna_target_data_frame,
            extra_covariates,
            directory_to_write
        )
    }
}

#' Read Cell Ranger batches into in-memory response and guide matrices
#' @param directories Batch directories in desired cell order. Feature tables must
#'   match across batches. Matrix Market entries may arrive unsorted, and response
#'   and guide rows may be interleaved in the feature file. Reject repeated
#'   (feature, cell) entries before building the sparse matrix.
#' @param moi,grna_target_data_frame,extra_covariates Inputs to import_data_memory();
#'   extra covariates follow concatenated batch/cell order.
#' @return Initialized object with cells in directory order and features in their
#'   feature-file order. Supply barcode files for every batch or none. For multiple
#'   batches, prefix retained barcodes with batch IDs. Add or replace the batch
#'   covariate and report progress.
#' @noRd
import_data_from_cellranger_memory <- function(
    directories,
    moi,
    grna_target_data_frame,
    extra_covariates
) {
    for (curr_directory in directories) {
        if (!dir.exists(curr_directory)) {
            stop("The directory ", curr_directory, " does not exist.")
        }
    }

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

    out_list <- vector(mode = "list", length = length(directories))
    n_cells_per_matrix <- vector(mode = "integer", length = length(directories))
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
        for (col in colnames(feature_df)) {
            if (
                any(
                    dplyr::pull(feature_df, dplyr::all_of(col)) !=
                        dplyr::pull(curr_feature_df, dplyr::all_of(col))
                )
            ) {
                stop("The features.tsv files must match across directories.")
            }
        }

        # process the matrix
        matrix_fp <- matrix_fps[i]
        matrix_metadata <- get_mtx_metadata(matrix_fp)
        n_cells_per_matrix[i] <- matrix_metadata$n_cells
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
            if (length(curr_barcodes) != matrix_metadata$n_cells) {
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
                    matrix_metadata$n_cells,
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
        dt <- data.table::fread(
            file = matrix_fp,
            skip = matrix_metadata$n_to_skip,
            col.names = c("feature_idx", "cell_idx", "x"),
            colClasses = c("integer", "integer", "double"),
            showProgress = FALSE,
            nThread = 2L
        )
        subsetted_mats <- lapply(modalities, function(modality) {
            curr_feature_idxs <- modality_idxs[[modality]]
            # subset
            dt_sub <- dt[dt$feature_idx %in% curr_feature_idxs, ]
            # Genes and guides may be interleaved in the feature file. Find
            # each feature's row within this type's own list, preserving file
            # order, then subtract one for Matrix's internal row numbers.
            dt_sub$feature_idx <- match(dt_sub$feature_idx, curr_feature_idxs) -
                1L
            dt_sub$cell_idx <- dt_sub$cell_idx - 1L
            # create a sparse matrix in csc format
            out <- create_csc_matrix(
                dt = dt_sub,
                n_cells = matrix_metadata$n_cells,
                n_features = length(curr_feature_idxs)
            )
            return(out)
        })
        rm(dt)
        gc() |> invisible()
        out_list[[i]] <- subsetted_mats
        message(crayon::green(" \u2713"))
    }

    message("Combining matrices across directories.")
    out_mats <- lapply(modalities, function(modality) {
        l <- lapply(out_list, function(mat) mat[[modality]])
        combined_mat <- do.call(what = "cbind", args = l)
        return(combined_mat)
    })
    rm(out_list)
    gc() |> invisible()
    message(crayon::green(" \u2713"))

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
    for (modality in modalities) {
        rownames(out_mats[[modality]]) <- feature_df$feature_id[modality_idxs[[
            modality
        ]]]
        if (retain_cell_barcodes) {
            colnames(out_mats[[modality]]) <- cell_ids
        }
    }
    gene_names <- feature_df$feature_name[modality_idxs[["Gene Expression"]]]

    batch <- if (length(directories) >= 2L) {
        rep(paste0("b", seq_along(n_cells_per_matrix)), n_cells_per_matrix) |>
            factor()
    } else {
        NULL
    }
    if (nrow(extra_covariates) == 0L) {
        extra_covariates <- data.frame(batch = batch)
    } else {
        extra_covariates$batch <- batch
    }

    message("Creating the sceptre object.")
    sceptre_object <- import_data_memory(
        response_matrix = out_mats[["Gene Expression"]],
        grna_matrix = out_mats[["CRISPR Guide Capture"]],
        grna_target_data_frame = grna_target_data_frame,
        moi = moi,
        extra_covariates = extra_covariates,
        response_names = gene_names
    )
    gc() |> invisible()
    message(crayon::green(" \u2713"))
    return(sceptre_object)
}


#' Read Cell Ranger batches using the ondisc importer
#' @param directories Batch directories passed to ondisc in supplied order.
#' @param moi,grna_target_data_frame,extra_covariates Import settings and metadata;
#'   extra covariate rows must follow the cell order returned by ondisc.
#' @param directory_to_write Output directory, created if necessary.
#' @return Object with response and guide backing files sharing one import ID.
#'   ondisc writes the files, computes cell covariates and combines vector_id
#'   counts; update guide metadata to match its output. ondisc also parses and
#'   checks file coordinates itself. SCEPTRE's duplicate-coordinate rejection in
#'   the in-memory importer does not govern this separate parser.
#' @noRd
import_data_from_cellranger_use_ondisc <- function(
    directories,
    moi,
    grna_target_data_frame,
    extra_covariates,
    directory_to_write
) {
    directory_to_write <- prepare_directory_to_write(directory_to_write)

    vector_supplied <- "vector_id" %in% colnames(grna_target_data_frame)
    out <- ondisc::create_odm_from_cellranger(
        directories_to_load = directories,
        directory_to_write = directory_to_write,
        write_cellwise_covariates = FALSE,
        grna_target_data_frame = if (vector_supplied) {
            grna_target_data_frame
        } else {
            NULL
        }
    )

    if (vector_supplied) {
        grna_target_data_frame <- collapse_grna_target_data_frame(
            grna_target_data_frame
        )
    }

    check_import_data_inputs(
        out$gene,
        out$grna,
        grna_target_data_frame,
        moi,
        extra_covariates
    ) |>
        invisible()

    covariate_df <- out$cellwise_covariates
    colnames(covariate_df) <- gsub(
        pattern = "gene",
        replacement = "response",
        fixed = TRUE,
        x = colnames(covariate_df)
    )
    if (nrow(extra_covariates) >= 1L) {
        covariate_df <- cbind(extra_covariates, covariate_df)
    }

    sceptre_object <- init_sceptre_object(
        response_matrix = out$gene,
        grna_matrix = out$grna,
        covariate_data_frame = covariate_df,
        moi = moi,
        grna_target_data_frame = grna_target_data_frame
    )

    return(sceptre_object)
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
#' The example data shipped with the package (used in the examples below) are a small subset of Parse Biosciences' "Technical Performance of CRISPR Detect in Cell Lines" demo dataset, used under the Creative Commons Attribution 4.0 International (CC BY 4.0) license, with attribution to Parse Biosciences. See \url{https://www.parsebiosciences.com/datasets/technical-performance-of-crispr-detect-in-cell-lines/}.
#'
#' Count-file entries may be unsorted, but each cell-feature coordinate must
#' occur only once. Repeated coordinates cause an import error rather than being
#' summed; repeated guide-metadata rows are a separate matter.
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

    # Parse lists each entry as (cell, feature, count). Read into a new table
    # and subtract one from file positions for Matrix's internal numbering.
    # create_csc_matrix() then sorts by cell/feature and rejects duplicate
    # entries. The resulting matrix has features in rows and cells in columns.
    load_matrix <- function(mat_fp) {
        matrix_metadata <- get_mtx_metadata(
            mat_fp,
            c("n_cells", "n_features", "n_nonzero")
        )
        dt <- data.table::fread(
            file = mat_fp,
            skip = matrix_metadata$n_to_skip,
            col.names = c("cell_idx", "feature_idx", "x"),
            colClasses = c("integer", "integer", "double"),
            showProgress = FALSE,
            nThread = 2L
        )
        # make the cell and feature idxs zero-based
        dt$cell_idx <- dt$cell_idx - 1L
        dt$feature_idx <- dt$feature_idx - 1L

        # create a sparse matrix in csc format
        out <- create_csc_matrix(
            dt = dt,
            n_cells = matrix_metadata$n_cells,
            n_features = matrix_metadata$n_features
        )
        return(out)
    }
    gene_mat <- load_matrix(gene_mat_fp)
    grna_mat <- load_matrix(grna_mat_fp)

    gene_feature_df <- data.table::fread(
        file = all_genes_fp,
        colClasses = c("character", "character", "character"),
        col.names = c("feature_id", "feature_name", "genome"),
        header = TRUE
    )
    rownames(gene_mat) <- gene_feature_df$feature_id[seq(1L, nrow(gene_mat))]
    gene_names <- gene_feature_df$feature_name[seq(1L, nrow(gene_mat))]

    grna_feature_df <- data.table::fread(
        file = all_grnas_fp,
        colClasses = c("character", "character", "character"),
        col.names = c("feature_id", "feature_name", "genome"),
        header = TRUE
    )
    rownames(grna_mat) <- grna_feature_df$feature_name

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
#' Create the saved sceptre object from prepared matrices and covariates
#' @param response_matrix,grna_matrix Validated feature-by-cell matrices with cells
#'   in matching order, stored for row access in memory or in ondisc.
#' @param covariate_data_frame Covariates in the same cell order. Move the guide
#'   maximum/fraction summaries into import_grna_assignment_info.
#' @param moi `"low"` or `"high"`.
#' @param grna_target_data_frame Guide metadata. Convert guide IDs and targets to
#'   character while retaining row order and other columns.
#' @return New object with its existing field layout and initial workflow flags.
#'   Save the supplied matrices without reordering cells or loading all counts.
#' @noRd
init_sceptre_object <- function(
    response_matrix,
    grna_matrix,
    covariate_data_frame,
    moi,
    grna_target_data_frame,
    response_names = character()
) {
    # create sceptre_object
    sceptre_object <- methods::new("sceptre_object")

    sceptre_object@data@response_names <- as.character(response_names)
    # process the covariate data frame
    sceptre_object@data@import_grna_assignment_info <- list(
        max_grna = covariate_data_frame$grna_feature_w_max_expression,
        max_grna_frac_umis = covariate_data_frame$grna_frac_umis_max_feature
    )
    covariate_data_frame$grna_feature_w_max_expression <- covariate_data_frame$grna_frac_umis_max_feature <- NULL

    # set data fields
    sceptre_object <- set_response_matrix(sceptre_object, response_matrix)
    sceptre_object <- set_grna_matrix(sceptre_object, grna_matrix)
    sceptre_object@data@covariate_data_frame <- covariate_data_frame
    sceptre_object@data@low_moi <- (moi == "low")
    sceptre_object@data@grna_target_data_frame <- grna_target_data_frame |>
        dplyr::mutate(
            grna_id = as.character(grna_id),
            grna_target = as.character(grna_target)
        )
    sceptre_object <- initialize_screen_identity(sceptre_object)
    validObject(sceptre_object)
    return(sceptre_object)
}


#' Build a sparse feature-by-cell matrix from file entries
#'
#' Each row of dt gives cell_idx, feature_idx and value x. Positions already start
#' at zero. Make new output vectors sorted by cell, then feature, without changing
#' dt. Keep every value, including explicitly stored zeros. Reject repeated
#' (feature, cell) positions rather than summing them or choosing one value.
#' Check that dimensions and entry counts fit C++ integers before building the
#' column pointers. Return a numeric dgCMatrix without row or column names;
#' features and cells keep their file order.
#' @noRd
create_csc_matrix <- function(dt, n_cells, n_features) {
    check_native_scalar_count(n_cells, "n_cells")
    check_native_scalar_count(n_features, "n_features")
    if (n_cells == .Machine$integer.max) {
        stop(
            "The imported cell count leaves no room for its final sparse pointer."
        )
    }
    check_native_scalar_count(nrow(dt), "Stored coordinate count")
    if (
        !is.integer(dt$cell_idx) ||
            !is.integer(dt$feature_idx) ||
            anyNA(dt$cell_idx) ||
            anyNA(dt$feature_idx) ||
            any(dt$cell_idx < 0L | dt$cell_idx >= n_cells) ||
            any(dt$feature_idx < 0L | dt$feature_idx >= n_features)
    ) {
        stop(
            "Imported cell or feature coordinates are outside the matrix dimensions."
        )
    }
    scan_order <- order(dt$cell_idx, dt$feature_idx, method = "radix")
    cell_idxs <- dt$cell_idx[scan_order]
    feature_idxs <- dt$feature_idx[scan_order]
    if (any(diff(cell_idxs) == 0L & diff(feature_idxs) == 0L)) {
        stop(
            "Repeated (cell, feature) coordinates in imported count file; each coordinate must be unique."
        )
    }
    methods::new(
        "dgCMatrix",
        i = feature_idxs,
        p = obtain_pointer_vector(i = cell_idxs, dim = as.integer(n_cells)),
        x = dt$x[scan_order],
        Dim = as.integer(c(n_features, n_cells))
    )
}


#' Read the dimensions at the start of a Matrix Market file
#'
#' Skip header/comment lines and read the first dimensions line, using the
#' existing space-separated parser for plain or gzip input. Return the three
#' integer dimensions named by col_id and n_to_skip, the total header/comment/
#' dimensions lines to skip when reading values. Close the file on success or
#' error. Other malformed-file checks are not added here.
#' @noRd
get_mtx_metadata <- function(
    mtx_file,
    col_id = c("n_features", "n_cells", "n_nonzero")
) {
    con <- file(mtx_file, "r")
    on.exit(close(con), add = TRUE)
    n_to_skip <- 0L
    while (TRUE) {
        n_to_skip <- n_to_skip + 1L
        line <- readLines(con, n = 1)
        if (substr(line, 0, 1) != "%") break()
    }
    metrics <- strsplit(line, split = " ", fixed = TRUE)[[1]] |> as.integer()
    out <- list(metrics[1], metrics[2], metrics[3], n_to_skip)
    names(out) <- c(col_id, "n_to_skip")
    return(out)
}


#' Write response and guide counts with the same import ID
#' @param response_matrix,grna_matrix Feature-by-cell matrices with matching cells.
#' @param directory_to_write Output directory from prepare_directory_to_write().
#' @param integer_id ID shared by the saved object and both count files.
#' @return Open response and guide matrices from storage_create_backed_pair(),
#'   after writing response.odm followed by grna.odm.
#' @noRd
write_matrices_to_disk <- function(
    response_matrix,
    grna_matrix,
    directory_to_write,
    integer_id
) {
    directory_to_write <- prepare_directory_to_write(directory_to_write)
    storage_create_backed_pair(
        response_matrix,
        grna_matrix,
        directory_to_write,
        integer_id
    )
}


#' Create the output directory and return its absolute path
#' @param directory_to_write One nonmissing path. Expand a leading ~ for the home
#'   directory and create parent directories as needed.
#' @return Existing absolute path with forward slashes. Reject any remaining ~
#'   before passing the path to ondisc, whose path handling could misread it.
#' @noRd
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


#' Expand and check an existing object or count-file path
#' @param file_path One nonmissing path that must already exist.
#' @param arg_name Argument name used in error messages.
#' @return Absolute path with forward slashes and the home-directory abbreviation
#'   expanded. Reject remaining ~ characters. Do not open or modify the file.
#' @noRd
prepare_ondisc_file_path <- function(file_path, arg_name) {
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


#' Reject remaining ~ characters before passing a path to ondisc
#' @param file_path Path whose leading home-directory abbreviation is already expanded.
#' @param arg_name Argument name used in the error message.
#' @return NULL if no ~ remains, otherwise an error. This checks the text only.
#' @noRd
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


#' Compute cell covariates and prepare count matrices for row access
#' @param response_matrix,grna_matrix In-memory feature-by-cell counts.
#' @param extra_covariates Extra covariates in the same cell order, or an empty table.
#' @param response_names Response labels matching rows, or NA_character_ for row names.
#' @return covariate_data_frame, response_matrix and grna_matrix. Keep dimensions,
#'   feature names and cell order while converting count storage as needed.
#'   The backed R-matrix importer also calls this before writing files.
#' @noRd
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
        response_names = if (identical(response_names, NA_character_)) rownames(response_matrix) else response_names,
        covariate_data_frame = covariate_data_frame,
        response_matrix = response_matrix,
        grna_matrix = grna_matrix
    ))
}


#' Replace individual guide IDs with the vector IDs used to combine their counts
#' @param grna_target_data_frame Metadata with grna_id, vector_id and grna_target.
#'   Every vector must map to exactly one target.
#' @return Metadata with vector_id renamed to grna_id and the old guide IDs removed.
#'   Keep the first occurrence of each distinct remaining row. Different extra
#'   annotations can therefore leave repeated vector IDs in the output.
#' @noRd
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


#' Sum counts of guides belonging to the same vector during in-memory import
#' @param grna_matrix Guide-by-cell count matrix with named guide rows.
#' @param grna_target_data_frame Metadata identifying rows by grna_id and their
#'   vector groups by vector_id.
#' @return Dense vector-by-cell matrix in first-seen vector order and original
#'   cell order. Repeated metadata rows select the same guide more than once.
#' @details Prepare sparse rows using Matrix's rules for duplicate numeric/logical
#'   entries. Sum bounded blocks of complete columns using dense colSums(), to
#'   preserve its platform-specific rounding. The vector-level output and transpose
#'   workspace are dense, but the entire sparse guide input is not made dense.
#'   Caller inputs remain unchanged.
#' @noRd
collapse_grna_matrix <- function(grna_matrix, grna_target_data_frame) {
    vector_ids <- grna_target_data_frame$vector_id |> unique()
    # Prepare rows once using Matrix's conversion, which combines duplicate
    # triplet entries. Here logical FALSE contributes zero to count sums;
    # do not use storage's separate rule that counts all stored positions.
    if (methods::is(grna_matrix, "sparseMatrix")) {
        grna_matrix <- methods::as(grna_matrix, "RsparseMatrix")
    } else {
        # Keep the original non-sparse coercion/dispatch behavior.
        grna_matrix <- as.matrix(grna_matrix)
    }
    n_cells <- ncol(grna_matrix)
    grna_matrix_collapsed <- vapply(
        vector_ids,
        function(curr_vector_id) {
            curr_grna_ids <- grna_target_data_frame |>
                dplyr::filter(vector_id == curr_vector_id) |>
                dplyr::pull(grna_id)
            selected <- grna_matrix[curr_grna_ids, , drop = FALSE]
            if (methods::is(selected, "sparseMatrix")) {
                sum_sparse_vector_rows(selected)
            } else {
                Matrix::colSums(selected)
            }
        },
        FUN.VALUE = numeric(n_cells)
    ) |>
        t()
    return(grna_matrix_collapsed)
}

#' Sum one vector's selected guide rows using dense column blocks
#' @param counts Sparse matrix with selected guides in rows and original cells in
#'   columns. Repeated metadata selections appear as repeated rows, in order.
#' @return Named column sums. Each temporary dense block targets at most 2^20
#'   values (8 MiB for doubles). A whole column is the minimum, so a vector with
#'   more guides than that can exceed the target. Output memory is additional.
#' @noRd
sum_sparse_vector_rows <- function(counts) {
    n_cells <- ncol(counts)
    # Most vectors contain only a few guides and already fit the block target.
    # Avoid another sparse column subset and numeric output copy in that case.
    if (as.double(nrow(counts)) * n_cells <= 2^20) {
        return(Matrix::colSums(as.matrix(counts)))
    }
    sums <- numeric(n_cells)
    names(sums) <- colnames(counts)
    if (n_cells == 0L) {
        return(sums)
    }

    # On Linux, dense colSums() can use extra precision while sparse sums use
    # doubles. Sum every selected guide for a cell together to keep the same
    # rounding as the original dense calculation. Split only between cells
    # (columns), not between the guides contributing to one cell.
    block_width <- max(1, floor(2^20 / max(1L, nrow(counts))))
    for (first in seq.int(1, n_cells, by = block_width)) {
        columns <- seq.int(first, min(n_cells, first + block_width - 1))
        sums[columns] <- Matrix::colSums(as.matrix(counts[,
            columns,
            drop = FALSE
        ]))
    }
    return(sums)
}
