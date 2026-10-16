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
#' sceptre_object <- import_data(
#'   response_matrix = lowmoi_example_data$response_matrix,
#'   grna_matrix = lowmoi_example_data$grna_matrix,
#'   grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'   extra_covariates = lowmoi_example_data$extra_covariates,
#'   moi = "low",
#'   use_ondisc = TRUE,
#'   directory_to_write = tempdir()
#' )
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
#' Import aligned in-memory counts and initialize the saved S4 facade
#'
#' @param response_matrix,grna_matrix Feature-by-source-cell count matrices with
#'   unique feature IDs; cell columns and any supplied cell names must align.
#' @param grna_target_data_frame Guide/target metadata; `vector_id`, when
#'   present,
#'   groups guide rows for count summation before covariates are computed.
#' @param moi Character `"low"` or `"high"`.
#' @param extra_covariates Data frame in source-cell order, or an empty data
#'   frame.
#' @param response_names Character feature labels aligned with response rows,
#'   or `NA_character_` to use response row names for mitochondrial
#'   identification.
#' @return Initialized `sceptre_object` with row-accessible sparse counts and
#'   source-ordered covariates. Validation/conversion can allocate new matrices.
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
        moi = moi,
        grna_target_data_frame = grna_target_data_frame
    )
}

#' Import R counts through the existing backed-store construction path
#'
#' @param response_matrix,grna_matrix,grna_target_data_frame,moi,extra_covariates,response_names
#'   Aligned import inputs as for `import_data_memory()`; guide rows are written
#'   as supplied, without `vector_id` aggregation in this path.
#' @param directory_to_write Single output directory, created if necessary.
#' @return Initialized object with backed matrices sharing a newly sampled
#'   `integer_id`. Writes `response.odm` and `grna.odm`, advancing the R RNG for
#'   the ID and following the existing ondisc construction/RNG order.
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
        moi = moi,
        grna_target_data_frame = grna_target_data_frame
    )

    # add integer id
    sceptre_object@integer_id <- integer_id
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
#' sceptre_object <- import_data_from_cellranger(
#'   directories = directories,
#'   moi = "high",
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   use_ondisc = TRUE,
#'   directory_to_write = tempdir()
#' )
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

#' Read Cell Ranger batches into aligned in-memory modalities
#'
#' @param directories Character paths in batch order, with matching feature
#'   tables and coordinate Matrix Market counts. Each modality follows its
#'   feature-file row order, even when modalities are interleaved; unique
#'   coordinates may be unsorted. Repeated coordinates fail before CSC use.
#' @param moi,grna_target_data_frame,extra_covariates Inputs to
#'   `import_data_memory()`;
#'   extra covariates follow concatenated batch/cell order.
#' @return Initialized object in input-directory/column order. Barcode files
#'   must
#'   be present for all batches or none; retained barcodes are batch-prefixed
#'   for
#'   multiple inputs. Adds/replaces the `batch` covariate and emits progress.
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
            # Feature-file rows need not form contiguous modality blocks. Map
            # into this modality's original row order, then use zero-based
            # slots.
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


#' Read Cell Ranger batches through ondisc and initialize the S4 facade
#'
#' @param directories Character paths passed in order to the ondisc importer.
#' @param moi,grna_target_data_frame,extra_covariates Import
#'   metadata/covariates;
#'   extra covariates must follow the ondisc importer's returned cell order.
#' @param directory_to_write Single output directory, created if necessary.
#' @return Object with backed gene/gRNA modalities and their shared import ID.
#'   ondisc writes backing files and computes cell covariates; `vector_id`
#'   aggregation is delegated to it and then reflected in guide metadata.
#'   Coordinate validation/duplicate behavior belongs to the external parser;
#'   the owned in-memory CSC policy does not qualify this route.
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

    sceptre_object@integer_id <- out$gene@integer_id
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

    # Parse files are cell-first triplets. Each read owns its table, converts
    # file positions to zero-based coordinates, then delegates complete stable
    # CSC ordering/duplicate rejection to create_csc_matrix. Output rows are
    # features and columns are cells; this is the in-memory import route.
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
#' Initialize persistent S4 slots without changing source-cell order
#'
#' @param response_matrix,grna_matrix Aligned feature-by-cell stores
#'   (row-accessible
#'   sparse matrices or ondisc objects), already validated by the import path.
#' @param covariate_data_frame Source-cell-ordered data frame with computed
#'   guide
#'   maximum/fraction columns, moved into `import_grna_assignment_info`.
#' @param moi Character `"low"` or `"high"`.
#' @param grna_target_data_frame Guide metadata; IDs/targets are coerced to
#'   character
#'   while row order and additional columns are preserved.
#' @return New `sceptre_object` with the legacy slot layout and initial workflow
#'   flags. Stores supplied matrices without reordering or eagerly reading them.
#' @noRd
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
    sceptre_object <- set_response_matrix(sceptre_object, response_matrix)
    sceptre_object <- set_grna_matrix(sceptre_object, grna_matrix)
    sceptre_object@covariate_data_frame <- covariate_data_frame
    sceptre_object@low_moi <- (moi == "low")
    sceptre_object@covariate_names <- sort(colnames(covariate_data_frame))
    sceptre_object@grna_target_data_frame <- grna_target_data_frame |>
        dplyr::mutate(
            grna_id = as.character(grna_id),
            grna_target = as.character(grna_target)
        )
    sceptre_object@elements_to_analyze <- NA_character_
    sceptre_object@nf_pipeline <- FALSE
    sceptre_object@nuclear <- FALSE
    sceptre_object@last_function_called <- "import_data"
    sceptre_object@functs_called <- workflow_initial_state()
    return(sceptre_object)
}


#' Construct complete CSC slots from a parsed file's coordinates
#'
#' `dt` contains zero-based integer `cell_idx` and `feature_idx`, and numeric
#' `x`. Sort owned output vectors by cell then feature; the input table is not
#' mutated. Keep every value, including explicit zeros, and reject repeated
#' file coordinates rather than selecting or summing them. Dimensions and slot
#' lengths must fit the native integer pointer builder before it is called.
#' Returns a numeric dgCMatrix in feature-file/cell-file order, without
#' dimnames.
#' No dense assay, RNG draws or file writes occur.
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


#' Read the first noncomment Matrix Market dimension line
#'
#' Returns three integer dimensions named in `col_id` order and `n_to_skip`,
#' counting comment/header/dimension lines. Retains the existing space-separated
#' parser for plain/gzip files. Owns its read connection and closes it on error
#' as well as success; this helper does not redesign malformed-file validation.
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


#' Write the two imported modalities under a common backing identity
#'
#' @param response_matrix,grna_matrix Aligned feature-by-cell count matrices.
#' @param directory_to_write Output directory accepted by
#'   `prepare_directory_to_write()`.
#' @param integer_id Integer identity shared by the saved object and both
#'   stores.
#' @return The `storage_create_backed_pair()` list of response/gRNA handles
#'   after
#'   writing `response.odm` and `grna.odm` in that order.
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


#' Create and canonicalize an ondisc output directory
#'
#' @param directory_to_write Nonmissing scalar character path; a leading home
#'   abbreviation is expanded before recursive directory creation.
#' @return Existing absolute path using forward slashes, with embedded tildes
#'   rejected before it reaches ondisc's legacy path expansion.
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


#' Canonicalize an existing object or backing-file path
#'
#' @param file_path Nonmissing scalar character path that must already exist.
#' @param arg_name Character argument label used in validation errors.
#' @return Absolute path using forward slashes, with home expansion applied and
#'   remaining tildes rejected; does not open or modify the file.
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


#' Reject embedded tildes before invoking ondisc's legacy path handling
#'
#' @param file_path Canonical character path with leading home expansion
#'   complete.
#' @param arg_name Character argument label used in the error.
#' @return `NULL` on success; errors if any tilde remains. Does not access
#'   files.
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


#' Compute aligned cell covariates and prepare row-accessible R matrices
#'
#' @param response_matrix,grna_matrix In-memory feature-by-cell count matrices.
#' @param extra_covariates Empty data frame or covariates in the same cell
#'   order.
#' @param response_names Aligned response labels, or `NA_character_` for row
#'   names.
#' @return List with `covariate_data_frame`, `response_matrix`, and
#'   `grna_matrix`.
#'   Dimensions, feature IDs and cell order are retained; count storage may be
#'   converted. This runs before the backed import writes its matrices to disk.
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
        covariate_data_frame = covariate_data_frame,
        response_matrix = response_matrix,
        grna_matrix = grna_matrix
    ))
}


#' Replace guide metadata rows by their vector identities
#'
#' @param grna_target_data_frame Data frame with `grna_id`, `vector_id`, and
#'   `grna_target`; every vector must map to exactly one target.
#' @return Metadata with `vector_id` renamed to `grna_id` and original guide IDs
#'   removed. Deduplicates complete remaining rows in first-seen order, so
#'   differing
#'   additional annotations can retain repeated vector IDs.
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


#' Sum guide counts within each vector for the in-memory import path
#'
#' @param grna_matrix Guide-by-cell count matrix with named guide rows.
#' @param grna_target_data_frame Metadata identifying matrix rows by `grna_id`
#'   and their vector groups by `vector_id`.
#' @return Dense numeric vector-by-cell matrix in first-seen vector order and
#'   original cell order. Metadata repetitions contribute repeated row
#'   selections.
#'   Sparse row preparation retains Matrix's numeric/logical duplicate
#'   semantics.
#'   Sums use bounded dense column blocks to preserve the original platform's
#'   dense accumulation. Vector-by-cell assembly/output and its existing
#'   transpose
#'   workspace are dense; the full sparse input is not. Inputs remain unchanged;
#'   there are no RNG or file effects.
#' @noRd
collapse_grna_matrix <- function(grna_matrix, grna_target_data_frame) {
    vector_ids <- grna_target_data_frame$vector_id |> unique()
    # Prepare row access once, including valid triplet duplicate normalization.
    # Matrix keeps logical FALSE at zero here; the storage adapter's distinct
    # stored-position conversion is not the aggregation count contract.
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

#' Sum a selected sparse guide group with the original dense accumulator
#'
#' @param counts Sparse matrix whose rows already include metadata repetitions
#'   in their original order; columns retain source-cell order.
#' @return Named numeric column sums. Each temporary dense block targets at most
#'   2^20 elements (8 MiB for doubles); one complete column is the minimum, so a
#'   larger guide group can exceed that target. Output storage is additional.
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

    # Linux dense colSums can accumulate in extended precision while Matrix's
    # sparse method uses doubles. Keep every selected row in a single sum;
    # adding
    # partial row sums would change rounding. Column blocks only bound
    # workspace.
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
