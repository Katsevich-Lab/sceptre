############
# R MATRICES
############
#' Import data
#'
#' `import_data()` imports data from a collection of R objects to create a `sceptre_object`. See \href{https://timothy-barry.github.io/sceptre-book/import-data.html#import-data-from-a-collection-of-r-objects}{Chapter 1 of the manual} for detailed information about this function.
#'
#' @param response_matrix a matrix of response UMI counts, with responses in rows and cells in columns
#' @param grna_matrix a matrix of gRNA UMI counts, with gRNAs in rows and cells in columns
#' @param grna_target_data_frame a data frame containing columns `grna_id` and `grna_target` mapping each individual gRNA to its target
#' @param moi a string indicating the MOI of the dataset, either "low" or "high"
#' @param extra_covariates (optional) a data frame containing extra covariates (e.g., batch, biological replicate) beyond those that `sceptre` can compute
#' @param response_names (optional) a vector of human-readable response names; names with the prefix "MT-" are assumed to be mitochondrial genes and are used to compute the covariate `response_p_mito`.
#' @param use_ondisc (default `FALSE`) a logical value (i.e., `TRUE` or `FALSE`) indicating whether to create an `ondisc`-backed `sceptre_object` (`TRUE`) or a standard `sceptre_object` (`FALSE`)
#' @param directory_to_write (optional) a file path to a directory in which to write the backing response and gRNA expression matrices. Must be supplied if `use_ondisc` is set to `TRUE`.
#'
#' @return an initialized `sceptre_object`
#' @export
#' @examples
#' # see example via ?sceptre
import_data <- function(response_matrix, grna_matrix, grna_target_data_frame, moi, extra_covariates = data.frame(),
                        response_names = NA_character_, use_ondisc = FALSE, directory_to_write = NULL) {
  if (use_ondisc) { # use ondisc
    import_data_use_ondisc(response_matrix, grna_matrix, grna_target_data_frame, moi,
                           extra_covariates, response_names, directory_to_write)
  } else { # do not use ondisc
    import_data_memory(response_matrix, grna_matrix, grna_target_data_frame, moi, extra_covariates, response_names)
  }
}

# memory
import_data_memory <- function(response_matrix, grna_matrix, grna_target_data_frame,
                               moi, extra_covariates, response_names) {
  #  perform initial check
  check_import_data_inputs(response_matrix, grna_matrix, grna_target_data_frame, moi, extra_covariates) |> invisible()

  # collapse gRNA matrix (if vector_id supplied)
  if ("vector_id" %in% colnames(grna_target_data_frame)) {
    grna_matrix <- collapse_grna_matrix(grna_matrix = grna_matrix,
                                        grna_target_data_frame = grna_target_data_frame)
    grna_target_data_frame <- collapse_grna_target_data_frame(grna_target_data_frame)
  }

  # process covariates and matrices
  processed_inputs <- process_covariates_and_matrices(response_matrix = response_matrix,
                                                      grna_matrix = grna_matrix,
                                                      extra_covariates = extra_covariates,
                                                      response_names = response_names)
  # initialize the sceptre_object
  sceptre_object <- init_sceptre_object(response_matrix = processed_inputs$response_matrix,
                                        grna_matrix = processed_inputs$grna_matrix,
                                        covariate_data_frame = processed_inputs$covariate_data_frame,
                                        moi = moi,
                                        grna_target_data_frame = grna_target_data_frame)
}

# 1. ondisc
import_data_use_ondisc <- function(response_matrix, grna_matrix, grna_target_data_frame, moi,
                                   extra_covariates, response_names, directory_to_write) {
  #  perform initial check
  check_import_data_inputs(response_matrix, grna_matrix, grna_target_data_frame, moi, extra_covariates) |> invisible()
  # check directory_to_write
  if (is.null(directory_to_write)) stop("`directory_to_write` cannot be `NULL`.")

  # process covariates and matrices
  processed_inputs <- process_covariates_and_matrices(response_matrix = response_matrix,
                                                      grna_matrix = grna_matrix,
                                                      extra_covariates = extra_covariates,
                                                      response_names = response_names)
  response_matrix <- processed_inputs$response_matrix
  grna_matrix <- processed_inputs$grna_matrix
  covariate_data_frame <- processed_inputs$covariate_data_frame

  # write matrices to disk
  integer_id <- sample(x = seq(0L, .Machine$integer.max), size = 1L)
  out <- write_matrices_to_disk(response_matrix = response_matrix,
                                grna_matrix = grna_matrix,
                                directory_to_write = directory_to_write,
                                integer_id = integer_id)

  # initialize the sceptre_object
  sceptre_object <- init_sceptre_object(response_matrix = out$response_matrix,
                                        grna_matrix = out$grna_matrix,
                                        covariate_data_frame = covariate_data_frame,
                                        moi = moi,
                                        grna_target_data_frame = grna_target_data_frame)

  # add integer id
  sceptre_object@integer_id <- integer_id
  return(sceptre_object)
}


############
# CELLRANGER
############
#' Import data from cellranger
#'
#' `import_data_from_cellranger()` imports data from the output of one or more calls to cellranger count. See \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-whole_game_import_data}{Section 1 of the introductory chapter in the manual} for more information about this function.
#'
#' @param directories a character vector of directories containing the output of one or more calls to cellranger count. Each directory should contain the files "matrix.mtx.gz" and "features.tsv.gz" (and optionally "barcodes.tsv.gz").
#' @param moi a string indicating the MOI of the dataset, either "low" or "high"
#' @param grna_target_data_frame a data frame containing columns `grna_id` and `grna_target` mapping each individual gRNA to its target. Optionally, `grna_target_data_frame` can contain columns `chr`, `start`, and `end`, giving the chromosome, start coordinate, and end coordiante, respectively, of each gRNA. Additionally, `grna_target_data_frame` can contain the column `vector_id` specifying the vector to which a given gRNA belongs.
#' @param extra_covariates (optional) a data frame containing extra covariates (e.g., batch, biological replicate) beyond those that `sceptre` can compute
#' @param use_ondisc (optional; default `FALSE`) a logical indicating whether to store the expression data in a disk-backed `ondisc` matrix (TRUE) or an in-memory sparse matrix (FALSE)
#' @param directory_to_write (optional) a string indicating the directory in which to write the backing `.odm` files (must be specified if `use_ondisc` is set to `TRUE`)
#'
#' @return an initialized `sceptre_object`
#' @export
#' @examples
#' library(sceptredata)
#' # 1. point to directories containing cellranger output
#' directories <- paste0(system.file("extdata", package = "sceptredata"),
#'                      "/highmoi_example/gem_group_", 1:2)
#'
#' # 2. simulate an additional covariate
#' cell_type <- sample(x = paste0("type_", 1:3),
#'                     size = 45919, replace = TRUE) |> factor()
#' extra_covariates <- data.frame(cell_type = cell_type)
#'
#' # 3. initialize the sceptre_object
#' #' library(sceptredata)
#' data(grna_target_data_frame_highmoi)
#' sceptre_object <- import_data_from_cellranger(
#'   directories = directories,
#'   moi = "high",
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   extra_covariates = extra_covariates
#' )
import_data_from_cellranger <- function(directories, moi, grna_target_data_frame, extra_covariates = data.frame(), use_ondisc = FALSE, directory_to_write = NULL) {
  # take cases on use_ondisc
  if (!use_ondisc) {
    import_data_from_cellranger_memory(directories, moi, grna_target_data_frame, extra_covariates)
  } else {
    import_data_from_cellranger_use_ondisc(directories, moi, grna_target_data_frame, extra_covariates, directory_to_write)
  }
}

# 1. memory
import_data_from_cellranger_memory <- function(directories, moi, grna_target_data_frame, extra_covariates) {
  # 1. check that directories exist
  for (curr_directory in directories) {
    if (!dir.exists(curr_directory)) stop(paste0("The directory ", curr_directory, " does not exist."))
  }

  # 2. create the vector of matrix, features, and barcode files
  input_files <- lapply(X = directories, function(curr_directory) {
    fs <- list.files(curr_directory)
    grep_strs <- c("*features.tsv($|.gz)", "*matrix.mtx($|.gz)")
    out <- vapply(grep_strs, function(grep_str) {
      file_names <- grep(pattern = grep_str, x = fs, value = TRUE)
      if (length(file_names) >= 2L) {
        stop(paste0("There are multiple ", grep_str, " files within the directory ", curr_directory, "."))
      }
      if (length(file_names) == 0L) {
        stop(paste0("The directory ", curr_directory, " contains zero ", grep_str, " files."))
      }
      return(paste0(curr_directory, "/", file_names))
    }, FUN.VALUE = character(1)) |> stats::setNames(c("features", "matrix"))
  }) |> stats::setNames(NULL)

  # 3. obtain the barcodes, features, and matrix fp vectors
  feature_fps <- sapply(X = input_files, FUN = function(i) i[["features"]])
  matrix_fps <- sapply(X = input_files, FUN = function(i) i[["matrix"]])

  # 4. obtain the feature data frame; determine the ids of each feature
  feature_df <- data.table::fread(file = feature_fps[1],
                                  colClasses = c("character", "character", "character"),
                                  col.names = c("feature_id", "feature_name", "modality"), header = FALSE)
  if (!("CRISPR Guide Capture" %in% feature_df$modality)) {
    stop("The features.tsv file should contain the modality `CRISPR Guide Capture`.")
  }
  modalities <- unique(feature_df$modality)
  modalities <- modalities |> stats::setNames(modalities)
  modality_idxs <- lapply(modalities, function(modality) which(feature_df$modality == modality))

  # 5. loop over directories, loading the matrix corresponding to each
  out_list <- vector(mode = "list", length = length(directories))
  n_cells_per_matrix <- vector(mode = "integer", length = length(directories))
  for (i in seq_along(directories)) {
    cat(paste0("Processing directory ", i, "."))
    # check that the features df matches
    features_fp <- feature_fps[i]
    curr_feature_df <- data.table::fread(file = features_fp,
                                         colClasses = c("character", "character", "character"),
                                         col.names = c("feature_id", "feature_name", "modality"), header = FALSE)
    for (col in colnames(feature_df)) {
      if (any(dplyr::pull(feature_df, dplyr::all_of(col)) !=
              dplyr::pull(curr_feature_df, dplyr::all_of(col)))) {
        stop("The features.tsv files must match across directories.")
      }
    }

    # process the matrix
    matrix_fp <- matrix_fps[i]
    matrix_metadata <- get_mtx_metadata(matrix_fp)
    n_cells_per_matrix[i] <- matrix_metadata$n_cells
    dt <- data.table::fread(file = matrix_fp, skip = matrix_metadata$n_to_skip,
                            col.names = c("feature_idx", "cell_idx", "x"),
                            colClasses = c("integer", "integer", "double"),
                            showProgress = FALSE, nThread = 2L)
    subsetted_mats <- lapply(modalities, function(modality) {
      curr_feature_idxs <- modality_idxs[[modality]]
      # subset
      dt_sub <- dt[dt$feature_idx %in% curr_feature_idxs,]
      # offset the feature idxs such that they start at 1
      increment_feature_vector_value <- -min(curr_feature_idxs)
      increment_vector(dt_sub$feature_idx, increment_feature_vector_value)
      increment_vector(dt_sub$cell_idx, -1L)
      # create a sparse matrix in csc format
      out <- create_csc_matrix(dt = dt_sub, n_cells = matrix_metadata$n_cells,
                               n_features = length(curr_feature_idxs))
      return(out)
    })
    rm(dt); gc() |> invisible()
    out_list[[i]] <- subsetted_mats
    cat(crayon::green(' \u2713\n'))
  }

  # 6. combine the matrices via do.call cbind
  cat("Combining matrices across directories.")
  out_mats <- lapply(modalities, function(modality) {
    l <- lapply(out_list, function(mat) mat[[modality]])
    combined_mat <- do.call(what = "cbind", args = l)
    return(combined_mat)
  })
  rm(out_list); gc() |> invisible()
  cat(crayon::green(' \u2713\n'))

  # 7. collect metadata pieces
  for (modality in modalities) {
    rownames(out_mats[[modality]]) <- feature_df$feature_id[modality_idxs[[modality]]]
  }
  gene_names <- feature_df$feature_name[modality_idxs[["Gene Expression"]]]

  # 8. set up the extra covariates
  batch <- if (length(directories) >= 2L) {
    rep(paste0("b", seq_along(n_cells_per_matrix)), n_cells_per_matrix) |> factor()
  } else {
    NULL
  }
  if (nrow(extra_covariates) == 0L) {
    extra_covariates <- data.frame(batch = batch)
  } else {
    extra_covariates$batch <- batch
  }

  # 9. initialize the sceptre object
  cat("Creating the sceptre object.")
  sceptre_object <- import_data_memory(response_matrix = out_mats[["Gene Expression"]],
                                       grna_matrix = out_mats[["CRISPR Guide Capture"]],
                                       grna_target_data_frame = grna_target_data_frame,
                                       moi = moi,
                                       extra_covariates = extra_covariates,
                                       response_names = gene_names)
  gc() |> invisible()
  cat(crayon::green(' \u2713'))
  return(sceptre_object)
}


# 2. disk
import_data_from_cellranger_use_ondisc <- function(directories, moi, grna_target_data_frame, extra_covariates, directory_to_write) {
  # 0. check that directory_to_write has been supplied
  if (is.null(directory_to_write)) stop("`directory_to_write` must be supplied.")

  # 1. call the corresponding ondisc function
  vector_supplied <- "vector_id" %in% colnames(grna_target_data_frame)
  out <- ondisc::create_odm_from_cellranger(directories_to_load = directories,
                                            directory_to_write = directory_to_write,
                                            write_cellwise_covariates = FALSE,
                                            grna_target_data_frame = if (vector_supplied) grna_target_data_frame else NULL)

  # 1.5. update grna_target_data_frame, if vector supplied
  if (vector_supplied) {
    grna_target_data_frame <- collapse_grna_target_data_frame(grna_target_data_frame)
  }

  # 2. check data imports
  check_import_data_inputs(out$gene, out$grna, grna_target_data_frame, moi, extra_covariates) |> invisible()

  # 3. process the cellwise covariates
  covariate_df <- out$cellwise_covariates
  colnames(covariate_df) <- gsub(pattern = "gene", replacement = "response", fixed = TRUE, x = colnames(covariate_df))
  if (nrow(extra_covariates) >= 1L) covariate_df <- cbind(extra_covariates, covariate_df)

  # 4. initialize the sceptre_object
  sceptre_object <- init_sceptre_object(response_matrix = out$gene,
                                        grna_matrix = out$grna,
                                        covariate_data_frame = covariate_df,
                                        moi = moi,
                                        grna_target_data_frame = grna_target_data_frame)

  # 5. add integer id
  sceptre_object@integer_id <- out$gene@integer_id
  return(sceptre_object)
}


#######
# PARSE
#######
#' Import data from Parse (experimental)
#'
#' `import_data_from_parse()` imports data from the output of the Parse count matrix generation program. See \href{https://timothy-barry.github.io/sceptre-book/import-data.html#import-from-the-parse-program-experimental}{Chapter 1 of the manual} for more information about this function. It is assumed that the data are stored in a single set of files (as opposed to multiple sets of files corresponding to, e.g., different samples).
#'
#' `import_data_from_parse()` is experimental, and the API of this function is subject to change. We expect the API to solidify as we learn more about the Parse platform and the structure of the Parse count matrix generation program output.
#'
#' @param gene_mat_fp file path to the gene count_matrix.mtx file
#' @param grna_mat_fp file path to the gRNA count_matrix.mtx file
#' @param all_genes_fp file path to the all_genes.csv file containing the gene IDs
#' @param all_grnas_fp file path to the all_guides.csv file containing the gRNA IDs. The gRNA IDs are assumed to be in the second column (i.e., the "gene_name" column) of this file.
#' @param moi a string indicating the MOI of the dataset, either "low" or "high"
#' @param grna_target_data_frame a data frame containing columns `grna_id` and `grna_target` mapping each individual gRNA to its target
#' @param extra_covariates (optional) a data frame containing extra covariates (e.g., batch, biological replicate) beyond those that `sceptre` can compute
#'
#' @return an initialized `sceptre_object`
#' @export
#'
#' @examples
#' directory <- paste0(system.file("extdata", package = "sceptredata"), "/parse_example/")
#' gene_mat_fp <- paste0(directory, "gene_mat.mtx")
#' grna_mat_fp <- paste0(directory, "grna_mat.mtx")
#' all_genes_fp <- paste0(directory, "all_genes.csv")
#' all_grnas_fp <- paste0(directory, "all_grnas.csv")
#' grna_target_data_frame <- data.frame(
#'    grna_id = c("guide_A", "guide_B", "guide_C"),
#'    grna_target = c("target-A", "target-B", "non-targeting")
#' )
#' sceptre_object <- import_data_from_parse(
#'   gene_mat_fp = gene_mat_fp,
#'   grna_mat_fp = grna_mat_fp,
#'   all_genes_fp = all_genes_fp,
#'   all_grnas_fp = all_grnas_fp,
#'   moi = "low",
#'   grna_target_data_frame = grna_target_data_frame
#' )
import_data_from_parse <- function(gene_mat_fp, grna_mat_fp, all_genes_fp, all_grnas_fp,
                                   moi, grna_target_data_frame, extra_covariates = data.frame()) {
  # 1. check that files exist
  fps <- c(gene_mat_fp = gene_mat_fp, grna_mat_fp = grna_mat_fp,
           all_genes_fp = all_genes_fp, all_grnas_fp = all_grnas_fp)
  for (fp_name in names(fps)) {
    fp <- fps[[fp_name]]
    if (!file.exists(fp)) stop(paste0("File ", fp, " does not exist. Set `", fp_name, "` to a different value."))
  }

  # 2. load the gene.mtx and grna.mtx data
  load_matrix <- function(mat_fp) {
    matrix_metadata <- get_mtx_metadata(mat_fp, c("n_cells", "n_features", "n_nonzero"))
    dt <- data.table::fread(file = mat_fp, skip = matrix_metadata$n_to_skip,
                            col.names = c("cell_idx", "feature_idx", "x"),
                            colClasses = c("integer", "integer", "double"),
                            showProgress = FALSE, nThread = 2L)
    # make the cell and feature idxs zero-based
    increment_vector(dt$cell_idx, value = -1L)
    increment_vector(dt$feature_idx, value = -1L)

    # create a sparse matrix in csc format
    out <- create_csc_matrix(dt = dt, n_cells = matrix_metadata$n_cells,
                             n_features = matrix_metadata$n_features)
    return(out)
  }
  gene_mat <- load_matrix(gene_mat_fp)
  grna_mat <- load_matrix(grna_mat_fp)

  # 3. load the gene ids; update rownames of gene mat
  gene_feature_df <- data.table::fread(file = all_genes_fp,
                                       colClasses = c("character", "character", "character"),
                                       col.names = c("feature_id", "feature_name", "genome"), header = TRUE)
  rownames(gene_mat) <- gene_feature_df$feature_id[seq(1L, nrow(gene_mat))]
  gene_names <- gene_feature_df$feature_name[seq(1L, nrow(gene_mat))]

  # 4. load the grna ids; update rownames of grna mat
  grna_feature_df <- data.table::fread(file = all_grnas_fp,
                                       colClasses = c("character", "character", "character"),
                                       col.names = c("feature_id", "feature_name", "genome"), header = TRUE)
  rownames(grna_mat) <- grna_feature_df$feature_name

  # 5. create the sceptre_object
  sceptre_object <- import_data_memory(response_matrix = gene_mat,
                                       grna_matrix = grna_mat,
                                       grna_target_data_frame = grna_target_data_frame,
                                       moi = moi,
                                       extra_covariates = extra_covariates,
                                       response_names = gene_names)
  return(sceptre_object)
}


##################
# HELPER FUNCTIONS
##################
init_sceptre_object <- function(response_matrix, grna_matrix, covariate_data_frame, moi, grna_target_data_frame) {
  # create sceptre_object
  sceptre_object <- methods::new("sceptre_object")

  # process the covariate data frame
  sceptre_object@import_grna_assignment_info <- list(max_grna = covariate_data_frame$grna_feature_w_max_expression,
                                                     max_grna_frac_umis = covariate_data_frame$grna_frac_umis_max_feature)
  covariate_data_frame$grna_feature_w_max_expression <- covariate_data_frame$grna_frac_umis_max_feature <- NULL

  # set data fields
  sceptre_object <- set_response_matrix(sceptre_object, response_matrix)
  sceptre_object <- set_grna_matrix(sceptre_object, grna_matrix)
  sceptre_object@covariate_data_frame <- covariate_data_frame
  sceptre_object@low_moi <- (moi == "low")
  sceptre_object@covariate_names <- sort(colnames(covariate_data_frame))
  sceptre_object@grna_target_data_frame <- grna_target_data_frame |> dplyr::mutate(grna_id = as.character(grna_id),
                                                                                   grna_target = as.character(grna_target))
  # 5. initialize flags
  sceptre_object@elements_to_analyze <- NA_character_
  sceptre_object@nf_pipeline <- FALSE
  sceptre_object@nuclear <- FALSE
  sceptre_object@last_function_called <- "import_data"
  sceptre_object@functs_called <- c(import_data = TRUE, set_analysis_parameters = FALSE,
                                    assign_grnas = FALSE, run_qc = FALSE, run_calibration_check = FALSE,
                                    run_power_check = FALSE, run_discovery_analysis = FALSE)
  return(sceptre_object)
}


create_csc_matrix <- function(dt, n_cells, n_features) {
  p <- obtain_pointer_vector(i = dt$cell_idx, dim = n_cells)
  out <- Matrix::sparseMatrix(i = 1, p = c(0, 1), x = 1, repr = "C")
  out@i <- dt$feature_idx
  out@p <- p
  out@x <- dt$x
  out@Dim <- c(n_features, n_cells)
  return(out)
}


get_mtx_metadata <- function(mtx_file, col_id = c("n_features", "n_cells", "n_nonzero")) {
  con <- file(mtx_file, "r")
  n_to_skip <- 0L
  while (TRUE) {
    n_to_skip <- n_to_skip + 1L
    line <- readLines(con, n = 1)
    if (substr(line, 0, 1) != "%") break()
  }
  close(con)
  metrics <- strsplit(line, split = " ", fixed = TRUE)[[1]] |> as.integer()
  out <- list(metrics[1], metrics[2], metrics[3], n_to_skip)
  names(out) <- c(col_id, "n_to_skip")
  return(out)
}


write_matrices_to_disk <- function(response_matrix, grna_matrix, directory_to_write, integer_id) {
  response_matrix <- ondisc::create_odm_from_r_matrix(mat = response_matrix,
                                                      file_to_write = paste0(directory_to_write, "/response.odm"),
                                                      integer_id = integer_id)
  grna_matrix <- ondisc::create_odm_from_r_matrix(mat = grna_matrix,
                                                  file_to_write = paste0(directory_to_write, "/grna.odm"),
                                                  integer_id = integer_id)
  return(list(response_matrix = response_matrix, grna_matrix = grna_matrix))
}


process_covariates_and_matrices <- function(response_matrix, grna_matrix, extra_covariates, response_names) {
  covariate_data_frame <- auto_compute_cell_covariates(response_matrix = response_matrix,
                                                       grna_matrix = grna_matrix,
                                                       extra_covariates = extra_covariates,
                                                       response_names = if (identical(response_names, NA_character_)) rownames(response_matrix) else response_names)
  response_matrix <- set_matrix_accessibility(response_matrix, make_row_accessible = TRUE)
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
  return(list(covariate_data_frame = covariate_data_frame, response_matrix = response_matrix, grna_matrix = grna_matrix))
}


collapse_grna_target_data_frame <- function(grna_target_data_frame) {
  vector_nested_within_target <- (grna_target_data_frame |>
    dplyr::group_by(vector_id) |>
    dplyr::summarize(n_distinct_targets = length(unique(grna_target))) |>
    dplyr::pull(n_distinct_targets) == 1) |> all()
  if (!vector_nested_within_target) stop("`vector_id` must be nested within `grna_target`.")

  grna_target_data_frame |>
    dplyr::select(-grna_id) |>
    dplyr::rename(grna_id = vector_id) |>
    dplyr::distinct()
}


collapse_grna_matrix <- function(grna_matrix, grna_target_data_frame) {
  # temporary function; we may rewrite this function to improve its
  # speed and memory efficiency if multiguide vector data become more prevalent.
  vector_ids <- grna_target_data_frame$vector_id |> unique()
  suppressWarnings(grna_matrix <- as.matrix(grna_matrix))
  grna_matrix_collapsed <- sapply(vector_ids, function(curr_vector_id) {
    curr_grna_ids <- grna_target_data_frame |>
      dplyr::filter(vector_id == curr_vector_id) |>
      dplyr::pull(grna_id)
    Matrix::colSums(grna_matrix[curr_grna_ids,,drop=FALSE])
  }) |> t()
  return(grna_matrix_collapsed)
}
