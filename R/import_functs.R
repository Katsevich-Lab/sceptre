#' Import data from cellranger
#'
#' `import_data_from_cellranger()` imports data from the output of one or more calls to cellranger count. See \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-whole_game_import_data}{Section 1 of the introductory chapter in the manual} for more information about this function.
#'
#' @param directories a character vector of directories containing the output of one or more calls to cellranger count. Each directory should contain the files "matrix.mtx.gz" and "features.tsv.gz" (and optionally "barcodes.tsv.gz").
#' @param moi a string indicating the MOI of the dataset, either "low" or "high"
#' @param grna_target_data_frame a data frame containing columns `grna_id` and `grna_target` mapping each individual gRNA to its target
#' @param extra_covariates (optional) a data frame containing extra covariates (e.g., batch, biological replicate) beyond those that `sceptre` can compute
#'
#' @return an initialized `sceptre_object`
#' @export
#'
#' @examples
#' \dontrun{
#' # 1. point to directories containing cellranger output
#' directories <- paste0(system.file("extdata", package = "sceptre"),
#'                      "/highmoi_example/gem_group_", 1:2)
#'
#' # 2. simulate an additional covariate
#' cell_type <- sample(x = paste0("type_", 1:3),
#'                     size = 45919, replace = TRUE) |> factor()
#' extra_covariates <- data.frame(cell_type = cell_type)
#'
#' # 3. initialize the sceptre_object
#' data(grna_target_data_frame_highmoi)
#' sceptre_object <- import_data_from_cellranger(
#'   directories = directories,
#'   moi = "high",
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   extra_covariates = extra_covariates
#' )
#' }
import_data_from_cellranger <- function(directories, moi, grna_target_data_frame, extra_covariates = data.frame(), use_ondisc = FALSE, directory_to_write = NULL) {
  # take cases on use_ondisc
  if (!use_ondisc) {
    import_data_from_cellranger_memory(directories, moi, grna_target_data_frame, extra_covariates)
  } else {
    import_data_from_cellranger_disk(directories, moi, grna_target_data_frame, extra_covariates, directory_to_write)
  }
}


import_data_from_cellranger_memory <- function(directories, moi, grna_target_data_frame, extra_covariates) {
  # 1. check that directories exist
  for (curr_directory in directories) {
    if (!dir.exists(curr_directory)) stop(paste0("The directory ", curr_directory, " does not exist."))
  }

  # 2. create the vector of matrix, features, and barcode files
  input_files <- sapply(X = directories, function(curr_directory) {
    fs <- list.files(curr_directory)
    grep_strs <- c("*features.tsv($|.gz)", "*matrix.mtx($|.gz)")
    out <- sapply(grep_strs, function(grep_str) {
      file_names <- grep(pattern = grep_str, x = fs, value = TRUE)
      if (length(file_names) >= 2L) {
        stop(paste0("There are multiple ", grep_str, " files within the directory ", curr_directory, "."))
      }
      if (length(file_names) == 0L) {
        stop(paste0("The directory ", curr_directory, " contains zero ", grep_str, " files."))
      }
      return(paste0(curr_directory, "/", file_names))
    }) |> stats::setNames(c("features", "matrix"))
  }, simplify = FALSE) |> stats::setNames(NULL)

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
  sceptre_object <- import_data(response_matrix = out_mats[["Gene Expression"]],
                                grna_matrix = out_mats[["CRISPR Guide Capture"]],
                                grna_target_data_frame = grna_target_data_frame,
                                moi = moi,
                                extra_covariates = extra_covariates,
                                response_names = gene_names)
  gc() |> invisible()
  cat(crayon::green(' \u2713'))
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


write_sceptre_object_to_cellranger_format <- function(sceptre_object, directory) {
  # 0. create directory
  if (dir.exists(directory)) unlink(directory, recursive = TRUE)
  dir.create(directory)

  # 1. combine matrices across modalities
  response_matrix <- sceptre_object@response_matrix |> set_matrix_accessibility(make_row_accessible = FALSE)
  grna_matrix <- sceptre_object@grna_matrix |> set_matrix_accessibility(make_row_accessible = FALSE)
  combined_mat <- rbind(response_matrix, grna_matrix)

  # 2. construct the features df
  response_ids <- rownames(response_matrix)
  response_names <- sceptre_object@response_names
  grna_ids <- rownames(grna_matrix)
  feature_df <- data.frame(response_id = c(response_ids, grna_ids),
                           response_name = c(response_names, grna_ids),
                           modality = c(rep("Gene Expression", length(response_ids)),
                                        rep("CRISPR Guide Capture", length(grna_ids))))

  # 3. split the matrices according to batch; loop over the batches and save the matrix and features file
  batch_v <- sceptre_object@covariate_data_frame$batch
  batch_levels_v <- as.character(unique(batch_v))
  for (i in seq_along(batch_levels_v)) {
    batch_level <- batch_levels_v[i]
    mat_sub <- combined_mat[feature_df$response_id, batch_level == batch_v]
    dir_name <- paste0(directory, "/gem_group_", i)
    dir.create(dir_name)
    Matrix::writeMM(obj = mat_sub, file = paste0(dir_name, "/matrix.mtx"))
    readr::write_tsv(x = feature_df, file = paste0(dir_name, "/features.tsv"), col_names = FALSE)
    readr::write_tsv(x = data.frame(), file = paste0(dir_name, "/barcodes.tsv"))
    curr_files <- list.files(dir_name, full.names = TRUE)
    for (curr_file in curr_files) {
      R.utils::gzip(filename = curr_file, destname = paste0(curr_file, ".gz"))
    }
  }
  return(NULL)
}


#' Import data from Parse (experimental)
#'
#' `import_data_from_parse()` imports data from the output of the Parse count matrix generation program. See \href{https://timothy-barry.github.io/sceptre-book/import-data.html#import-from-the-parse-program-experimental}{Chapter 1 of the manual} for more information about this function. It is assumed that the data are stored in a single set of files (as opposed to multiple sets of files corresponding to, e.g., different samples).
#'
#' `import_data_from_parse()` is experimental, and the API of this function is subject to change. We expect the API to solidify as we learn more about the Parse platform and the structure of the Parse count matrix generation program output.
#'
#' @param gene_mat_fp file path to the gene count_matrix.mtx file
#' @param grna_mat_fp file path to the gRNA count_matrix.mtx file
#' @param all_genes_fp file path to the all_genes.csv file containing the gene IDs
#' @param all_grnas_fp file path to the all_guides.csv file containing the gRNA IDs
#' @param moi a string indicating the MOI of the dataset, either "low" or "high"
#' @param grna_target_data_frame a data frame containing columns `grna_id` and `grna_target` mapping each individual gRNA to its target
#' @param extra_covariates (optional) a data frame containing extra covariates (e.g., batch, biological replicate) beyond those that `sceptre` can compute
#'
#' @return an initialized `sceptre_object`
#' @export
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
  sceptre_object <- import_data(response_matrix = gene_mat,
                                grna_matrix = grna_mat,
                                grna_target_data_frame = grna_target_data_frame,
                                moi = moi,
                                extra_covariates = extra_covariates,
                                response_names = gene_names)
  return(sceptre_object)
}
