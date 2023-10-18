#' Import data from cellranger
#'
#' `import_data_from_cellranger()` imports data from the output of one or more calls to cellranger count. \href{https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-whole_game_import_data}{manual}
#'
#' @param directories TBD
#' @param moi TBD
#' @param grna_target_data_frame TBD
#' @param extra_covariates TBD
#'
#' @return TBD
#' @export
#'
#' @examples
#' \dontrun{
#' # 1. create the sceptre object from CellRanger output
#' directories <- paste0(system.file("extdata", package = "sceptre"),
#'                      "/highmoi_example/gem_group_", 1:2)
#' data(grna_target_data_frame_highmoi)
#'
#' # 2. simulate an additional covariate
#' cell_type <- sample(x = paste0("type_", 1:3),
#'                     size = 45919,
#'                     replace = TRUE) |> factor()
#' extra_covariates <- data.frame(cell_type = cell_type)
#'
#' # 3. initialize the sceptre_object
#' sceptre_object <- import_data_from_cellranger(directories = directories,
#'                                               moi = "high",
#'                                               grna_target_data_frame = grna_target_data_frame_highmoi,
#'                                               extra_covariates = extra_covariates)
#'}
import_data_from_cellranger <- function(directories, moi, grna_target_data_frame, extra_covariates = data.frame()) {
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
      p <- obtain_pointer_vector(i = dt_sub$cell_idx, dim = matrix_metadata$n_cells)
      out <- Matrix::sparseMatrix(i = 1, p = c(0, 1), x = 1, repr = "C")
      out@i <- dt_sub$feature_idx
      out@p <- p
      out@x <- dt_sub$x
      out@Dim <- c(length(curr_feature_idxs), matrix_metadata$n_cells)
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


get_mtx_metadata <- function(mtx_file) {
  con <- file(mtx_file, "r")
  n_to_skip <- 0L
  while (TRUE) {
    n_to_skip <- n_to_skip + 1L
    line <- readLines(con, n = 1)
    if (substr(line, 0, 1) != "%") break()
  }
  close(con)
  metrics <- strsplit(line, split = " ", fixed = TRUE)[[1]] |> as.integer()
  list(n_features = metrics[1], n_cells = metrics[2], n_nonzero = metrics[3], n_to_skip = n_to_skip)
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
