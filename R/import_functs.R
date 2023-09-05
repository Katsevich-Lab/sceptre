#' Import data from cellranger
#'
#' @param directories a vector of directories containing the count matrices
#'
#' @return a sceptre_object initialized
#' @export
#'
#' @examples
#' base_dir <- "/Users/timbarry/research_offsite/external/replogle-2022/raw/test_data/"
#' # load the grna group df
#' grna_group_data_frame <- readRDS(paste0(base_dir, "grna_group_df.rds"))
#' # set paths to the directories containing the cellranger count outputs
#' directories <- paste0(base_dir, "gemgroup", 1:3)
#' sceptre_object <- import_data_from_cellranger(directories, "low", grna_group_data_frame)
import_data_from_cellranger <- function(directories, moi, grna_group_data_frame) {
  # 1. check that directories exist
  for (curr_directory in directories) {
    if (!dir.exists(curr_directory)) stop(paste0("The directory ", curr_directory, " does not exist."))
  }

  # 2. create the vector of matrix, features, and barcode files
  input_files <- sapply(X = directories, function(curr_directory) {
    fs <- list.files(curr_directory)
    grep_strs <- c("*barcodes.tsv($|.gz)", "*features.tsv($|.gz)", "*matrix.mtx($|.gz)")
    out <- sapply(grep_strs, function(grep_str) {
      file_names <- grep(pattern = grep_str, x = fs, value = TRUE)
      if (length(file_names) >= 2L) {
        stop(paste0("There are multiple ", grep_str, " files within the directory ", curr_directory, "."))
      }
      if (length(file_names) == 0L) {
        stop(paste0("The directory ", curr_directory, " contains zero ", grep_str, " files."))
      }
      return(paste0(curr_directory, "/", file_names))
    }) |> stats::setNames(c("barcodes", "features", "matrix"))
  }, simplify = FALSE) |> stats::setNames(NULL)

  # 3. obtain the barcodes, features, and matrix fp vectors
  barcode_fps <- sapply(X = input_files, FUN = function(i) i[["barcodes"]])
  feature_fps <- sapply(X = input_files, FUN = function(i) i[["features"]])
  matrix_fps <- sapply(X = input_files, FUN = function(i) i[["matrix"]])

  # 4. obtain the feature data frame; determine the ids of each feature
  feature_df <- data.table::fread(file = feature_fps[1],
                                  colClasses = c("character", "character", "character"),
                                  col.names = c("feature_id", "feature_name", "modality"))
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
  batch <- rep(paste0("b", seq_along(n_cells_per_matrix)), n_cells_per_matrix) |> factor()

  # 8. initialize the sceptre object
  cat("Creating the sceptre object.")
  sceptre_object <- import_data(response_matrix = out_mats[["Gene Expression"]],
                                grna_matrix = out_mats[["CRISPR Guide Capture"]],
                                grna_group_data_frame = grna_group_data_frame,
                                moi = moi,
                                extra_covariates = data.frame(batch = batch))
  gc() |> invisible()
  cat(crayon::green(' \u2713\n'))
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
