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


data(response_matrix_highmoi)
data(grna_matrix_highmoi)
data(extra_covariates_highmoi)
data(grna_target_data_frame_highmoi)
data(gene_names_highmoi)

sceptre_object <- import_data(
  response_matrix = response_matrix_highmoi,
  grna_matrix = grna_matrix_highmoi,
  grna_target_data_frame = grna_target_data_frame_highmoi,
  moi = "high",
  extra_covariates = extra_covariates_highmoi,
  response_names = gene_names_highmoi)

write_sceptre_object_to_cellranger_format(sceptre_object = sceptre_object,
                                          directory = "~/research_code/sceptre/inst/extdata/highmoi_example")
