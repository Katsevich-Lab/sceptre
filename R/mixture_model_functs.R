assign_grnas_to_cells_mixture <- function(grna_matrix, grna_group_data_frame) {
  out <- list()

  # 1. make the matrix row-accessible
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
  grna_ids <- rownames(grna_matrix)

  # 2. iterate over grna ids, fitting the mixture model
}
