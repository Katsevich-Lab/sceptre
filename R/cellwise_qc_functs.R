determine_cells_to_retain <- function(sceptre_object, response_n_umis_range, additional_cells_to_remove) {
  cells_to_exclude <- integer()

  # 1. compute the cells to retain based on n_umis_range
  n_umis_vector <- sceptre_object@covariate_data_frame$response_n_umis
  expression_cutoffs <- stats::quantile(x = n_umis_vector, probs = response_n_umis_range)
  cells_to_exlude_1 <- which(n_umis_vector < expression_cutoffs[1] | n_umis_vector > expression_cutoffs[2])
  cells_to_exclude <- union(cells_to_exclude, cells_to_exlude_1)

  # 2. compute cells to retain based on cells containing multiple grnas (if in low MOI)
  if (sceptre_object@low_moi) {
    cells_to_exclude <- union(cells_to_exclude, sceptre_object@multiple_grnas)
  }

  # 3. remove additional cells specified by the user
  cells_to_exclude <- union(cells_to_exclude, additional_cells_to_remove)

  # 4. finally, determine the set of cells to retain, and update the sceptre_object
  n_cells <- ncol(sceptre_object@response_matrix)
  cells_to_retain <- seq(1L, n_cells)[-cells_to_exclude]
  sceptre_object@cells_in_use <- cells_to_retain

  return(sceptre_object)
}
