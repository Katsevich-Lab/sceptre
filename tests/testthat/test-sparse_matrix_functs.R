data(response_matrix_lowmoi)
test_that("load_csr_column", {
  response_matrix_lowmoi_r <- set_matrix_accessibility(response_matrix_lowmoi, TRUE)
  samp <- sample(x = seq(1, nrow(response_matrix_lowmoi_r)), size = 10, replace = FALSE)
  for (idx in samp) {
    x1 <- load_csr_row(j = response_matrix_lowmoi_r@j,
                       p = response_matrix_lowmoi_r@p,
                       x = response_matrix_lowmoi_r@x,
                       row_idx = idx,
                       n_cells = ncol(response_matrix_lowmoi_r))
    x2 <- response_matrix_lowmoi[idx,]
    expect_true(all(x1 == x2))
  }
})


combine_perturbations_v2 <- function(perturbation_matrix, grna_group_data_frame) {
  # convert the data frame to a named list
  grna_groups <- as.character(unique(grna_group_data_frame$grna_group))
  grna_groups <- stats::setNames(grna_groups, grna_groups)
  grna_grp_list <- lapply(X = grna_groups, FUN = function(grna_group) {
    dplyr::filter(grna_group_data_frame, grna_group == !!grna_group) |> dplyr::pull(grna_id)
  })
  # Determine which grnas are not contained within grna_grps
  out_grped <- sapply(X = names(grna_grp_list), FUN = function(grp_name) {
    mat_sub <- perturbation_matrix[grna_grp_list[[grp_name]],]
    Matrix::colSums(mat_sub)
  }) |> Matrix::t()
  out_grped <- out_grped >= 1
  return(out_grped)
}


test_that("group_and_threshold_matrix", {
  nrow <- 12
  ncol <- 100
  threshold <- 3
  grna_matrix <- matrix(data = rpois(n = nrow * ncol, lambda = 1),
                        nrow = nrow, ncol = ncol)
  grna_ids <- paste0("grna_", sample(seq(1, nrow)))
  grna_groups <- paste0("group_", seq(1, 4))
  rownames(grna_matrix) <- grna_ids
  grna_group_data_frame <- data.frame(grna_id = grna_ids,
                                      grna_group = rep(grna_groups, each = 3))
  grna_assignments <- assign_grnas_to_cells_highmoi(grna_matrix = grna_matrix,
                                                    threshold = threshold,
                                                    grna_group_data_frame = grna_group_data_frame)
  combined_pert_matrix <- (grna_matrix >= threshold) |>
    combine_perturbations_v2(grna_group_data_frame = grna_group_data_frame)
  for (grna_group in grna_groups) {
    expect_equal(sort(grna_assignments[[grna_group]]), which(combined_pert_matrix[grna_group,]))
  }
})
