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
