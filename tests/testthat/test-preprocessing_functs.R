test_that("set_matrix_accessibility", {
  library(Matrix)
  n_row <- sample(seq(40, 80), 1L) 
  n_col <- sample(seq(40, 80), 1L)
  m <- matrix(data = rpois(n = n_row * n_col, lambda = 0.5), nrow = n_row, ncol = n_col)
  m_r <- set_matrix_accessibility(matrix_in = m, make_row_accessible = TRUE)
  expect_true(is(m_r, "dgRMatrix"))
  expect_true(all(m_r == m))
  
  m <- matrix(data = rpois(n = n_row * n_col, lambda = 0.5), nrow = n_row, ncol = n_col)
  m_c <- set_matrix_accessibility(matrix_in = m, make_row_accessible = FALSE)
  expect_true(is(m_c, "dgCMatrix"))
  expect_true(all(m_c == m))
})


test_that("assign_grnas_to_cells_lowmoi", {
  samp <- sample(x = seq(1, ncol(grna_matrix_lowmoi)))
  grna_matrix_lowmoi <- grna_matrix_lowmoi[,samp]
  grna_assignments <- assign_grnas_to_cells_lowmoi_v2(grna_matrix = grna_matrix_lowmoi,
                                                      calibration_check = TRUE,
                                                      grna_group_data_frame = grna_group_data_frame_lowmoi,
                                                      n_calibration_pairs = NULL)
  # test grna groups
  grna_group_idxs <- grna_assignments$grna_group_idxs
  grna_groups <- names(grna_group_idxs)
  for (grna_group in grna_groups) {
    if (length(grna_group_idxs[[grna_group]]) == 0L) next
    idx <- sample(grna_group_idxs[[grna_group]], 1)
    grna_id <- names(which.max(grna_matrix_lowmoi[,idx]))
    test_grna_group <- grna_group_data_frame_lowmoi |>
      dplyr::filter(grna_id == !!grna_id) |>
      dplyr::pull(grna_group) |> as.character()
    expect_equal(grna_group, test_grna_group)
  }
  
  # test indiv nt
  indiv_nt_idxs <- grna_assignments$indiv_nt_grna_idxs
  indiv_nts <- names(indiv_nt_idxs)
  for (indiv_nt in indiv_nts) {
    idx <- sample(indiv_nt_idxs[[indiv_nt]], 1L)
    grna_id <- names(which.max(grna_matrix_lowmoi[,grna_assignments$all_nt_idxs[idx]]))
    expect_equal(indiv_nt, grna_id)
  }
})


test_that("response precomputation", {
  formula_object <- formula(~log(response_n_umis) + log(response_n_nonzero) + bio_rep + p_mito)
  data(covariate_data_frame_lowmoi)
  X <- model.matrix(data = covariate_data_frame_lowmoi, object = formula_object)
  betas <- c(-5, 0.5, 0.2, 1.1, 0.8, 0.2)
  mus <- exp(as.numeric(X %*% betas))
  theta <- 9
  y <- sapply(X = mus, FUN = function(mu) MASS::rnegbin(n = 1, mu = mu, theta = theta))
  precomp <- perform_response_precomputation(expressions = y,
                                             covariate_matrix = X)
  expect_true(all(abs(betas - precomp$fitted_coefs) < 0.5))
  abs(precomp$theta - theta) < 0.5
})
