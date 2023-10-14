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


test_that("response precomputation", {
  formula_object <- formula(~log(response_n_umis) + log(response_n_nonzero) + bio_rep + p_mito)
  data(covariate_data_frame_lowmoi)
  X <- model.matrix(data = covariate_data_frame_lowmoi, object = formula_object)
  betas <- c(-5, 0.5, 0.2, 1.1, 0.8, 0.2)
  mus <- exp(as.numeric(X %*% betas))
  theta <- 9
  y <- sapply(X = mus, FUN = function(mu) MASS::rnegbin(n = 1, mu = mu, theta = theta))
  precomp <- perform_response_precomputation(expressions = y,
                                             covariate_matrix = X,
                                             regression_method = "nb_glm")
  expect_true(all(abs(betas - precomp$fitted_coefs) < 0.5))
  abs(precomp$theta - theta) < 0.5
})
