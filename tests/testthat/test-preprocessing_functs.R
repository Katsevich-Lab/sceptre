test_that("set_matrix_accessibility", {
  set.seed(115)
  n_row <- 23
  n_col <- 38
  m <- matrix(data = rpois(n = n_row * n_col, lambda = 0.5), nrow = n_row, ncol = n_col)

  ## testing row accessible when the input is a matrix
  m_r <- set_matrix_accessibility(matrix_in = m, make_row_accessible = TRUE)
  expect_true(is(m_r, "dgRMatrix"))
  expect_true(all(m_r == m))
  expect_false("i" %in% names(attributes(m_r)))

  ## testing row accessible when the input is a TsparseMatrix
  m_r <- m |> as("TsparseMatrix") |> set_matrix_accessibility(make_row_accessible = TRUE)
  expect_true(is(m_r, "dgRMatrix"))
  expect_true(all(m_r == m))
  expect_false("i" %in% names(attributes(m_r)))

  ## testing row accessible when the input is already a dgRMatrix
  m_r <- m_r |> set_matrix_accessibility(make_row_accessible = TRUE)
  expect_true(is(m_r, "dgRMatrix"))
  expect_true(all(m_r == m))
  expect_false("i" %in% names(attributes(m_r)))

  ## testing col accessible when the input is a matrix
  m_c <- set_matrix_accessibility(matrix_in = m, make_row_accessible = FALSE)
  expect_true(is(m_c, "dgCMatrix"))
  expect_true(all(m_c == m))
  expect_false("j" %in% names(attributes(m_c)))

  ## testing col accessible when the input is a TsparseMatrix
  m_c <- m |> as("TsparseMatrix") |> set_matrix_accessibility(make_row_accessible = FALSE)
  expect_true(is(m_c, "dgCMatrix"))
  expect_true(all(m_c == m))
  expect_false("j" %in% names(attributes(m_c)))

  ## testing col accessible when the input is already a dgCMatrix
  m_c <- m_c |> set_matrix_accessibility(make_row_accessible = FALSE)
  expect_true(is(m_c, "dgCMatrix"))
  expect_true(all(m_c == m))
  expect_false("j" %in% names(attributes(m_c)))
})

