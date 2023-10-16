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

test_that("convert_covariate_df_to_design_matrix", {
  set.seed(187)
  n <- 12
  covariate_data_frame <- data.frame(
    x = rnorm(n), y = 1:n, z = factor(rep(0:1, each = n/2), levels = 0:1)
  )
  fmla <- formula("~ x*z + log(y + 1)")

  ## testing errors with Inf or NA values
  FAIL_bad_values_covariate_data_frame <- covariate_data_frame
  FAIL_bad_values_covariate_data_frame[1,1] <- Inf

  expect_error(
    convert_covariate_df_to_design_matrix(FAIL_bad_values_covariate_data_frame,
                                          formula_object = fmla),
    regex = "contains entries that are -Inf, Inf, or NA"
  )

  FAIL_bad_values_covariate_data_frame[1,1] <- NA
  expect_error(
    convert_covariate_df_to_design_matrix(FAIL_bad_values_covariate_data_frame,
                                          formula_object = fmla),
    regex = "Some rows of `covariate_data_frame` were lost"
  )

  ## now testing if it's the formula that causes these
  FAIL_fmla <- formula("~ log(y - 1)") # should lead to an  -Inf from log(0)
  expect_error(
    convert_covariate_df_to_design_matrix(covariate_data_frame,
                                          formula_object = FAIL_fmla),
    regex = "has been applied contains entries that are -Inf, Inf, or NA"
  )

  ## testing low rank
  FAIL_low_Rank_covariate_data_frame <- covariate_data_frame |>
    dplyr::mutate(low_rank_col = 2*x - log(y + 1))
  # the formula doesn't use the low rank column so there should be no error currently
  expect_no_error(
    convert_covariate_df_to_design_matrix(FAIL_low_Rank_covariate_data_frame,
                                          formula_object = fmla)
  )
  expect_error(
    convert_covariate_df_to_design_matrix(
      FAIL_low_Rank_covariate_data_frame,
      formula_object = formula("~ x*z + log(y + 1) + low_rank_col")
    ),
    regex = "contains redundant information"
  )

  ## testing correctness
  results <- convert_covariate_df_to_design_matrix(covariate_data_frame,
                                                   formula_object = fmla)
  expect_equal(
    ncol(results), 5 # with intercept
  )
  expect_true(
    all(results[,1] == 1)
  )
  expect_equal(
    results[,"x"] |> as.numeric(), covariate_data_frame$x
  )
  expect_equal(
    results[,"z1"] |> as.numeric(), covariate_data_frame$z |> as.character() |> as.numeric()
  )
  expect_equal(
    results[,"log(y + 1)"] |> as.numeric(), covariate_data_frame$y |> log1p()
  )
  expect_equal(
    results[,"x:z1"] |> as.numeric(),
    (covariate_data_frame$z |> as.character() |> as.numeric()) * covariate_data_frame$x
  )
})

