
test_that("auto_compute_cell_covariates", {
  set.seed(11123)
  num_cells <- 43
  num_responses <- 21

  # only needed to make grna_matrix
  grna_target_data_frame <- make_mock_grna_target_data(c(1,2,4), 1, 1, 5)
  fixed_response_matrix <- make_mock_response_matrices(num_responses, num_cells, patterns = "column")
  fixed_grna_matrix <- make_mock_grna_matrices(grna_target_data_frame, num_cells, non_nt_patterns = "column", nt_pattern = "row") |>
    `rownames<-`(grna_target_data_frame$grna_id)
  extra_covariates_big <- make_mock_extra_covariates_data_frames(num_cells, patterns = "many_columns")
  extra_covariates_empty <- data.frame()

  ## testing response_matrix
  # only `response_n_nonzero` and `response_n_umis` are tested for this part
  for(response_matrix in make_mock_response_matrices(num_responses, num_cells, patterns = "all")) {
    results <- auto_compute_cell_covariates(
      response_matrix = response_matrix,
      grna_matrix = fixed_grna_matrix,
      extra_covariates = extra_covariates_big,
      response_names = NA_character_  # this is what it is by default in `import_data`
    )

    expect_equal(
      results$response_n_nonzero,
      response_matrix |> apply(2, function(cell) sum(cell > 0))
    )
    expect_equal(
      results$response_n_umis,
      response_matrix |> colSums()
    )
  }

  ## testing grna_matrix
  # only `response_n_nonzero` and `response_n_umis` are tested for this part
  for(grna_matrix in make_mock_grna_matrices(grna_target_data_frame, num_cells, non_nt_patterns = "all", nt_patterns = "column")) {
    results <- auto_compute_cell_covariates(
      response_matrix = fixed_response_matrix,
      grna_matrix = grna_matrix,
      extra_covariates = extra_covariates_big,
      response_names = NA_character_
    )

    expect_equal(
      results$grna_n_nonzero,
      grna_matrix |> apply(2, function(cell) sum(cell > 0))
    )
    expect_equal(
      results$grna_n_umis,
      grna_matrix |> colSums()
    )
  }

  ## testing extra_covariates
  results_big <- auto_compute_cell_covariates(
    response_matrix = fixed_response_matrix,
    grna_matrix = fixed_grna_matrix,
    extra_covariates = extra_covariates_big,
    response_names = NA_character_
  )
  expect_equal(
    results_big[,names(extra_covariates_big)],
    extra_covariates_big
  )

  results_empty <- auto_compute_cell_covariates(
    response_matrix = fixed_response_matrix,
    grna_matrix = fixed_grna_matrix,
    extra_covariates = extra_covariates_empty,
    response_names = NA_character_
  )
  expect_equal(
    colnames(results_empty),
    c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis")
  )

  ## testing response_names
  results_no_names <- auto_compute_cell_covariates(
    response_matrix = fixed_response_matrix,
    grna_matrix = fixed_grna_matrix,
    extra_covariates = extra_covariates_big,
    response_names = paste0("aaa", 1:num_responses)
  )

  expect_equal(
    results_no_names$response_n_nonzero,
    fixed_response_matrix |> apply(2, function(cell) sum(cell > 0))
  )
  expect_equal(
    results_no_names$response_n_umis,
    fixed_response_matrix |> colSums()
  )
  expect_equal(
    results_no_names$grna_n_nonzero,
    fixed_grna_matrix |> apply(2, function(cell) sum(cell > 0))
  )
  expect_equal(
    results_no_names$grna_n_umis,
    fixed_grna_matrix |> colSums()
  )
})

# this function is never called directly; `sceptre_object@covariate_data_frame`
# will always have at least 4 columns from the response and grna matrices
# from `auto_compute_cell_covariates` in `import_data`
test_that("auto_construct_formula_object", {
  covariate_data_frame <- data.frame(
    aaa_n_umis = 0:5, n_nonzero_aaa = 1:6, x = rnorm(6),
    response_p_mito = 1.23, grna_n_umis = 999
  )

  fmla_with_grna <- auto_construct_formula_object(covariate_data_frame, include_grna_covariates = TRUE)
  fmla_without_grna <- auto_construct_formula_object(covariate_data_frame, include_grna_covariates = FALSE)

  expect_equal(
    as.character(fmla_with_grna)[2],
    "log(aaa_n_umis + 1) + log(n_nonzero_aaa) + x + log(grna_n_umis)"
  )
  expect_equal(
    as.character(fmla_without_grna)[2],
    "log(aaa_n_umis + 1) + log(n_nonzero_aaa) + x"
  )
})

