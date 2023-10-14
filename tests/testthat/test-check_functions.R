

test_that("check_import_data_inputs", {

  ##### setting up data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  valid_grna_target_data_frame <- make_mock_grna_target_data(
    num_guides_per_target = c(2,3), chr_distances = 1, chr_starts = 1,
    num_nt_guides = 2
  )
  num_cells <- 23
  num_responses <- 18
  num_targets <- nrow(valid_grna_target_data_frame)
  # response and grna matrices are all 1
  valid_response_matrix <- make_mock_response_matrices(num_responses, num_cells, patterns = "one")
  valid_grna_matrix <- make_mock_grna_matrices(
    valid_grna_target_data_frame, num_cells,
    non_nt_patterns = "one", nt_patterns = "one"
  )
  # this is a nice one with basically balanced levels
  valid_extra_covariates <- make_mock_extra_covariates_data_frames(num_cells, patterns = "many_levels")

  ##### running tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##### 0. confirming no problems with all valid data
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ))

  # it should also work if `valid_response_matrix` has no row names
  expect_no_error(
    import_data(
      response_matrix = valid_response_matrix |>
        magrittr::set_rownames(NULL),
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    )
  )

  ##### 1. colnames of `grna_target_data_frame`
  # no column named "grna_id"
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame |>
        dplyr::rename(not_grna_id = grna_id),
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The data frame `grna_target_data_frame` should have columns `grna_id` and `grna_target`"
  )

  # no column named "grna_target"
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame |>
        dplyr::rename(not_grna_target = grna_target),
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The data frame `grna_target_data_frame` should have columns `grna_id` and `grna_target`"
  )

  ##### 2. unique row names for response and grna matrices

  # duplicated rownames in response matrix
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix |>
        magrittr::set_rownames(c("a", "a", paste0("b", 1:(num_responses - 2)))),
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The rownames of the `response_matrix` must be unique"
  )

  # duplicated row names in grna matrix
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix |>
        magrittr::set_rownames(c("a", "a", paste0("b", 1:(num_targets - 2)))),
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The rownames of the `grna_matrix` must be unique"
  )

  # missing row names in grna matrix
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix |>
        magrittr::set_rownames(NULL),
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The column `grna_id` of the `grna_target_data_frame` must be a subset of the row names of the grna expression matrix"
  )

  ##### 5. checking that "&" does not appear in any grna ids and that no grna is named "non-targeting"

  FAIL_id_contains_ampersand_target_df <- valid_grna_target_data_frame |>
    magrittr::inset(1,1,"bad&id")
  FAIL_id_contains_ampersand_grna_matrix <- valid_grna_matrix |>
    magrittr::set_rownames(FAIL_id_contains_ampersand_target_df$grna_id)

  FAIL_id_contains_non_targeting_target_df <- valid_grna_target_data_frame |>
    magrittr::inset(2,1,"non-targeting")
  FAIL_id_contains_non_targeting_grna_matrix <- valid_grna_matrix |>
    magrittr::set_rownames(FAIL_id_contains_non_targeting_target_df$grna_id)

  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = FAIL_id_contains_ampersand_grna_matrix,
      grna_target_data_frame = FAIL_id_contains_ampersand_target_df,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The ampersand character \\(&\\) cannot be present in the gRNA IDs"
  )
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = FAIL_id_contains_non_targeting_grna_matrix,
      grna_target_data_frame = FAIL_id_contains_non_targeting_target_df,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "No individual gRNA can have the ID `non-targeting`"
  )

  ##### 6. ids of grna target data.frame should be a subset of the ids of the grna matrix

  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame |>
        magrittr::inset(3,1,"this grna id is not in the grna matrix"),
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The column `grna_id` of the `grna_target_data_frame` must be a subset of the row names of the grna expression matrix"
  )

  ##### 7. response and grna matrices must have an allowed class

  allowed_classes <- c("matrix", "dgTMatrix", "dgCMatrix", "dgRMatrix")
  # not currently testing for"lgTMatrix", "lgCMatrix", "lgRMatrix" for grna_matrix

  # currently both response_matrix and grna_matrix are dgTMatrix
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix |> as("matrix"),
      grna_matrix = valid_grna_matrix |> as("matrix"),
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    )
  )
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix |> as("CsparseMatrix"),
      grna_matrix = valid_grna_matrix |> as("CsparseMatrix"),
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    )
  )
  expect_no_error(
    check_import_data_inputs(
      response_matrix = set_matrix_accessibility(valid_response_matrix, make_row_accessible = TRUE),
      grna_matrix = set_matrix_accessibility(valid_grna_matrix, make_row_accessible = TRUE),
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    )
  )

  FAIL_dataframe_response_matrix <- valid_response_matrix |>
    as.matrix() |>
    as.data.frame()
  FAIL_dataframe_grna_matrix <- valid_grna_matrix |>
    as.matrix() |>
    as.data.frame()

  expect_error(
    check_import_data_inputs(
      response_matrix = FAIL_dataframe_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "`response_matrix` must be an object of class matrix, dgTMatrix, dgCMatrix, dgRMatrix"
  )
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = FAIL_dataframe_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "`grna_matrix` must be an object of class matrix, dgTMatrix, dgCMatrix, dgRMatrix, lgTMatrix, lgCMatrix, lgRMatrix"
  )

  ##### 8. agreement in number of cells
  FAIL_ncol_resonse_matrix <- valid_response_matrix[,-1]
  FAIL_ncol_grna_matrix <- valid_grna_matrix[,2:num_cells]
  FAIL_nrow_extra_covariates <- valid_extra_covariates[-1,, drop = FALSE]

  expect_error(
    check_import_data_inputs(
      response_matrix = FAIL_ncol_resonse_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The number of cells in the `response_matrix`, `grna_matrix`, and `extra_covariates` \\(if supplied\\) must coincide"
  )
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = FAIL_ncol_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The number of cells in the `response_matrix`, `grna_matrix`, and `extra_covariates` \\(if supplied\\) must coincide"
  )
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = FAIL_nrow_extra_covariates
    ),
    regex = "The number of cells in the `response_matrix`, `grna_matrix`, and `extra_covariates` \\(if supplied\\) must coincide"
  )

  ##### 9. barcodes are valid and agree
  valid_response_matrix_mismatched_colnames <- valid_response_matrix |>
    `colnames<-`(paste0("a", 1:num_cells))
  valid_grna_matrix_mismatched_colnames <- valid_grna_matrix |>
    `colnames<-`(paste0("b", 1:num_cells))
  valid_extra_covariates_mismatched_rownames <- valid_extra_covariates |>
    `rownames<-`(paste0("c", 1:num_cells))

  # if just one of response_matrix, grna_matrix, and extra_covariates have cell names,
  # there will be no error
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix_mismatched_colnames,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = data.frame()
    )
  )
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix_mismatched_colnames,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = data.frame()
    )
  )
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates_mismatched_rownames
    )
  )

  # and there should be no error if they agree in names
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix_mismatched_colnames |>
        `colnames<-`(colnames(valid_grna_matrix_mismatched_colnames)),
      grna_matrix = valid_grna_matrix_mismatched_colnames,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = data.frame()
    )
  )
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix_mismatched_colnames,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates_mismatched_rownames |>
        `rownames<-`(colnames(valid_grna_matrix_mismatched_colnames))
    )
  )
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix_mismatched_colnames,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates_mismatched_rownames |>
        `rownames<-`(colnames(valid_response_matrix_mismatched_colnames))
    )
  )

  # but if at least two have names, then they should agree
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix_mismatched_colnames,
      grna_matrix = valid_grna_matrix_mismatched_colnames,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = data.frame()
    ),
    regex = "You have provided cell barcodes in the `response_matrix` and `grna_matrix`"
  )
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix_mismatched_colnames,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates_mismatched_rownames
    ),
    regex = "You have provided cell barcodes in the `response_matrix` and `extra_covariates`"
  )
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix_mismatched_colnames,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates_mismatched_rownames
    ),
    regex = "You have provided cell barcodes in the `grna_matrix` and `extra_covariates`"
  )

  ##### 10. reserved extra_covariate names

  # just checking a single reserved col name to confirm any check is being run
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = make_mock_extra_covariates_data_frames(num_cells, patterns = "many_columns") |>
        dplyr::mutate(response_n_nonzero = "a")
    ),
    regex = "The covariate names"
  )

  ##### 11. extra_covariates column types
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = make_mock_extra_covariates_data_frames(num_cells, patterns = "many_columns") |>
        dplyr::mutate(char_col = "a", logical_col = FALSE)
    )
  )
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = make_mock_extra_covariates_data_frames(num_cells, patterns = "many_columns") |>
        dplyr::mutate(list_col = c(list(c(1,1,1)), as.list(2:num_cells)))
    ),
    regex = "of the `extra_covariates` data frame should be of type numeric"
  )

  ##### 12. moi values
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    )
  )
  expect_no_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "high",
      extra_covariates = valid_extra_covariates
    )
  )
  expect_error(
    check_import_data_inputs(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "aaa",
      extra_covariates = valid_extra_covariates
    ),
    regex = "`moi` should be either `low` or `high`"
  )
})
