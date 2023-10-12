#############################
## unit test `import_data` ##
#############################

# The tests in here are broken up by which subfunction of `import_data`
# they are really testing. In order, they are:
# 1. import_data-check_import_data_inputs
# 2. import_data-auto_compute_cell_covariates
# 3. import_data-set_matrix_accessibility
# 4. import_data-slots (this is just testing the final slots/attributes rather than a subfunction)

# this is testing `check_import_data_inputs` as called by `import_data`
# the section numberings in this test correspond to the numberings
# as they appear in that function
test_that("import_data-check_import_data_inputs", {

  ##### setting up data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  valid_grna_target_data_frame <- make_mock_grna_target_data(
    num_guides_per_target = c(2,3), chr_distances = 1, chr_starts = 1,
    num_nt_guides = 2
  )
  num_cells <- 23
  num_responses <- 18
  num_targets <- nrow(valid_grna_target_data_frame)
  # response and grna matrices are all 1
  valid_response_matrix <- make_mock_response_matrix_list(num_responses, num_cells)$all_one
  valid_grna_matrix <- make_mock_grna_matrix_list(
    valid_grna_target_data_frame, num_cells
  )$non_nt_all_one_and_nt_all_one
  # this is a nice one with basically balanced levels
  valid_extra_covariates <- make_mock_extra_covariates_list(num_cells)$many_columns |>
    dplyr::select(batch)

  ##### running tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ##### 0. confirming no problems with all valid data
  expect_no_error(
    import_data(
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
    import_data(
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
    import_data(
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
    import_data(
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
    import_data(
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
    import_data(
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
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = FAIL_id_contains_ampersand_grna_matrix,
      grna_target_data_frame = FAIL_id_contains_ampersand_target_df,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The ampersand character \\(&\\) cannot be present in the gRNA IDs"
  )
  expect_error(
    import_data(
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
    import_data(
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
    import_data(
      response_matrix = valid_response_matrix |> as("matrix"),
      grna_matrix = valid_grna_matrix |> as("matrix"),
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    )
  )
  expect_no_error(
    import_data(
      response_matrix = valid_response_matrix |> as("CsparseMatrix"),
      grna_matrix = valid_grna_matrix |> as("CsparseMatrix"),
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    )
  )
  expect_no_error(
    import_data(
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
    import_data(
      response_matrix = FAIL_dataframe_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "`response_matrix` must be an object of class matrix, dgTMatrix, dgCMatrix, dgRMatrix"
  )
  expect_error(
    import_data(
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
    import_data(
      response_matrix = FAIL_ncol_resonse_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The number of cells in the `response_matrix`, `grna_matrix`, and `extra_covariates` \\(if supplied\\) must coincide"
  )
  expect_error(
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = FAIL_ncol_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The number of cells in the `response_matrix`, `grna_matrix`, and `extra_covariates` \\(if supplied\\) must coincide"
  )
  expect_error(
    import_data(
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
    import_data(
      response_matrix = valid_response_matrix_mismatched_colnames,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low"
    )
  )
  expect_no_error(
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix_mismatched_colnames,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low"
    )
  )
  expect_no_error(
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates_mismatched_rownames
    )
  )

  # and there should be no error if they agree in names
  expect_no_error(
    import_data(
      response_matrix = valid_response_matrix_mismatched_colnames |>
        `colnames<-`(colnames(valid_grna_matrix_mismatched_colnames)),
      grna_matrix = valid_grna_matrix_mismatched_colnames,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low"
    )
  )
  expect_no_error(
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix_mismatched_colnames,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates_mismatched_rownames |>
        `rownames<-`(colnames(valid_grna_matrix_mismatched_colnames))
    )
  )
  expect_no_error(
    import_data(
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
    import_data(
      response_matrix = valid_response_matrix_mismatched_colnames,
      grna_matrix = valid_grna_matrix_mismatched_colnames,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low"
    ),
    regex = "You have provided cell barcodes in the `response_matrix` and `grna_matrix`"
  )
  expect_error(
    import_data(
      response_matrix = valid_response_matrix_mismatched_colnames,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates_mismatched_rownames
    ),
    regex = "You have provided cell barcodes in the `response_matrix` and `extra_covariates`"
  )
  expect_error(
    import_data(
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
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = make_mock_extra_covariates_list(num_cells)$many_columns |>
        dplyr::mutate(response_n_nonzero = "a")
    ),
    regex = "The covariate names"
  )

  ##### 11. extra_covariates column types
  expect_no_error(
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = make_mock_extra_covariates_list(num_cells)$many_columns |>
        dplyr::mutate(char_col = "a", logical_col = FALSE)
    )
  )
  expect_error(
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = make_mock_extra_covariates_list(num_cells)$many_columns |>
        dplyr::mutate(list_col = c(list(c(1,1,1)), as.list(2:num_cells)))
    ),
    regex = "of the `extra_covariates` data frame should be of type numeric"
  )

  ##### 12. moi values
  expect_no_error(
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    )
  )
  expect_no_error(
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "high",
      extra_covariates = valid_extra_covariates
    )
  )
  expect_error(
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "aaa",
      extra_covariates = valid_extra_covariates
    ),
    regex = "`moi` should be either `low` or `high`"
  )
})

# in `import_data`, the automatically computed covariates come from
# `auto_compute_cell_covariates`. This function first computes covariates from
# `response_matrix` and then computes covariates from `grna_matrix` before
# cbind'ing these together along with `extra_covariates`.
# This test first fixes `grna_matrix` and checks the calculations for a variety of
# `response_matrix` values, then fixes the latter and varies the former.
test_that("import_data-auto_compute_cell_covariates", {
  set.seed(11123)
  num_cells <- 43
  num_responses <- 21

  # having NT or not doesn't matter for this so it will be included
  fixed_grna_target_data_frame <- make_mock_grna_target_data(c(1,2,4), 1, 1, 5)

  fixed_response_matrix <- matrix(1, num_responses, num_cells)
  fixed_grna_matrix <- matrix(1, nrow(fixed_grna_target_data_frame), num_cells) |>
    `rownames<-`(fixed_grna_target_data_frame$grna_id)
  fixed_extra_covariates <- data.frame(x = rpois(num_cells, 5))

  ##### 1. testing the appending of `extra_covariates`
  # no extra_covariates so we should just have the default 4 columns in
  # @covariate_data_frame
  sceptre_object <- import_data(
    response_matrix = fixed_response_matrix,
    grna_matrix = fixed_grna_matrix,
    grna_target_data_frame = fixed_grna_target_data_frame,
    moi = "low"  # moi doesn't matter for this test
  )
  expect_equal(
    names(sceptre_object@covariate_data_frame),
    c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis")
  )

  # with extra_covariates those column names and values should appear now too
  sceptre_object <- import_data(
    response_matrix = fixed_response_matrix,
    grna_matrix = fixed_grna_matrix,
    grna_target_data_frame = fixed_grna_target_data_frame,
    moi = "low",  # moi doesn't matter for this test
    extra_covariates = fixed_extra_covariates
  )
  expect_equal(
    names(sceptre_object@covariate_data_frame),
    c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis", "x")
  )
  expect_equal(
    sceptre_object@covariate_data_frame$x,
    fixed_extra_covariates$x
  )

  ##### 2. testing computations for `response_matrix`
  for(response_matrix in make_mock_response_matrix_list(num_responses, num_cells)) {
    sceptre_object <- import_data(
      response_matrix = response_matrix,
      grna_matrix = fixed_grna_matrix,
      grna_target_data_frame = fixed_grna_target_data_frame,
      moi = "low"  # moi doesn't matter for this test
    )
    expect_equal(
      sceptre_object@covariate_data_frame$response_n_nonzero,
      response_matrix |> apply(2, function(cell) sum(cell > 0))
    )
    expect_equal(
      sceptre_object@covariate_data_frame$response_n_umis,
      response_matrix |> colSums()
    )
  }

  ##### 3. testing computations for `grna_matrix`
  for(grna_matrix in make_mock_grna_matrix_list(fixed_grna_target_data_frame, num_cells)) {
    sceptre_object <- import_data(
      response_matrix = fixed_response_matrix,
      grna_matrix = grna_matrix,
      grna_target_data_frame = fixed_grna_target_data_frame,
      moi = "low"  # moi doesn't matter for this test
    )
    expect_equal(
      sceptre_object@covariate_data_frame$grna_n_nonzero,
      grna_matrix |> apply(2, function(cell) sum(cell > 0))
    )
    expect_equal(
      sceptre_object@covariate_data_frame$grna_n_umis,
      grna_matrix |> colSums()
    )
  }
})

# just confirming that sceptre_object@response_matrix has the right values, class, and attributes
# after `set_matrix_accessibility`
test_that("import_data-set_matrix_accessibility", {
  set.seed(132)
  num_cells <- 27
  num_responses <- 12
  grna_target_data_frame <- make_mock_grna_target_data(c(1,2,3), 1, 1, 5)
  grna_matrix <- matrix(rpois(nrow(grna_target_data_frame) * num_cells, 1), ncol = num_cells) |>
    `rownames<-`(grna_target_data_frame$grna_id)
  response_matrix <- matrix(rpois(num_responses * num_cells, 1), ncol = num_cells) |>
    as("TsparseMatrix")

  sceptre_object <- import_data(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low"  # moi doesn't matter for this test
  )

  expect_equal(sceptre_object@response_matrix@x, response_matrix@x) # values ok?
  expect_equal(class(sceptre_object@response_matrix) |> as.character(), "dgRMatrix") # class ok?
  expect_true(all(c("p", "j") %in% names(attributes(sceptre_object@response_matrix)))) # attr ok?
  expect_false("i" %in% names(attributes(sceptre_object@response_matrix)))
})


# testing that the various attributes of the returned sceptre object all came out ok
test_that("import_data-slots", {
  set.seed(12321)
  num_cells <- 29
  num_responses <- 17
  grna_target_data_frame <- make_mock_grna_target_data(c(1,3,1), 1, 1, 6)
  grna_matrix <- matrix(rpois(nrow(grna_target_data_frame) * num_cells, 1), ncol = num_cells) |>
    `rownames<-`(grna_target_data_frame$grna_id)
  response_matrix <- matrix(rpois(num_responses * num_cells, 1), ncol = num_cells) |>
    as("TsparseMatrix")
  extra_covariates <- data.frame(x = rep("aaa", num_cells))

  sceptre_object_low_no_ec_with_response_names <- import_data(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low",
    response_names = paste0("rrr", 1:num_responses)
  )
  sceptre_object_high_with_ec <- import_data(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "high",
    extra_covariates = extra_covariates
  )

  # most of the slots are so simple that only cursory tests are done, mostly to confirm
  # that future changes do not change the API

  all_slots <- c("response_matrix", "grna_matrix", "covariate_data_frame", "covariate_matrix",
                 "grna_target_data_frame", "low_moi", "user_specified_covariates",
                 "response_names", "discovery_pairs", "positive_control_pairs", "formula_object",
                 "side_code", "fit_parametric_curve", "control_group_complement",
                 "run_permutations", "n_nonzero_trt_thresh", "n_nonzero_cntrl_thresh",
                 "B1", "B2", "B3", "grna_grouping_strategy", "grna_assignment_method",
                 "grna_assignment_hyperparameters", "multiple_testing_alpha", "multiple_testing_method",
                 "cell_removal_metrics", "mitochondrial_gene", "M_matrix", "n_nonzero_tot_vector",
                 "discovery_pairs_with_info", "positive_control_pairs_with_info",
                 "negative_control_pairs", "initial_grna_assignment_list", "grna_assignments_raw",
                 "grna_assignments", "grnas_per_cell", "cells_w_multiple_grnas",
                 "cells_in_use", "n_ok_discovery_pairs", "n_ok_positive_control_pairs",
                 "calibration_group_size", "n_calibration_pairs", "response_precomputations",
                 "last_function_called", "functs_called", "calibration_result", "power_result",
                 "discovery_result")
  expect_equal(slotNames(sceptre_object_low_no_ec_with_response_names), all_slots)
  expect_equal(slotNames(sceptre_object_high_with_ec), all_slots)

  ##### slots set in section 4 of `import_data`

  expect_equal(sceptre_object_high_with_ec@response_matrix, set_matrix_accessibility(response_matrix))
  expect_equal(sceptre_object_high_with_ec@grna_matrix, grna_matrix)
  expect_equal(sceptre_object_low_no_ec_with_response_names@response_names, paste0("rrr", 1:num_responses))

  expect_true(sceptre_object_low_no_ec_with_response_names@low_moi)
  expect_false(sceptre_object_high_with_ec@low_moi)

  expect_equal(sceptre_object_low_no_ec_with_response_names@user_specified_covariates, character(0))
  expect_equal(sceptre_object_high_with_ec@user_specified_covariates, "x")

  ##### slots set in section 5 of `import_data`

  expect_equal(sceptre_object_high_with_ec@last_function_called, "import_data")
  expect_equal(names(sceptre_object_high_with_ec@functs_called), c("import_data", "set_analysis_parameters",
                                                          "assign_grnas", "run_qc", "run_calibration_check",
                                                          "run_power_check", "run_discovery_analysis"))
  expect_equal(sceptre_object_high_with_ec@functs_called |> as.logical(), rep(c(TRUE, FALSE), c(1, 6)))
})
