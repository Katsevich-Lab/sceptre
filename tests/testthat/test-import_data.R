#############################
## unit test `import_data` ##
#############################

# TODO modify helper functions so i can get specific individual matrix returns rather
# than always generating the whole list and then subsetting

rm(list=ls())
source("~/academics/katsevich-lab/sceptre/tests/testthat/helper-mock-data.R")

## pretty comprehensive set of grna_target_data_frames aside from varying `start` and `end`

## these do not explore different values for `start` and `end`
# and `chr` is also not the focus.
grna_target_data_frame_list <- list(
  ## single target datasets, with and without NT ~~~~~~~~~~
  # just one grna total, no NT
  make_mock_grna_target_data(
    num_guides_per_target = 1, chr_distances = 2,
    chr_starts = 1, num_nt_guides = 0
  ) |> dplyr::filter(grna_target == grna_target[1]),
  # 10 grnas all still for just one target, no NT
  make_mock_grna_target_data(
    num_guides_per_target = 10, chr_distances = 2,
    chr_starts = 1, num_nt_guides = 0
  ) |> dplyr::filter(grna_target == grna_target[1]),
  # one target, one NT
  make_mock_grna_target_data(
    num_guides_per_target = 1, chr_distances = 2,
    chr_starts = 1, num_nt_guides = 1
  ) |> dplyr::filter(grna_target %in% c(grna_target[1], "non-targeting")),
  # one target, several NT
  make_mock_grna_target_data(
    num_guides_per_target = 1, chr_distances = 2,
    chr_starts = 1, num_nt_guides = 5
  ) |> dplyr::filter(grna_target %in% c(grna_target[1], "non-targeting")),

  ## two targets, with and without NT ~~~~~
  # one guide per target, no NT
  make_mock_grna_target_data(
    num_guides_per_target = 1, chr_distances = 2,
    chr_starts = 1, num_nt_guides = 0
  ),
  # one guide per target, one NT
  make_mock_grna_target_data(
    num_guides_per_target = 1, chr_distances = 2,
    chr_starts = 1, num_nt_guides = 1
  ),
  # one guide per target, several NT
  make_mock_grna_target_data(
    num_guides_per_target = 1, chr_distances = 2,
    chr_starts = 1, num_nt_guides = 6
  ),
  # 5 guides per target, no NT
  make_mock_grna_target_data(
    num_guides_per_target = 5, chr_distances = 2,
    chr_starts = 1, num_nt_guides = 0
  ),
  # 5 guides per target, one NT
  make_mock_grna_target_data(
    num_guides_per_target = 5, chr_distances = 2,
    chr_starts = 1, num_nt_guides = 1
  ),
  # 5 guides per target, several NT
  make_mock_grna_target_data(
    num_guides_per_target = 5, chr_distances = 2,
    chr_starts = 1, num_nt_guides = 7
  ),
  # one target has 1 guide, the other has several, no NT
  make_mock_grna_target_data(
    num_guides_per_target = c(1, 10), chr_distances = 2,
    chr_starts = 1, num_nt_guides = 0
  ) |> dplyr::filter(chr != tail(chr)[1]),
  # one target has 1 guide, the other has several, 1 NT
  make_mock_grna_target_data(
    num_guides_per_target = c(1, 10), chr_distances = 2,
    chr_starts = 1, num_nt_guides = 1
  )[-c(12:22),],
  # one target has 1 guide, the other has several, several NT
  make_mock_grna_target_data(
    num_guides_per_target = c(1, 10), chr_distances = 2,
    chr_starts = 1, num_nt_guides = 8
  )[-c(12:22),],

  ## many targets, with and without NT ~~~~~
  # many targets, one has just one guide, no NT
  make_mock_grna_target_data(
    num_guides_per_target = c(1, 2, 3), chr_distances = 2,
    chr_starts = 1, num_nt_guides = 0
  ),
  # many targets, one has just one guide, 1 NT
  make_mock_grna_target_data(
    num_guides_per_target = c(1, 2, 3), chr_distances = 2,
    chr_starts = 1, num_nt_guides = 1
  ),
  # many targets, one has just one guide, many NT
  make_mock_grna_target_data(
    num_guides_per_target = c(1, 2, 3), chr_distances = 2,
    chr_starts = 1, num_nt_guides = 10
  ),
  # many targets, all have at least two guides, no NT
  make_mock_grna_target_data(
    num_guides_per_target = c(3, 4, 3), chr_distances = 2,
    chr_starts = 1, num_nt_guides = 0
  ),
  # many targets, all have at least two guides, 1 NT
  make_mock_grna_target_data(
    num_guides_per_target = c(3, 4, 3), chr_distances = 2,
    chr_starts = 1, num_nt_guides = 1
  ),
  # many targets, all have at least two guides, many NT
  make_mock_grna_target_data(
    num_guides_per_target = c(3, 4, 3), chr_distances = 2,
    chr_starts = 1, num_nt_guides = 11
  )
)



# really this is just testing `check_import_data_inputs`
test_that("import_data-invalid-inputs", {

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
  expect_no_error(import_data(
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
  # TODO what about "lgTMatrix", "lgCMatrix", "lgRMatrix" for grna_matrix?

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










  expect_error(
    import_data(
      response_matrix = valid_response_matrix,
      grna_matrix = valid_grna_matrix,
      grna_target_data_frame = valid_grna_target_data_frame,
      moi = "low",
      extra_covariates = valid_extra_covariates
    ),
    regex = "The ampersand character (&) cannot be present in the gRNA IDs"
  )

})
# function(response_matrix, grna_matrix, grna_target_data_frame, moi, extra_covariates) {
#
#   # 8. check for agreement in number of cells
#   check_ncells <- ncol(response_matrix) == ncol(grna_matrix)
#   if (!is.null(extra_covariates)) {
#     check_ncells <- check_ncells && (ncol(response_matrix) == nrow(extra_covariates))
#   }
#   if (!check_ncells) {
#     stop("The number of cells in the `response_matrix`, `grna_matrix`, and `extra_covariates` (if supplied) must coincide.")
#   }
#
#   # 9. if applicable, check that the cell barcodes match across the grna, gene, and covariate matrices
#   check_barcodes_provided <- function(barcodes) {
#     !is.null(barcodes) && !all(grepl(pattern = "^[0-9]+$", x = barcodes))
#   }
#   response_cell_barcodes <- colnames(response_matrix)
#   grna_cell_barcodes <- colnames(grna_matrix)
#   covariate_cell_barcodes <- rownames(extra_covariates)
#   if (check_barcodes_provided(response_cell_barcodes) &&
#       check_barcodes_provided(grna_cell_barcodes) &&
#       check_barcodes_provided(covariate_cell_barcodes)) {
#     if (!(identical(response_cell_barcodes, grna_cell_barcodes) &&
#           identical(response_cell_barcodes, covariate_cell_barcodes))) {
#       stop("You have provided cell barcodes in the `response_matrix`, `grna_matrix`, and `extra_covariates`. These cell barcodes must have the same ordering across objects.")
#     }
#   }
#
#   # 10. check that column names of extra_covariates are not already taken
#   reserved_covariate_names <- c("response_n_nonzero", "response_n_umis", "response_p_mito", "grna_n_nonzero", "grna_n_umis")
#   extra_covariate_names <- colnames(extra_covariates)
#   if (any(extra_covariate_names %in% reserved_covariate_names)) {
#     stop("The covariate names `response_n_nonzero`, `response_n_umis`, `response_p_mito`, `grna_n_nonzero`, and `grna_n_umis` are reserved. Change the column names of the `extra_covariates` data frame.")
#   }
#
#   # 11. verify that the types of the extra covariates are acceptable
#   for (extra_covariate_name in extra_covariate_names) {
#     v <- extra_covariates[,extra_covariate_name]
#     accept_type <- methods::is(v, "numeric") || methods::is(v, "character") || methods::is(v, "factor")
#     if (!accept_type) {
#       stop(paste0("The column `", extra_covariate_name, "` of the `extra_covariates` data frame should be of type numeric, character, or factor."))
#     }
#   }
#
#   # 12. verify that moi is specified
#   if (!(moi %in% c("low", "high"))) {
#     stop("`moi` should be either `low` or `high`.")
#   }
#
#   return(NULL)
# }


# for testing computation, only grna_matrix and response_matrix change here
# for testing inputs, we need more variation
test_that("import_data-valid-inputs", {








  num_cells <- 50
  num_responses <- 12

  set.seed(112)
  response_matrix <- matrix(rpois(num_responses * num_cells, 1), num_responses, num_cells) |>
    as("TsparseMatrix")

  # makes a matrix of iid Popis(1) in the right size
  random_grna_matrix_from_target_data_frame <- function(grna_target_data_frame, num_cells) {
    matrix(rpois(nrow(grna_target_data_frame) * num_cells, 1), nrow(grna_target_data_frame), num_cells) |>
      `dimnames<-`(list(grna_target_data_frame$grna_id)) |>
      as("TsparseMatrix")
  }

  extra_covariates <- make_mock_extra_covariates_list(num_cells)$many_columns

  for(i in seq_along(grna_target_data_frame_list)) {
    set.seed(i)
    grna_target_data_frame <- grna_target_data_frame_list[[i]]
    grna_matrix <- random_grna_matrix_from_target_data_frame(grna_target_data_frame, num_cells)

    sceptre_object_low <- import_data(
      response_matrix = response_matrix,
      grna_matrix = grna_matrix,
      grna_target_data_frame = grna_target_data_frame,
      moi = "low",
      extra_covariates = extra_covariates
    )
    sceptre_object_high <- import_data(
      response_matrix = response_matrix,
      grna_matrix = grna_matrix,
      grna_target_data_frame = grna_target_data_frame,
      moi = "high",
      extra_covariates = extra_covariates
    )

    # 1. inputs
    # - for now, all these are valid inputs so all should pass

    # 2. covariates
    # the focus of this test block is grna_target_data_frame so not super relevant but I'll get them started here

    # test for low MOI
    expect_equal(
      sceptre_object_low@covariate_data_frame$response_n_nonzero,
      response_matrix |> apply(2, function(cell) sum(cell > 0))
    )
    expect_equal(
      sceptre_object_low@covariate_data_frame$response_n_umis,
      response_matrix |> colSums()
    )
    expect_equal(
      sceptre_object_low@covariate_data_frame$grna_n_nonzero,
      grna_matrix |> apply(2, function(cell) sum(cell > 0))
    )
    expect_equal(
      sceptre_object_low@covariate_data_frame$grna_n_umis,
      grna_matrix |> colSums()
    )
    expect_equal(
      sceptre_object_low@covariate_data_frame[,names(extra_covariates)],
      extra_covariates
    )

    # test for high MOI
    expect_equal(
      sceptre_object_high@covariate_data_frame$response_n_nonzero,
      response_matrix |> apply(2, function(cell) sum(cell > 0))
    )
    expect_equal(
      sceptre_object_high@covariate_data_frame$response_n_umis,
      response_matrix |> colSums()
    )
    expect_equal(
      sceptre_object_high@covariate_data_frame$grna_n_nonzero,
      grna_matrix |> apply(2, function(cell) sum(cell > 0))
    )
    expect_equal(
      sceptre_object_high@covariate_data_frame$grna_n_umis,
      grna_matrix |> colSums()
    )
    expect_equal(
      sceptre_object_high@covariate_data_frame[,names(extra_covariates)],
      extra_covariates
    )


    # 3. response matrix
    expect_equal(sceptre_object_low@response_matrix@x, response_matrix@x)
    expect_equal(class(sceptre_object_low@response_matrix) |> as.character(), "dgRMatrix")

    expect_equal(sceptre_object_high@response_matrix@x, response_matrix@x)
    expect_equal(class(sceptre_object_high@response_matrix) |> as.character(), "dgRMatrix")

    # 4. attributes (that haven't already been tested)
    expect_equal(sceptre_object_low@grna_matrix, grna_matrix)
    expect_equal(sceptre_object_high@grna_matrix, grna_matrix)

    grna_target_data_frame_with_chars <- dplyr::mutate(
      grna_target_data_frame,
      grna_id = as.character(grna_id), grna_target = as.character(grna_target)
    )
    expect_equal(sceptre_object_low@grna_target_data_frame, grna_target_data_frame_with_chars)
    expect_equal(sceptre_object_high@grna_target_data_frame, grna_target_data_frame_with_chars)

    expect_true(is.na(sceptre_object_low@response_names))
    expect_true(is.na(sceptre_object_high@response_names))

    expect_true(sceptre_object_low@low_moi)
    expect_true(!sceptre_object_high@low_moi)

    expect_equal(sceptre_object_low@user_specified_covariates, colnames(extra_covariates))
    expect_equal(sceptre_object_high@user_specified_covariates, colnames(extra_covariates))

    # 5. flags
    expect_equal(sceptre_object_low@last_function_called, "import_data")
    expect_equal(sceptre_object_high@last_function_called, "import_data")

    expect_equal(names(sceptre_object_low@functs_called), c("import_data", "set_analysis_parameters",
                                                            "assign_grnas", "run_qc", "run_calibration_check",
                                                            "run_power_check", "run_discovery_analysis"))
    expect_equal(names(sceptre_object_high@functs_called), c("import_data", "set_analysis_parameters",
                                                            "assign_grnas", "run_qc", "run_calibration_check",
                                                            "run_power_check", "run_discovery_analysis"))

    expect_equal(sceptre_object_low@functs_called |> as.logical(), rep(c(TRUE, FALSE), c(1, 6)))
    expect_equal(sceptre_object_high@functs_called |> as.logical(), rep(c(TRUE, FALSE), c(1, 6)))
  }
})
