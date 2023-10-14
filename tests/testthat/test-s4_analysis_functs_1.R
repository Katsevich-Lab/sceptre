
# most of the content of `import_data` is subfunctions
# which are tested elsewhere
test_that("import_data", {
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

  expect_equal(sceptre_object_low_no_ec_with_response_names@covariate_data_frame |> rownames(),
               as.character(seq_len(num_cells)))

  expect_equal(
    names(sceptre_object_low_no_ec_with_response_names@covariate_data_frame),
    c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis")
  )
  expect_equal(
    names(sceptre_object_high_with_ec@covariate_data_frame),
    c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis", "x")
  )
  expect_equal(
    sceptre_object_high_with_ec@covariate_data_frame$x,
    extra_covariates$x
  )

  ##### slots set in section 5 of `import_data`

  expect_equal(sceptre_object_high_with_ec@last_function_called, "import_data")
  expect_equal(names(sceptre_object_high_with_ec@functs_called), c("import_data", "set_analysis_parameters",
                                                                   "assign_grnas", "run_qc", "run_calibration_check",
                                                                   "run_power_check", "run_discovery_analysis"))
  expect_equal(sceptre_object_high_with_ec@functs_called |> as.logical(), rep(c(TRUE, FALSE), c(1, 6)))
})
