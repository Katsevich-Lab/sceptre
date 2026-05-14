# Exercise the parallel = TRUE code path on platforms where it is supported
# (Mac fork-based parallelism). Skipped on Windows because
# check_parallel_supported() refuses parallel = TRUE there. CI catches
# regressions that silently break parallel execution on macOS / Linux.

test_that("full pipeline runs end-to-end with parallel = TRUE", {
  skip_on_os("windows")

  data(highmoi_example_data)
  data(grna_target_data_frame_highmoi)

  scep <- import_data(
    response_matrix = highmoi_example_data$response_matrix,
    grna_matrix = highmoi_example_data$grna_matrix,
    grna_target_data_frame = grna_target_data_frame_highmoi,
    moi = "high",
    extra_covariates = highmoi_example_data$extra_covariates,
    response_names = highmoi_example_data$gene_names
  )
  positive_control_pairs <- construct_positive_control_pairs(scep)
  discovery_pairs <- construct_cis_pairs(
    scep,
    positive_control_pairs = positive_control_pairs,
    distance_threshold = 5e6
  )

  scep <- scep |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs,
      positive_control_pairs = positive_control_pairs
    ) |>
    assign_grnas(method = "thresholding") |>
    run_qc(
      response_n_umis_range = c(0, 1),
      response_n_nonzero_range = c(0, 1),
      n_nonzero_trt_thresh = 0,
      n_nonzero_cntrl_thresh = 0
    ) |>
    run_calibration_check(
      calibration_group_size = 1,
      parallel = TRUE,
      n_processors = 2
    ) |>
    run_power_check(parallel = TRUE, n_processors = 2) |>
    run_discovery_analysis(parallel = TRUE, n_processors = 2)

  expect_s4_class(scep, "sceptre_object")
  expect_true(scep@functs_called[["run_calibration_check"]])
  expect_true(scep@functs_called[["run_power_check"]])
  expect_true(scep@functs_called[["run_discovery_analysis"]])
  expect_gt(nrow(scep@calibration_result), 0)
  expect_gt(nrow(scep@power_result), 0)
  expect_gt(nrow(scep@discovery_result), 0)
})


test_that("assign_grnas mixture method runs with parallel = TRUE", {
  skip_on_os("windows")

  data(highmoi_example_data)
  data(grna_target_data_frame_highmoi)

  scep <- import_data(
    response_matrix = highmoi_example_data$response_matrix,
    grna_matrix = highmoi_example_data$grna_matrix,
    grna_target_data_frame = grna_target_data_frame_highmoi,
    moi = "high",
    extra_covariates = highmoi_example_data$extra_covariates,
    response_names = highmoi_example_data$gene_names
  )
  positive_control_pairs <- construct_positive_control_pairs(scep)
  discovery_pairs <- construct_cis_pairs(
    scep,
    positive_control_pairs = positive_control_pairs,
    distance_threshold = 5e6
  )

  scep <- scep |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs,
      positive_control_pairs = positive_control_pairs
    ) |>
    assign_grnas(method = "mixture", parallel = TRUE, n_processors = 2)
  expect_true(scep@functs_called[["assign_grnas"]])
  expect_true(length(scep@initial_grna_assignment_list) > 0)
})
