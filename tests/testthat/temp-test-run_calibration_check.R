
test_that("run_calibration_check", {
  set.seed(2)

  num_cells <- 100
  num_responses <- 40

  grna_target_data_frame <- make_mock_grna_target_data(c(2,2,2), 1, 1, 10)
  grna_matrix <- rpois(num_cells * nrow(grna_target_data_frame), 5) |>
    matrix(nrow = nrow(grna_target_data_frame), ncol = num_cells) |>
    `rownames<-`(grna_target_data_frame$grna_id)

  response_matrix <- rpois(num_cells * num_responses, 5) |>
    matrix(nrow = num_responses, ncol = num_cells) |>
    `rownames<-`(paste0("response_", 1:num_responses))

  discovery_pairs <- data.frame(
    grna_target = "t1_c1_d1",
    response_id = rownames(response_matrix)[1:10]
  )

  scep_pre_calib <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs,
      control_group = "nt_cells"
    ) |>
    assign_grnas(method = "thresholding", threshold = 5) |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1),
           n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0)

  scep_calib_1 <- scep_pre_calib |>
    run_calibration_check(calibration_group_size=1)
  scep_calib_3 <- scep_pre_calib |>
    run_calibration_check(calibration_group_size=3)

  expect_equal(nrow(scep_calib_1@calibration_result), nrow(discovery_pairs))
  expect_equal(nrow(scep_calib_3@calibration_result), nrow(discovery_pairs))

  expect_false(any(grepl(pattern="&", x=scep_calib_1@calibration_result$grna_target, fixed = TRUE)))
  expect_equal(strsplit(scep_calib_3@calibration_result$grna_target, "&") |> sapply(length), rep(3, nrow(discovery_pairs)))

  # with this seed all nulls are false for both objects
  expect_false(any(scep_calib_1@calibration_result$significant))
  expect_false(any(scep_calib_3@calibration_result$significant))
})
