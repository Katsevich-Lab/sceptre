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



test_that("run_calibration_check negative control pairs complement set", {
  grna_target_data_frame <- data.frame(
    grna_id = c("id1", "id2", "id3", "nt1", "nt2"),
    grna_target = c("t1", "t2", "t3", "non-targeting", "non-targeting"),
    chr = "", start = 0, end = 1
  )
  num_grna <- nrow(grna_target_data_frame)
  num_cells <- 50
  num_responses <- 10

  set.seed(1)
  # using sample(0:1) so no entries can accidentally cross the threshold
  grna_matrix <- matrix(sample(0:1, num_grna * num_cells, replace=TRUE), num_grna, num_cells) |>
    `rownames<-`(grna_target_data_frame$grna_id)
  cells_expressing_t1 <- 1:10
  cells_expressing_t2 <- 11:20
  cells_expressing_t3 <- 21:30
  cells_expressing_nt1 <- 31:40
  cells_expressing_nt2 <- 41:50
  all_cells <- 1:num_cells

  grna_matrix["id1", cells_expressing_t1] <- 50
  grna_matrix["id2", cells_expressing_t2] <- 50
  # grna_matrix["id3", cells_expressing_t3] <- 50
  grna_matrix["nt1", cells_expressing_nt1] <- 50
  grna_matrix["nt2", cells_expressing_nt2] <- 50

  response_matrix <- matrix(rpois( num_responses * num_cells, 1), num_responses, num_cells) |>
    `rownames<-`(c("t1", "t2", "t3", paste0("response_", 4:num_responses)))

  response_matrix["response_4", cells_expressing_t1] <- 100
  response_matrix["response_4", cells_expressing_nt1] <- 0
  # these two only matter for complement set
  response_matrix["response_4", cells_expressing_t2] <- 0
  response_matrix["response_4", cells_expressing_t3] <- 100

  response_matrix["response_5", cells_expressing_t2] <- 0
  response_matrix["response_5", all_cells[-cells_expressing_t2]] <- 100

  response_matrix["response_6", ] <- 100

  cells_to_remove_low_umi <- c(1,2,4,5,6,11,12,31)
  cells_to_remove_high_umi <- c(3, 13,32,33,34)
  response_matrix[,cells_to_remove_low_umi] <- 0
  response_matrix[,cells_to_remove_high_umi] <- 100000


  discovery_pairs <- data.frame(
    grna_target = c("t1",         "t2",         "t3"),
    response_id = c("response_4", "response_5", "response_6")
  )

  ## testing `control_group = "nt_cells"` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scep <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs,
      control_group = "complement"
    ) |>
    assign_grnas(method = "thresholding", threshold = 40) |>
    run_qc(response_n_umis_range = c(0, .90), response_n_nonzero_range = c(.15, 1),
           n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 11) |>
    run_calibration_check(calibration_group_size=1)

  # making sure the correct cells were removed
  remaining_cells <- all_cells[-c(cells_to_remove_low_umi, cells_to_remove_high_umi)]
  expect_setequal(scep@cells_in_use, remaining_cells)

  neg_df <- scep@negative_control_pairs
  for(i in 1:nrow(neg_df)) {
    trt_cells <- remaining_cells[scep@grna_assignments$indiv_nt_grna_idxs[[neg_df$grna_group[i]]]]
    expect_equal(
      neg_df$n_nonzero_trt[i],
      sum(response_matrix[neg_df$response_id[i], trt_cells] > 0)
    )

    # complement control group
    cntrl_cells <- setdiff(scep@cells_in_use, trt_cells)
    expect_equal(
      neg_df$n_nonzero_cntrl[i],
      sum(response_matrix[neg_df$response_id[i], cntrl_cells] > 0)
    )
  }
})



test_that("run_power_check", {
  grna_target_data_frame <- data.frame(
    grna_id = c("id1", "id2", "id3", "nt1"),
    grna_target = c("t1", "t2", "t3", "non-targeting"),
    chr = 0, start = 0, end = 1
  )
  num_grna <- nrow(grna_target_data_frame)
  num_cells <- 40
  num_responses <- 10

  set.seed(1)
  # using sample(0:1) so no entries can accidentally cross the threshold
  grna_matrix <- matrix(sample(0:1, num_grna * num_cells, replace=TRUE), num_grna, num_cells) |>
    `rownames<-`(grna_target_data_frame$grna_id)
  cells_expressing_t1 <- 1:10
  cells_expressing_t2 <- 11:20
  cells_expressing_t3 <- 21:30
  cells_expressing_nt1 <- 31:40
  all_cells <- 1:num_cells

  grna_matrix["id1", cells_expressing_t1] <- 50
  grna_matrix["id2", cells_expressing_t2] <- 50
  # grna_matrix["id3", cells_expressing_t3] <- 50
  grna_matrix["nt1", cells_expressing_nt1] <- 50

  response_matrix <- matrix(rpois( num_responses * num_cells, 1), num_responses, num_cells) |>
    `rownames<-`(c("t1", "t2", "t3", paste0("response_", 4:num_responses)))

  response_matrix["t1", cells_expressing_t1] <- 100 # should be highly significant

  response_matrix["t2", ] <- 100 # should not be significant at all

  positive_control_pairs = data.frame(
    grna_target = c("t1", "t2", "t3"),
    response_id = c("t1", "t2", "t3")
  )
  discovery_pairs <- data.frame(
    grna_target = c("t1",        "t1",        "t2",         "t2"),
    response_id = c("response_4", "response_5",  "response_4", "response_6")
  )

  ## testing `control_group = "nt_cells"` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scep_power <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      positive_control_pairs = positive_control_pairs,
      discovery_pairs = discovery_pairs,
      control_group = "nt_cells"
    ) |>
    assign_grnas(method = "thresholding", threshold = 40) |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1), # don't want to remove any cells here
           n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0) |>
    run_power_check()

  expect_equal(nrow(scep_power@power_result), nrow(positive_control_pairs))
  # this test is extremely significant
  expect_true(dplyr::filter(scep_power@power_result, response_id == "t1") |> dplyr::pull(p_value) < 1e-10)
  # this test is not significant at all
  expect_true(dplyr::filter(scep_power@power_result, response_id == "t2") |> dplyr::pull(p_value) > 0.5)
  # log fold change should be 0 here because the response_matrix values are constant
  expect_true(dplyr::filter(scep_power@power_result, response_id == "t3") |> dplyr::pull(log_2_fold_change) |> is.na())
})
