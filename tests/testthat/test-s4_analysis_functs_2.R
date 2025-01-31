test_that("run_calibration_check negative control pairs complement set with cellwise and pairwise qc", {
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
  grna_matrix <- matrix(sample(0:1, num_grna * num_cells, replace = TRUE), num_grna, num_cells) |>
    `rownames<-`(grna_target_data_frame$grna_id)
  cells_expressing_t1 <- 1:10
  cells_expressing_t2 <- 11:20
  cells_expressing_no_grna <- 21:30
  cells_expressing_nt1 <- 31:40
  cells_expressing_nt2 <- 41:50
  all_cells <- 1:num_cells

  grna_matrix["id1", cells_expressing_t1] <- 50
  grna_matrix["id2", cells_expressing_t2] <- 50
  grna_matrix["nt1", cells_expressing_nt1] <- 50
  grna_matrix["nt2", cells_expressing_nt2] <- 50

  response_matrix <- matrix(rpois(num_responses * num_cells, 1), num_responses, num_cells) |>
    `rownames<-`(c("t1", "t2", "t3", paste0("response_", 4:num_responses)))

  cells_to_remove_low_umi <- c(1, 2, 4, 5, 6, 11, 12, 31)
  cells_to_remove_high_umi <- c(3, 13, 32, 33, 34)
  response_matrix[, cells_to_remove_low_umi] <- 0
  response_matrix[, cells_to_remove_high_umi] <- 100000

  discovery_pairs <- data.frame(
    grna_target = c("t1", "t2", "t3"),
    response_id = c("response_4", "response_5", "response_6")
  )

  ## testing `control_group = "complement"` and `calibration_group_size=1` ~~~~~~~~~~~~~~~~~~~~~
  scep_pre <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs,
      control_group = "complement"
    ) |>
    assign_grnas(method = "thresholding", threshold = 40)

  # set.seed(5)
  scep_complement_size_1 <- scep_pre |>
    run_qc(
      response_n_umis_range = c(0, .90), response_n_nonzero_range = c(.15, 1),
      # with `n_nonzero_cntrl_thresh = 17` one discovery pair fails
      n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 17
    ) |>
    run_calibration_check(calibration_group_size = 1)

  # making sure the correct cells were removed
  remaining_cells <- all_cells[-c(cells_to_remove_low_umi,
                                  cells_to_remove_high_umi, cells_expressing_no_grna)]
  expect_setequal(scep_complement_size_1@cells_in_use, remaining_cells)

  neg_df <- scep_complement_size_1@negative_control_pairs

  # a single neg ctrl pair with this seed fails QC so there are 2 left
  expect_equal(nrow(neg_df), sum(scep_complement_size_1@discovery_pairs_with_info$pass_qc))

  expect_true(all(neg_df$pass_qc))

  for (i in 1:nrow(neg_df)) {
    trt_cells <- remaining_cells[scep_complement_size_1@grna_assignments$indiv_nt_grna_idxs[[neg_df$grna_group[i]]]]
    expect_equal(
      neg_df$n_nonzero_trt[i],
      sum(response_matrix[neg_df$response_id[i], trt_cells] > 0)
    )

    # complement control group
    cntrl_cells <- setdiff(scep_complement_size_1@cells_in_use, trt_cells)
    expect_equal(
      neg_df$n_nonzero_cntrl[i],
      sum(response_matrix[neg_df$response_id[i], cntrl_cells] > 0)
    )
  }

  ## testing `control_group = "complement"` and `calibration_group_size=2` ~~~~~~~~~~~~~~~~~~~~~
  set.seed(2)
  scep_complement_size_2 <- scep_pre |>
    run_qc(
      response_n_umis_range = c(0, .90), response_n_nonzero_range = c(.15, 1),
      # with `n_nonzero_trt_thresh = 1` a single discovery pair fails pairwise QC
      n_nonzero_trt_thresh = 1, n_nonzero_cntrl_thresh = 0
    ) |>
    run_calibration_check(calibration_group_size = 2)

  # making sure the correct cells were removed
  remaining_cells <- all_cells[-c(cells_to_remove_low_umi,
                                  cells_to_remove_high_umi, cells_expressing_no_grna)]
  expect_setequal(scep_complement_size_2@cells_in_use, remaining_cells)

  neg_df <- scep_complement_size_2@negative_control_pairs

  # with this seed we lose one pair but its due to pairwise QC on the actual discovery pairs
  expect_equal(nrow(neg_df), sum(scep_complement_size_2@discovery_pairs_with_info$pass_qc))
  expect_true(all(neg_df$pass_qc))

  for (i in 1:nrow(neg_df)) {
    # all nt cells each time
    trt_cells <- remaining_cells[scep_complement_size_2@grna_assignments$indiv_nt_grna_idxs |> unlist()]
    expect_equal(
      neg_df$n_nonzero_trt[i],
      sum(response_matrix[neg_df$response_id[i], trt_cells] > 0)
    )

    # complement control group
    cntrl_cells <- setdiff(scep_complement_size_2@cells_in_use, trt_cells)
    expect_equal(
      neg_df$n_nonzero_cntrl[i],
      sum(response_matrix[neg_df$response_id[i], cntrl_cells] > 0)
    )
  }
})

test_that("run_calibration_check negative control pairs nt set with cellwise qc", {
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
  grna_matrix <- matrix(sample(0:1, num_grna * num_cells, replace = TRUE), num_grna, num_cells) |>
    `rownames<-`(grna_target_data_frame$grna_id)
  cells_expressing_t1 <- 1:10
  cells_expressing_t2 <- 11:20
  cells_expressing_no_grna <- 21:30
  cells_expressing_nt1 <- 31:40
  cells_expressing_nt2 <- 41:50
  all_cells <- 1:num_cells

  grna_matrix["id1", cells_expressing_t1] <- 50
  grna_matrix["id2", cells_expressing_t2] <- 50
  grna_matrix["nt1", cells_expressing_nt1] <- 50
  grna_matrix["nt2", cells_expressing_nt2] <- 50

  response_matrix <- matrix(rpois(num_responses * num_cells, 1), num_responses, num_cells) |>
    `rownames<-`(c("t1", "t2", "t3", paste0("response_", 4:num_responses)))

  cells_to_remove_low_umi <- c(1, 2, 4, 5, 6, 11, 12, 31)
  cells_to_remove_high_umi <- c(3, 13, 32, 33, 34)
  response_matrix[, cells_to_remove_low_umi] <- 0
  response_matrix[, cells_to_remove_high_umi] <- 100000

  discovery_pairs <- data.frame(
    grna_target = c("t1", "t2", "t3"),
    response_id = c("response_4", "response_5", "response_6")
  )

  ## testing `control_group = "complement"` and `calibration_group_size=1` ~~~~~~~~~~~~~~~~~~~~~
  scep_pre <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs,
      control_group = "nt_cells"
    ) |>
    assign_grnas(method = "thresholding", threshold = 40)

  scep_nt_size_1 <- scep_pre |>
    run_qc(
      response_n_umis_range = c(0, .90), response_n_nonzero_range = c(.15, 1),
      # with `n_nonzero_cntrl_thresh = 20` one discovery pair fails
      n_nonzero_trt_thresh = 1, n_nonzero_cntrl_thresh = 0
    ) |>
    run_calibration_check(calibration_group_size = 1)

  # making sure the correct cells were removed
  remaining_cells <- all_cells[-c(cells_to_remove_low_umi,
                                  cells_to_remove_high_umi, cells_expressing_no_grna)]
  expect_setequal(scep_nt_size_1@cells_in_use, remaining_cells)

  neg_df <- scep_nt_size_1@negative_control_pairs

  # a single neg ctrl pair with this seed fails QC so there are 2 left
  expect_equal(nrow(neg_df), sum(scep_nt_size_1@discovery_pairs_with_info$pass_qc))
  expect_true(all(neg_df$pass_qc))

  for (i in 1:nrow(neg_df)) {
    trt_cells <- remaining_cells[scep_nt_size_1@grna_assignments$all_nt_idxs[scep_nt_size_1@grna_assignments$indiv_nt_grna_idxs[[neg_df$grna_group[i]]]]]
    expect_equal(
      neg_df$n_nonzero_trt[i],
      sum(response_matrix[neg_df$response_id[i], trt_cells] > 0)
    )
    # nt control group
    other_nt <- setdiff(c("nt1", "nt2"), neg_df$grna_group[i])

    cntrl_cells <- remaining_cells[scep_nt_size_1@grna_assignments$all_nt_idxs[scep_nt_size_1@grna_assignments$indiv_nt_grna_idxs[[other_nt]]]]
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
  grna_matrix <- matrix(sample(0:1, num_grna * num_cells, replace = TRUE), num_grna, num_cells) |>
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

  response_matrix <- matrix(rpois(num_responses * num_cells, 1), num_responses, num_cells) |>
    `rownames<-`(c("t1", "t2", "t3", paste0("response_", 4:num_responses)))

  response_matrix["t1", cells_expressing_t1] <- 100 # should be highly significant

  response_matrix["t2", ] <- 100 # should not be significant at all

  positive_control_pairs <- data.frame(
    grna_target = c("t1", "t2", "t3"),
    response_id = c("t1", "t2", "t3")
  )
  discovery_pairs <- data.frame(
    grna_target = c("t1", "t1", "t2", "t2"),
    response_id = c("response_4", "response_5", "response_4", "response_6")
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
    run_qc(
      response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1), # don't want to remove any cells here
      n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0
    ) |>
    run_power_check()

  expect_equal(nrow(scep_power@power_result), nrow(positive_control_pairs))
  # this test is extremely significant
  expect_true(dplyr::filter(scep_power@power_result, response_id == "t1") |> dplyr::pull(p_value) < 1e-10)
  # this test is not significant at all
  expect_true(dplyr::filter(scep_power@power_result, response_id == "t2") |> dplyr::pull(p_value) > 0.5)
  # log fold change should be NA here because the response_matrix values are constant
  expect_true(dplyr::filter(scep_power@power_result, response_id == "t3") |> dplyr::pull(log_2_fold_change) |> is.na())
})

test_that("run_discovery_analysis", {
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
  grna_matrix <- matrix(sample(0:1, num_grna * num_cells, replace = TRUE), num_grna, num_cells) |>
    `rownames<-`(grna_target_data_frame$grna_id)
  cells_expressing_t1 <- 1:10
  cells_expressing_t2 <- 11:20
  cells_expressing_t3 <- 21:30
  cells_expressing_nt1 <- 31:40
  cells_expressing_nt2 <- 41:50
  all_cells <- 1:num_cells

  grna_matrix["id1", cells_expressing_t1] <- 50
  grna_matrix["id2", cells_expressing_t2] <- 50
  grna_matrix["id3", cells_expressing_t3] <- 50
  grna_matrix["nt1", cells_expressing_nt1] <- 50
  grna_matrix["nt2", cells_expressing_nt2] <- 50

  response_matrix <- matrix(rpois(num_responses * num_cells, 1), num_responses, num_cells) |>
    `rownames<-`(c("t1", "t2", "t3", paste0("response_", 4:num_responses)))

  response_matrix["t1", cells_expressing_t1] <- 100 # should be highly significant

  response_matrix["t2", ] <- 100 # should not be significant at all

  positive_control_pairs <- data.frame(
    grna_target = c("t1", "t2", "t3"),
    response_id = c("t1", "t2", "t3")
  )
  discovery_pairs <- data.frame(
    grna_target = c("t1", "t1", "t2", "t2"),
    response_id = c("response_4", "response_5", "response_4", "response_6")
  )

  ## testing `control_group = "nt_cells"` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  disc_results <- import_data(
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
    # don't want to remove any cells for this one
    run_qc(
      response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
      # making the response_5, t1 pair fail pairwise QC
      n_nonzero_trt_thresh = 7, n_nonzero_cntrl_thresh = 0
    ) |>
    run_calibration_check(calibration_group_size = 1) |>
    run_power_check() |>
    run_discovery_analysis() |>
    get_result(analysis = "run_discovery_analysis")

  # confirming everything is ok with the pair failing pairwise QC
  expect_false(disc_results[disc_results$grna_target == "t1" & disc_results$response_id == "response_5", "pass_qc"][[1]])

  # t1 should be significant
  expect_true(disc_results[disc_results$grna_target == "t1" & disc_results$response_id == "response_4", "significant"][[1]])

  # t2 tests are not significant
  expect_false(any(disc_results[disc_results$grna_target == "t2", "significant"]))

  # confirm that n_trt and n_cntrl are calculated correctly
  expect_true(all(disc_results$n_trt == 10))
  expect_true(all(disc_results$n_cntrl == 20))

  ## testing `control_group = "complement"` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  disc_results <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      positive_control_pairs = positive_control_pairs,
      discovery_pairs = discovery_pairs,
      control_group = "complement"
    ) |>
    assign_grnas(method = "thresholding", threshold = 40) |>
    # don't want to remove any cells for this one
    run_qc(
      response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
      n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0
    ) |>
    run_discovery_analysis() |>
    get_result(analysis = "run_discovery_analysis")

  # confirm that n_trt and n_cntrl are calculated correctly
  expect_true(all(disc_results$n_trt == 10))
  expect_true(all(disc_results$n_cntrl == 40))
})
