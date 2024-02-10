


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
