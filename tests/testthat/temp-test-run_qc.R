# NOTES
# - `run_qc()` fails if there are no cells with NT grnas expressed, and the error is not helpful:
#     "Error in seq.default(start[i], stop[i]) : 'to' must be a finite number"
# -


# TODO
# - confirm grna assignments updated correctly

make_mock_base_data_for_testing_run_qc <- function(num_cells) {
  # num_cells <- 24  # must be even for my choice of extra_covariates
  num_responses <- 7
  grna_target_data_frame <- make_mock_grna_target_data(c(2,2), 1, 1, 3)
  on_targets <- unique(grna_target_data_frame$grna_target)[1:2]
  num_grnas <- nrow(grna_target_data_frame)

  grna_matrix <- make_mock_grna_matrices(
    grna_target_data_frame, non_nt_patterns="zero",
    nt_patterns = "zero", num_cells = num_cells
  )

  response_matrix <- make_mock_response_matrices(
    num_responses = num_responses, num_cells = num_cells,
    patterns = "column"
  ) |>
    `rownames<-`(c(on_targets, paste0("response_", (length(on_targets) + 1):num_responses)))

  extra_covariates <- data.frame(x = rep(c("b1", "b2"), each = num_cells / 2))

  positive_control_pairs <- data.frame(
    grna_target = on_targets,
    response_id = on_targets
  )

  discovery_pairs <- data.frame(
    grna_target = rep(on_targets, times = 2),
    response_id = rep(paste0("response_", length(on_targets) + 1:2), each = 2)
  )

  list(
    grna_target_data_frame = grna_target_data_frame,
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    extra_covariates = extra_covariates,
    positive_control_pairs = positive_control_pairs,
    discovery_pairs = discovery_pairs
  )
}

test_that("run_qc remove cells with multiple grnas in low moi", {

  test_data_list <- make_mock_base_data_for_testing_run_qc(num_cells = 24)
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  set.seed(123)
  grna_matrix <- test_data_list$grna_matrix
  grna_matrix[] <- rpois(prod(dim(grna_matrix)), 1)
  for(i in 1:num_grnas) grna_matrix[i,i] <- 10 # for clear expression
  locs_of_extra_expression <- 2:4
  grna_matrix[1, locs_of_extra_expression] <- 10  # these cells now have 2 grnas clearly expressed

  response_matrix <- test_data_list$response_matrix
  response_matrix[] <- rep(c(0,0,1,0), times=prod(dim(response_matrix))/4)
  for(i in 1:nrow(response_matrix)) response_matrix[i,i] <- 2 # to avoid low rank covariate data frame

  scep_low <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      positive_control_pairs = test_data_list$positive_control_pairs,
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "thresholding", threshold = 9) |>
    run_qc()

  expect_equal(scep_low@cells_in_use, (1:num_cells)[-locs_of_extra_expression])
})


test_that("run_qc remove cells via `additional_cells_to_remove` high moi", {

  test_data_list <- make_mock_base_data_for_testing_run_qc(num_cells = 24)
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  set.seed(123)
  grna_matrix <- test_data_list$grna_matrix
  grna_matrix[] <- rpois(prod(dim(grna_matrix)), 2)
  for(i in 1:num_grnas) grna_matrix[i,i] <- 10 # for clear expression

  response_matrix <- test_data_list$response_matrix
  response_matrix[] <- rep(c(0,0,1,0), times=prod(dim(response_matrix))/4)
  for(i in 1:nrow(response_matrix)) response_matrix[i,i] <- 2 # to avoid low rank covariate data frame

  scep_high_pre_qc <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      positive_control_pairs = test_data_list$positive_control_pairs,
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "thresholding", threshold = 5)

  expect_equal(run_qc(scep_high_pre_qc)@cells_in_use, 1:num_cells) # all here

  cells_to_remove <- c(1,3,5)

  scep_qc <- run_qc(scep_high_pre_qc, additional_cells_to_remove = cells_to_remove)
  expect_equal(scep_qc@cells_in_use, (1:num_cells)[-cells_to_remove])
})


test_that("run_qc remove cells with `response_n_umis_range` and `response_n_nonzero_range` high moi", {

  test_data_list <- make_mock_base_data_for_testing_run_qc(num_cells = 100)
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  set.seed(123)
  grna_matrix <- test_data_list$grna_matrix
  grna_matrix[] <- rpois(prod(dim(grna_matrix)), 4)

  response_matrix <- test_data_list$response_matrix
  response_matrix[] <- rpois(prod(dim(response_matrix)), 2) + 10
  zero_umi_count_idx <- c(1,5)
  high_umi_count_idx <- 2:4
  response_matrix[,zero_umi_count_idx] <- 0
  response_matrix[,high_umi_count_idx] <- 100

  scep_high_pre_qc <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      positive_control_pairs = test_data_list$positive_control_pairs,
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "thresholding", threshold = 7)

  scep_all <- scep_high_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1))
  expect_identical(scep_all@cells_in_use, 1:100)  # all remain currently

  scep_qc <- scep_high_pre_qc |>
    run_qc(response_n_umis_range = c(0.02, 1), response_n_nonzero_range = c(0,1))
  expect_identical(scep_qc@cells_in_use, (1:100)[-zero_umi_count_idx])

  scep_qc <- scep_high_pre_qc |>
    run_qc(response_n_umis_range = c(0, .97), response_n_nonzero_range = c(0,1))
  expect_identical(scep_qc@cells_in_use, (1:100)[-high_umi_count_idx])

  scep_qc <- scep_high_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0.03,1))
  expect_identical(scep_qc@cells_in_use, (1:100)[-zero_umi_count_idx])

  scep_qc <- scep_high_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,.02))
  expect_equal(scep_qc@cells_in_use, zero_umi_count_idx)
})

test_that("run_qc remove cells with `p_mito_threshold` high moi", {

  test_data_list <- make_mock_base_data_for_testing_run_qc(num_cells = 100)
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  set.seed(123)
  grna_matrix <- test_data_list$grna_matrix
  grna_matrix[] <- rpois(prod(dim(grna_matrix)), 4)

  response_matrix <- test_data_list$response_matrix
  rownames(response_matrix)[6:7] <- paste0("MT-", rownames(response_matrix)[6:7])
  response_matrix[] <- rpois(prod(dim(response_matrix)), 2)
  high_p_mito_idx <- c(1,5)
  response_matrix[6:7,high_p_mito_idx] <- 100

  scep_high_pre_qc <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      positive_control_pairs = test_data_list$positive_control_pairs,
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "thresholding", threshold = 7)

  scep_qc <- scep_high_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1), p_mito_threshold = .93)

  expect_identical(scep_qc@cells_in_use, (1:100)[-high_p_mito_idx])
})

test_that("run_qc test positive control pairs", {
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

  # response_matrix <- matrix(sample(0:1, num_responses * num_cells, replace=TRUE, prob = c(.9, .1)), num_responses, num_cells) |>
  #   `rownames<-`(c("t1", "t2", "t3", paste0("response_", 4:num_responses)))
  response_matrix <- matrix(rpois( num_responses * num_cells, 1), num_responses, num_cells) |>
    `rownames<-`(c("t1", "t2", "t3", paste0("response_", 4:num_responses)))

  # # target t1: all control group are non-zero, all NT are non-zero
  # # whether or not the counts reflect this will depend on if there are "control" cells
  # # other than `cells_expressing_nt1`
  response_matrix["t1", cells_expressing_t1] <- 100
  response_matrix["t1", cells_expressing_nt1] <- 100
  # these next two only matter for complement control group: this makes it so
  # there will be 20 non-zero control cells in that situation
  response_matrix["t1", cells_expressing_t2] <- 0
  response_matrix["t1", cells_expressing_t3] <- 100

  # target t2: all control group are non-zero, all NT are 0
  response_matrix["t2", cells_expressing_t2] <- 100
  response_matrix["t2", cells_expressing_nt1] <- 0
  # these next two only matter for complement control group: this makes it so
  # there will be 0 non-zero control cells in that situation
  response_matrix["t2", cells_expressing_t1] <- 0
  response_matrix["t2", cells_expressing_t3] <- 0

  # target t3: all cells are non-zero
  response_matrix["t3", all_cells] <- 100

  positive_control_pairs = data.frame(
    grna_target = c("t1", "t2", "t3"),
    response_id = c("t1", "t2", "t3")
  )
  discovery_pairs <- data.frame(
    grna_target = c("t1",        "t1",        "t2",         "t2"),
    response_id = c("response_4", "response_5",  "response_4", "response_6")
  )

  ## testing `control_group = "nt_cells"` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scep_low_control_nt_pre_qc <- import_data(
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
    assign_grnas(method = "thresholding", threshold = 40)

  scep_low_control_nt_all_pass <- scep_low_control_nt_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1), # don't want to remove any cells here
           n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0)

  # confirming the values of `n_nonzero_trt` and `n_nonzero_cntrl`
  expect_equal(scep_low_control_nt_all_pass@positive_control_pairs_with_info$n_nonzero_trt, c(10, 10, 0))
  expect_equal(scep_low_control_nt_all_pass@positive_control_pairs_with_info$n_nonzero_cntrl, c(10, 0, 10))

  # confirming all pass QC with these thresholds
  expect_true(all(scep_low_control_nt_all_pass@positive_control_pairs_with_info$pass_qc))

  # for this one only t3 fails because it has no cells expressing the response, since no gRNAs express t3
  scep_low_control_nt_t3_fail <- scep_low_control_nt_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1), # don't want to remove any cells here
           n_nonzero_trt_thresh = 1, n_nonzero_cntrl_thresh = 0)
  expect_equal(scep_low_control_nt_t3_fail@positive_control_pairs_with_info$pass_qc, c(TRUE, TRUE, FALSE))

  # for this one only t2 fails, because that's the only one with 0 control group responses
  scep_low_control_nt_t2_fail <- scep_low_control_nt_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1), # don't want to remove any cells here
           n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 1)
  expect_equal(scep_low_control_nt_t2_fail@positive_control_pairs_with_info$pass_qc, c(TRUE, FALSE, TRUE))


  ## testing `control_group = "complement"` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scep_high_control_complement_pre_qc <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "high" # forces complement control group
  ) |>
    set_analysis_parameters(
      positive_control_pairs = positive_control_pairs,
      discovery_pairs = discovery_pairs,
    ) |>
    assign_grnas(method = "thresholding", threshold = 40)

  scep_high_control_complement_all_pass <- scep_high_control_complement_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1), # don't want to remove any cells here
           n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0)

  # confirming the values of `n_nonzero_trt` and `n_nonzero_cntrl`
  expect_equal(scep_high_control_complement_all_pass@positive_control_pairs_with_info$n_nonzero_trt, c(10, 10, 0))
  expect_equal(scep_high_control_complement_all_pass@positive_control_pairs_with_info$n_nonzero_cntrl, c(20, 0, 40))

  # confirming all pass QC with these thresholds
  expect_true(all(scep_high_control_complement_all_pass@positive_control_pairs_with_info$pass_qc))

  scep_high_control_complement_t3_fail <- scep_high_control_complement_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1),  # don't want to remove any cells here
           n_nonzero_trt_thresh = 1, n_nonzero_cntrl_thresh = 0)
  expect_equal(scep_high_control_complement_t3_fail@positive_control_pairs_with_info$pass_qc, c(TRUE, TRUE, FALSE))

  scep_high_control_complement_t1_t2_fail <- scep_high_control_complement_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1), # don't want to remove any cells here
           n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 21)
  expect_equal(scep_high_control_complement_t1_t2_fail@positive_control_pairs_with_info$pass_qc, c(FALSE, FALSE, TRUE))
})



test_that("run_qc test discovery pairs", {
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

  response_matrix["response_4", cells_expressing_t1] <- 100
  response_matrix["response_4", cells_expressing_nt1] <- 0
  # these two only matter for complement set
  response_matrix["response_4", cells_expressing_t2] <- 0
  response_matrix["response_4", cells_expressing_t3] <- 100

  response_matrix["response_5", cells_expressing_t2] <- 0
  response_matrix["response_5", all_cells[-cells_expressing_t2]] <- 100

  response_matrix["response_6", ] <- 100

  discovery_pairs <- data.frame(
    grna_target = c("t1",         "t2",         "t3"),
    response_id = c("response_4", "response_5", "response_6")
  )

  ## testing `control_group = "nt_cells"` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scep_low_control_nt_pre_qc <- import_data(
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


  scep_low_control_nt_all_pass <- scep_low_control_nt_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1),
           n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0)  # don't want to remove the cells I'm messing with

  # confirming the values of `n_nonzero_trt` and `n_nonzero_cntrl`
  expect_equal(scep_low_control_nt_all_pass@discovery_pairs_with_info$n_nonzero_trt, c(10, 0, 0))
  expect_equal(scep_low_control_nt_all_pass@discovery_pairs_with_info$n_nonzero_cntrl, c(0, 10, 10))

  # confirming all pass QC with these thresholds
  expect_true(all(scep_low_control_nt_all_pass@discovery_pairs_with_info$pass_qc))

  scep_low_control_nt_t2_t3_fail <- scep_low_control_nt_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1),
           n_nonzero_trt_thresh = 1, n_nonzero_cntrl_thresh = 0)  # don't want to remove the cells I'm messing with
  expect_equal(scep_low_control_nt_t2_t3_fail@discovery_pairs_with_info$pass_qc, c(TRUE, FALSE, FALSE))

  scep_low_control_nt_t1_fail <- scep_low_control_nt_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1),
           n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 1)  # don't want to remove the cells I'm messing with
  expect_equal(scep_low_control_nt_t1_fail@discovery_pairs_with_info$pass_qc, c(FALSE, TRUE, TRUE))

  ## testing `control_group = "complement"` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scep_high_control_complement_pre_qc <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "high" # forces complement control group
  ) |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs,
    ) |>
    assign_grnas(method = "thresholding", threshold = 40)

  scep_high_control_complement_all_pass <- scep_high_control_complement_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1),
           n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0)  # don't want to remove the cells I'm messing with

  # confirming the values of `n_nonzero_trt` and `n_nonzero_cntrl`
  expect_equal(scep_high_control_complement_all_pass@discovery_pairs_with_info$n_nonzero_trt, c(10, 0, 0))
  expect_equal(scep_high_control_complement_all_pass@discovery_pairs_with_info$n_nonzero_cntrl, c(10, 30, 40))

  # confirming all pass QC with these thresholds
  expect_true(all(scep_high_control_complement_all_pass@discovery_pairs_with_info$pass_qc))

  scep_high_control_complement_t2_t3_fail <- scep_high_control_complement_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1),
           n_nonzero_trt_thresh = 1, n_nonzero_cntrl_thresh = 0)  # don't want to remove the cells I'm messing with
  expect_equal(scep_high_control_complement_t2_t3_fail@discovery_pairs_with_info$pass_qc, c(TRUE, FALSE, FALSE))

  scep_high_control_complement_t1_fail <- scep_high_control_complement_pre_qc |>
    run_qc(response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0,1),
           n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 21)  # don't want to remove the cells I'm messing with
  expect_equal(scep_high_control_complement_t1_fail@discovery_pairs_with_info$pass_qc, c(FALSE, TRUE, TRUE))


})


