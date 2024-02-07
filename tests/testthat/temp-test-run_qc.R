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
  # for(i in 1:num_grnas) grna_matrix[i,i] <- 100 # for clear expression


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

# TODO finish controlling n_nonzero_trt and n_nonzero_cntrl for pass_qc in pos control
test_that("run_qc test positive control pairs", {
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
  response_matrix[] <- rpois(prod(dim(response_matrix)), 2)

  scep_high_all_pass <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      positive_control_pairs = test_data_list$positive_control_pairs,
      discovery_pairs = test_data_list$discovery_pairs,
      control_group = "nt_cells"
    ) |>
    assign_grnas(method = "thresholding", threshold = 7) |>
    run_qc()

  expect_true(all(scep_high_all_pass@positive_control_pairs_with_info$pass_qc))

  # TODO finish
  pos_ctrl_to_fail <- test_data_list$positive_control_pairs$grna_target[2]
  grna_ids_to_fail <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == pos_ctrl_to_fail])
  grna_matrix[grna_ids_to_fail,] <- 0
  grna_matrix[grna_ids_to_fail,1:2] <- 1  # not making it entirely 0

  scep_high_one_pos_fails <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      positive_control_pairs = test_data_list$positive_control_pairs,
      discovery_pairs = test_data_list$discovery_pairs,
      control_group = "nt_cells"
    ) |>
    assign_grnas(method = "thresholding", threshold = 7) |>
    run_qc()

  scep_high_one_pos_fails@positive_control_pairs_with_info


  setdiff(1:100, scep_high_one_pos_fails@cells_in_use)

})



test_that("run_qc test discovery pairs", {


})


