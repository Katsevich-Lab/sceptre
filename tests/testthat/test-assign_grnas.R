# mock data function just for the tests in this file
make_mock_base_data_for_testing_assign_grnas <- function() {
  num_cells <- 24 # must be even for my choice of extra_covariates
  num_responses <- 7
  grna_target_data_frame <- make_mock_grna_target_data(c(2, 2), 1, 1, 3)
  on_targets <- unique(grna_target_data_frame$grna_target)[1:2]
  num_grnas <- nrow(grna_target_data_frame)

  grna_matrix <- make_mock_grna_matrices(
    grna_target_data_frame,
    non_nt_patterns = "zero",
    nt_patterns = "zero", num_cells = num_cells
  )

  response_matrix <- make_mock_response_matrices(
    num_responses = num_responses, num_cells = num_cells,
    patterns = "column"
  ) |>
    `rownames<-`(c(on_targets, paste0("response_", (length(on_targets) + 1):num_responses)))

  discovery_pairs <- data.frame(
    grna_target = on_targets,
    response_id = rep(paste0("response_", length(on_targets) + 1), times = 2)
  )

  list(
    grna_target_data_frame = grna_target_data_frame,
    response_matrix = response_matrix,
    grna_matrix_all_0 = grna_matrix,
    discovery_pairs = discovery_pairs
  )
}

test_that("assign_grnas method=maximum moi=low grna_matrix all 1", {
  test_data_list <- make_mock_base_data_for_testing_assign_grnas()
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix_all_0)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  grna_matrix_all_1 <- test_data_list$grna_matrix_all_0
  grna_matrix_all_1[] <- 1

  scep_low_all_1 <- import_data(
    grna_matrix = grna_matrix_all_1,
    response_matrix = test_data_list$response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "maximum")

  # every cell here should look like it has multiple GRNAs
  expect_equal(scep_low_all_1@cells_w_zero_or_twoplus_grnas, 1:num_cells)

  # first element should have all idx
  expect_equal(
    scep_low_all_1@initial_grna_assignment_list,
    # every cell is assigned to the first grna_id that appears since they are all tied
    # plus all of these cells will be removed anyway in qc so the assignment doesn't matter
    lapply(test_data_list$grna_target_data_frame$grna_id, function(target_name) {
      if (target_name == test_data_list$grna_target_data_frame$grna_id[1]) 1:num_cells else integer(0)
    }) |>
      setNames(test_data_list$grna_target_data_frame$grna_id)
  )

  # first element should have all idx

  guide_part_all_1 <- lapply(unique_targets, function(target_name) {
    if (target_name == test_data_list$grna_target_data_frame$grna_target[1]) 1:num_cells else integer(0)
  }) |>
    setNames(unique_targets)
  nt_part_all_1 <- lapply(nt_guides, function(nt_guide) integer(0)) |>
    setNames(nt_guides)

  expect_equal(
    scep_low_all_1@grna_assignments_raw,
    list(
      grna_group_idxs = guide_part_all_1,
      indiv_nt_grna_idxs = nt_part_all_1,
      all_nt_idxs = integer(0)
    )
  )
  # should be empty for maximum assignment
  expect_equal(scep_low_all_1@grnas_per_cell, integer(0))
})

test_that("assign_grnas method=maximum moi=low grna_matrix clear max", {
  test_data_list <- make_mock_base_data_for_testing_assign_grnas()
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix_all_0)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  grna_matrix_clear_max <- test_data_list$grna_matrix_all_0
  for (i in 1:num_grnas) grna_matrix_clear_max[i, i] <- 100

  scep_low_clear_max <- import_data(
    grna_matrix = grna_matrix_clear_max,
    response_matrix = test_data_list$response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "maximum")

  # expected answer: cell i expresses grna i
  expected_initial_assignment_list <- lapply(1:num_grnas, function(i) {
    if(i == 1) c(1, (num_grnas + 1):num_cells) else i
  })  |>
    setNames(test_data_list$grna_target_data_frame$grna_id)
  expect_equal(
    scep_low_clear_max@initial_grna_assignment_list,
    expected_initial_assignment_list
  )
  # all cells from idx `num_grnas+1` onward have no UMI counts at all, so they
  # all are considered to have multiple grnas
  expect_equal(scep_low_clear_max@cells_w_zero_or_twoplus_grnas, (num_grnas + 1):num_cells)

  expect_equal(
    scep_low_clear_max@grna_assignments_raw,
    list(
      grna_group_idxs = list(
        t1_c1_d1 = c(1L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 21L, 22L, 23L, 24L, 2L),
        t1_c2_d1 = 3:4, t1_c3_d1 = 5:6, t2_c3_d1 = 7:8
      ),
      indiv_nt_grna_idxs = list(nt1 = 9L, nt2 = 10L, nt3 = 11L),
      all_nt_idxs = c(9L, 10L, 11L)
    )
  )
  # should be empty for maximum assignment
  expect_equal(scep_low_clear_max@grnas_per_cell, integer(0))
})

test_that("assign_grnas method=maximum moi=low grna_matrix varying umi_fraction_threshold", {
  test_data_list <- make_mock_base_data_for_testing_assign_grnas()
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix_all_0)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  grna_matrix_vary_frac <- test_data_list$grna_matrix_all_0
  grna_matrix_vary_frac[1, ] <- 100 # all cells express grna1 very strongly
  grna_matrix_vary_frac[2, 1] <- 25 # cell 1 also expresses some grna2
  # in total, for cell 1 we have grna1 as exactly 80% of the total UMI count

  # cell 1 IS NOT flagged as expressing 2 grna with low threshold
  scep_low_with_low_frac <- import_data(
    grna_matrix = grna_matrix_vary_frac,
    response_matrix = test_data_list$response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "maximum", umi_fraction_threshold = .8 - .001)

  expect_equal(scep_low_with_low_frac@cells_w_zero_or_twoplus_grnas, integer(0))

  # cell 1 IS flagged as expressing 2 grna with high threshold
  scep_low_with_high_frac <- import_data(
    grna_matrix = grna_matrix_vary_frac,
    response_matrix = test_data_list$response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "maximum", umi_fraction_threshold = .8 + .001)

  expect_equal(scep_low_with_high_frac@cells_w_zero_or_twoplus_grnas, 1)
})

test_that("assign_grnas method=threshold moi=low", {
  test_data_list <- make_mock_base_data_for_testing_assign_grnas()
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix_all_0)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  grna_matrix_vary_thresh <- test_data_list$grna_matrix_all_0
  expressed_grna_ids <- c(5, 6)
  grna_matrix_vary_thresh[expressed_grna_ids[1], ] <- 100 # all cells express this grna_id very strongly
  grna_matrix_vary_thresh[expressed_grna_ids[2], ] <- 25 # all cells also express some of this grna_id

  # every cell is flagged as expressing 2 grnas ~~~~~~~~~~~~~~~~~~~~~~~~~
  scep_low_with_low_thresh <- import_data(
    grna_matrix = grna_matrix_vary_thresh,
    response_matrix = test_data_list$response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "thresholding", threshold = 25)

  expect_equal(scep_low_with_low_thresh@cells_w_zero_or_twoplus_grnas, 1:num_cells)

  expect_equal(
    scep_low_with_low_thresh@initial_grna_assignment_list,
    lapply(1:num_grnas, function(i) if (i %in% expressed_grna_ids) 1:num_cells else integer(0)) |>
      setNames(test_data_list$grna_target_data_frame$grna_id)
  )

  # the check is `i == 3` because the target with expressed guides is t1_c3_d1, the 3rd one
  guide_part <- lapply(seq_along(unique_targets), function(i) if (i == 3) 1:num_cells else integer(0)) |>
    setNames(unique_targets)
  nt_part <- lapply(nt_guides, function(nt_guide) integer(0)) |>
    setNames(nt_guides)

  expect_equal(
    scep_low_with_low_thresh@grna_assignments_raw,
    list(
      grna_group_idxs = guide_part,
      indiv_nt_grna_idxs = nt_part,
      all_nt_idxs = integer(0)
    )
  )

  expect_equal(scep_low_with_low_thresh@grnas_per_cell, rep(2, num_cells))

  # no cell is flagged as expressing multiple grnas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scep_low_with_high_thresh <- import_data(
    grna_matrix = grna_matrix_vary_thresh,
    response_matrix = test_data_list$response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "low"
  ) |>
    set_analysis_parameters(
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "thresholding", threshold = 26)

  expect_equal(scep_low_with_high_thresh@cells_w_zero_or_twoplus_grnas, integer(0))

  # `i == 5` is used because g1_t1_c3_d1 is the 5th grna_id and is what all these cells are assigned to
  expect_equal(
    scep_low_with_high_thresh@initial_grna_assignment_list,
    lapply(1:num_grnas, function(i) if (i == 5) 1:num_cells else integer(0)) |>
      setNames(test_data_list$grna_target_data_frame$grna_id)
  )

  # as above, the check is `i == 3` because the target with expressed guides is t1_c3_d1, the 3rd one
  guide_part <- lapply(seq_along(unique_targets), function(i) if (i == 3) 1:num_cells else integer(0)) |>
    setNames(unique_targets)
  nt_part <- lapply(nt_guides, function(nt_guide) integer(0)) |>
    setNames(nt_guides)

  expect_equal(
    scep_low_with_high_thresh@grna_assignments_raw,
    list(
      grna_group_idxs = guide_part,
      indiv_nt_grna_idxs = nt_part,
      all_nt_idxs = integer(0)
    )
  )
  expect_equal(scep_low_with_high_thresh@grnas_per_cell, rep(1, num_cells))
})

test_that("assign_grnas method=threshold moi=high", {
  test_data_list <- make_mock_base_data_for_testing_assign_grnas()
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix_all_0)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  grna_matrix_vary_thresh <- test_data_list$grna_matrix_all_0
  grna_matrix_vary_thresh[1, ] <- 100 # all cells express grna1 very strongly
  grna_matrix_vary_thresh[2, ] <- 50 # all cells also express some of grna2
  grna_matrix_vary_thresh[3:6, 3:6] <- 1:16 # little bit extra to get numerically full rank covariate data frame

  # every cell is flagged as expressing 2 grnas ~~~~~~~~~~~~~~~~~~~~~~~~~
  scep_high_with_low_thresh <- import_data(
    grna_matrix = grna_matrix_vary_thresh,
    response_matrix = test_data_list$response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "thresholding", threshold = 25)

  expect_equal(scep_high_with_low_thresh@grnas_per_cell, rep(2, num_cells))

  expect_equal(
    scep_high_with_low_thresh@initial_grna_assignment_list,
    lapply(1:num_grnas, function(i) if (i <= 2) 1:num_cells else integer(0)) |>
      setNames(test_data_list$grna_target_data_frame$grna_id)
  )

  # no cell is flagged as expressing multiple grnas ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  scep_high_with_high_thresh <- import_data(
    grna_matrix = grna_matrix_vary_thresh,
    response_matrix = test_data_list$response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "thresholding", threshold = 55)

  expect_equal(scep_high_with_high_thresh@grnas_per_cell, rep(1, num_cells))

  expect_equal(
    scep_high_with_high_thresh@initial_grna_assignment_list,
    lapply(1:num_grnas, function(i) if (i == 1) 1:num_cells else integer(0)) |>
      setNames(test_data_list$grna_target_data_frame$grna_id)
  )
})

test_that("assign_grnas method=mixture moi=high", {
  set.seed(1312)
  num_guides <- 3
  num_nt <- 1
  num_cells <- 10
  num_responses <- 20

  grna_target_data_frame <- data.frame(
    grna_id = c(paste0("grna_", 1:num_guides), paste0("nt_", 1:num_nt)),
    grna_target = rep(c("target_1", "non-targeting"), c(num_guides, num_nt)),
    chr = "", start = 0, end = 1
  )

  grna_matrix <- sample(0:2, (num_guides + num_nt) * num_cells, TRUE) |>
    matrix(ncol = num_cells) |>
    `rownames<-`(grna_target_data_frame$grna_id)
  cells_getting_grna_1 <- c(1, 2, 3)
  cells_getting_grna_2 <- c(1, 4, 5)
  cells_getting_grna_3 <- c(1, 6:7)
  cells_getting_nt_1 <- c(2, 8:9)
  # cell 10 has nothing

  expressed_value <- 20
  grna_matrix["grna_1", cells_getting_grna_1] <- expressed_value
  grna_matrix["grna_2", cells_getting_grna_2] <- expressed_value
  grna_matrix["grna_3", cells_getting_grna_3] <- expressed_value
  grna_matrix["nt_1", cells_getting_nt_1] <- expressed_value

  response_matrix <- sample(0:2, num_responses * num_cells, TRUE) |>
    matrix(ncol = num_cells) |>
    `rownames<-`(paste0("response_", 1:num_responses))

  scep_high_mixture <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      discovery_pairs = data.frame(
        grna_target = "target_1", response_id = "response_1"
      )
    ) |>
    assign_grnas(method = "mixture")

  expect_equal(scep_high_mixture@initial_grna_assignment_list$grna_1, cells_getting_grna_1)
  expect_equal(scep_high_mixture@initial_grna_assignment_list$grna_2, cells_getting_grna_2)
  expect_equal(scep_high_mixture@initial_grna_assignment_list$grna_3, cells_getting_grna_3)
  expect_equal(scep_high_mixture@initial_grna_assignment_list$nt_1, cells_getting_nt_1)
})
