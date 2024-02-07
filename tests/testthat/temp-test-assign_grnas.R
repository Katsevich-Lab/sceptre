
# NOTES
# - anywhere with the tag <--discuss--> is something i want to discuss when presenting this
# - Is `@grnas_per_cell` used for method="maximum"?
# - for method="thresholding" there is some inconsistency with the slots. Is that fine?



# - what happens to  a cell with 0 expression for every response? Deleted in qc? (yeah should be)


# this function hard-codes a particular dataset that will be the base
# of the tests for `assign_grnas()`.
# It is called in every test, so a little repetitive, but it only takes
# ~2 milliseconds per call on my laptop
make_mock_base_data_for_testing_assign_grnas <- function() {
  num_cells <- 24  # must be even for my choice of extra_covariates
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

  # extra_covariates <- data.frame(x = rep(c("b1", "b2"), each = num_cells / 2))

  positive_control_pairs <- data.frame(
    grna_target = on_targets,
    response_id = on_targets
  )

  discovery_pairs <- data.frame(
    grna_target = on_targets,
    response_id = rep(paste0("response_", length(on_targets) + 1), times=2)
  )

  list(
    grna_target_data_frame = grna_target_data_frame,
    response_matrix = response_matrix,
    grna_matrix_all_0 = grna_matrix,
    # extra_covariates = extra_covariates,
    positive_control_pairs = positive_control_pairs,
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
  expect_equal(scep_low_all_1@cells_w_multiple_grnas, 1:num_cells)

  # first element should have all idx
  expect_equal(
    scep_low_all_1@initial_grna_assignment_list,
    lapply(test_data_list$grna_target_data_frame$grna_id, function(target_name)
      if(target_name == test_data_list$grna_target_data_frame$grna_id[1]) 1:num_cells else integer(0)) |>
      setNames(test_data_list$grna_target_data_frame$grna_id)
  )

  # first element should have all idx
  guide_part_all_1 <- lapply(unique_targets, function(target_name)
    if(target_name == test_data_list$grna_target_data_frame$grna_target[1]) 1:num_cells else integer(0)) |>
    setNames(unique_targets)
  nt_part_all_1 <- lapply(nt_guides, function(nt_guide) integer(0)) |>
    setNames(nt_guides)

  expect_equal(
    scep_low_all_1@grna_assignments_raw,
    list(
      grna_group_idxs = guide_part_all_1,
      indiv_nt_grna_idxs = nt_part_all_1
    )
  )

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
  for(i in 1:num_grnas) grna_matrix_clear_max[i,i] <- 100

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
  expected_initial_assignment_list <- lapply(1:num_grnas, function(i) i)  |>
    setNames(test_data_list$grna_target_data_frame$grna_id)
  expect_equal(
    scep_low_clear_max@initial_grna_assignment_list,
    expected_initial_assignment_list
  )
  # all cells from idx `num_grnas+1` onward have no UMI counts at all, so they
  # all are considered to have multiple grnas
  expect_equal(scep_low_clear_max@cells_w_multiple_grnas, (num_grnas+1):num_cells)

  # there are two guides per target, and 3 NT guides
  # this is assembling the expected answer which is a list of lists of indicecs
  guide_part <- lapply(unique_targets, function(target_name) which(test_data_list$grna_target_data_frame$grna_target == target_name)) |>
    setNames(unique_targets)
  nt_part <- lapply(nt_guides, function(nt_guide) which(test_data_list$grna_target_data_frame$grna_id == nt_guide)) |>
    setNames(nt_guides)
  expect_equal(
    scep_low_clear_max@grna_assignments_raw,
    list(
      grna_group_idxs = guide_part,
      indiv_nt_grna_idxs = nt_part
    )
  )

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
  grna_matrix_vary_frac[1,] <- 100 # all cells express grna1 very strongly
  grna_matrix_vary_frac[2,1] <- 25 # cell 1 also expresses some grna2
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
    assign_grnas(method = "maximum", umi_fraction_threshold=.8 - .001)

  expect_equal(scep_low_with_low_frac@cells_w_multiple_grnas, integer(0))

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
    assign_grnas(method = "maximum", umi_fraction_threshold=.8 + .001)

  expect_equal(scep_low_with_high_frac@cells_w_multiple_grnas, 1)
})

# TODO change tests?
test_that("assign_grnas method=threshold moi=low", {

  test_data_list <- make_mock_base_data_for_testing_assign_grnas()
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix_all_0)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  grna_matrix_vary_thresh <- test_data_list$grna_matrix_all_0
  grna_matrix_vary_thresh[1,] <- 100 # all cells express grna1 very strongly
  grna_matrix_vary_thresh[2,] <- 25  # all cells also express some of grna2

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

  expect_equal(scep_low_with_low_thresh@cells_w_multiple_grnas, 1:num_cells)

  # TODO should this fail?? <--discuss-->
  # some grnas are missing and I'm not sure why
  expect_equal(
    scep_low_with_low_thresh@initial_grna_assignment_list,
    lapply(1:num_grnas, function(i) if(i <= 2) 1:num_cells else integer(0)) |>
      setNames(test_data_list$grna_target_data_frame$grna_id)
  )

  guide_part <- lapply(seq_along(unique_targets), function(i) if(i == 1) 1:num_cells else integer(0)) |>
      setNames(unique_targets)
  nt_part <- lapply(nt_guides, function(nt_guide) integer(0)) |>
    setNames(nt_guides)

  # TODO there's a mix of integer(0) and NULL here that seems weird <--discuss-->
  expect_equal(
    scep_low_with_low_thresh@grna_assignments_raw,
    list(
      grna_group_idxs = guide_part,
      indiv_nt_grna_idxs = nt_part
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

  expect_equal(scep_low_with_high_thresh@cells_w_multiple_grnas, integer(0))

  # TODO again -- why are some names missing?? <--discuss-->
  expect_equal(
    scep_low_with_high_thresh@initial_grna_assignment_list,
    lapply(1:num_grnas, function(i) if(i == 1) 1:num_cells else integer(0)) |>
      setNames(test_data_list$grna_target_data_frame$grna_id)
  )

  guide_part <- lapply(seq_along(unique_targets), function(i) if(i == 1) 1:num_cells else integer(0)) |>
    setNames(unique_targets)
  nt_part <- lapply(nt_guides, function(nt_guide) integer(0)) |>
    setNames(nt_guides)

  # TODO there's a mix of integer(0) and NULL here that seems weird <--discuss-->
  expect_equal(
    scep_low_with_high_thresh@grna_assignments_raw,
    list(
      grna_group_idxs = guide_part,
      indiv_nt_grna_idxs = nt_part
    )
  )
  expect_equal(scep_low_with_high_thresh@grnas_per_cell, rep(1, num_cells))
})

# TODO uncomment tests?
test_that("assign_grnas method=threshold moi=high", {

  test_data_list <- make_mock_base_data_for_testing_assign_grnas()
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix_all_0)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  grna_matrix_vary_thresh <- test_data_list$grna_matrix_all_0
  grna_matrix_vary_thresh[1,] <- 100 # all cells express grna1 very strongly
  grna_matrix_vary_thresh[2,] <- 50  # all cells also express some of grna2
  grna_matrix_vary_thresh[3:6,3:6] <- 1:16 # little bit extra to get numerically full rank covariate data frame

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

  # TODO what should this be? <--discuss-->
  # expect_equal(scep_high_with_low_thresh@cells_w_multiple_grnas, 1:num_cells)

  # TODO should this fail? <--discuss-->
  # expect_equal(
  #   scep_high_with_low_thresh@initial_grna_assignment_list,
  #   lapply(1:num_grnas, function(i) if(i <= 2) 1:num_cells else integer(0)) |>
  #     setNames(test_data_list$grna_target_data_frame$grna_id)
  # )

  # guide_part <- lapply(seq_along(unique_targets), function(i) if(i == 1) 1:num_cells else integer(0)) |>
  #   setNames(unique_targets)
  # nt_part <- lapply(nt_guides, function(nt_guide) integer(0)) |>
  #   setNames(nt_guides)
  #
  # # TODO there's a mix of integer(0) and NULL here again <--discuss-->
  # expect_equal(
  #   scep_high_with_low_thresh@grna_assignments_raw,
  #   list(
  #     grna_group_idxs = guide_part,
  #     indiv_nt_grna_idxs = nt_part
  #   )
  # )



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

  expect_equal(scep_high_with_high_thresh@cells_w_multiple_grnas, integer(0))
  expect_equal(scep_high_with_high_thresh@grnas_per_cell, rep(1, num_cells))

  # TODO again -- why are some names missing? <--discuss-->
  # expect_equal(
  #   scep_high_with_high_thresh@initial_grna_assignment_list,
  #   lapply(1:num_grnas, function(i) if(i == 1) 1:num_cells else integer(0)) |>
  #     setNames(test_data_list$grna_target_data_frame$grna_id)
  # )

  # guide_part <- lapply(seq_along(unique_targets), function(i) if(i == 1) 1:num_cells else integer(0)) |>
  #   setNames(unique_targets)
  # nt_part <- lapply(nt_guides, function(nt_guide) integer(0)) |>
  #   setNames(nt_guides)

  # TODO there's a mix of integer(0) and NULL here again <--discuss-->
  # expect_equal(
  #   scep_high_with_high_thresh@grna_assignments_raw,
  #   list(
  #     grna_group_idxs = guide_part,
  #     indiv_nt_grna_idxs = nt_part
  #   )
  # )
})

# TODO write this test
test_that("assign_grnas method=mixture moi=high", {

  test_data_list <- make_mock_base_data_for_testing_assign_grnas()
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix_all_0)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  grna_matrix_mixture <- test_data_list$grna_matrix_all_0
  grna_matrix_mixture[1,] <- 100 # all cells express grna1 very strongly
  grna_matrix_mixture[2,] <- 100  # all cells also express some of grna2
  grna_matrix_mixture[3:6,3:6] <- (1:16)/10 # little bit extra to get numerically full rank covariate data frame

  # every cell is flagged as expressing 2 grnas ~~~~~~~~~~~~~~~~~~~~~~~~~
  scep_high_mixture <- import_data(
    grna_matrix = grna_matrix_mixture,
    response_matrix = test_data_list$response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "mixture")

  # expect_equal(scep_high_mixture@grnas_per_cell, ...)

  # expect_equal(scep_high_mixture@cells_w_multiple_grnas, ...)

  # TODO should this fail? <--discuss-->
  # expect_equal(
  #   scep_high_mixture@initial_grna_assignment_list,
  #   ...
  # )

  # # TODO there's a mix of integer(0) and NULL here again <--discuss-->
  # expect_equal(
  #   scep_high_mixture@grna_assignments_raw,
  #   ...
  # )
})



