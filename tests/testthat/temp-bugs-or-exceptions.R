# POSSIBLE BUGS / UNHANDLED EXCEPTIONS?


# 1. problem with thresholding in `assign_grnas()`
# 2. problem with `assign_grnas()` that I don't understand


bug_threshold_assign_grnas <- function(threshold) {

  test_data_list <- make_mock_base_data_for_testing_run_qc(24) # from temp-test-run_qc.R in tests/testthat
  num_cells <- ncol(test_data_list$response_matrix)
  num_grnas <- nrow(test_data_list$grna_matrix)
  num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
  unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
  nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])

  set.seed(123)
  grna_matrix <- test_data_list$grna_matrix
  grna_matrix[] <- rpois(prod(dim(grna_matrix)), 1.5)
  for(i in 1:num_grnas) grna_matrix[i,i] <- 10 # for clear expression

  scep_high <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = test_data_list$response_matrix,
    grna_target_data_frame = test_data_list$grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      positive_control_pairs = test_data_list$positive_control_pairs,
      discovery_pairs = test_data_list$discovery_pairs
    ) |>
    assign_grnas(method = "thresholding", threshold = threshold) # works when taking to 5

}

## run these
# bug_threshold_assign_grnas(threshold=5) # works
# bug_threshold_assign_grnas(threshold=6) # fails



bug_error_in_assign_grnas_involving_matrices <- function() {
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
  grna_matrix["id3", cells_expressing_t3] <- 50
  grna_matrix["nt1", cells_expressing_nt1] <- 50

  response_matrix <- matrix(sample(0:1, num_responses * num_cells, replace=TRUE, prob = c(.9, .1)), num_responses, num_cells) |>
    `rownames<-`(c("t1", "t2", "t3", paste0("response_", 4:num_responses)))

  # target t1: I will have all non-zero for both the treatment and control cells
  # whether or not the counts reflect this will depend on if there are "control" cells
  # other than `cells_expressing_nt1`
  response_matrix["t1", cells_expressing_t1] <- 100
  response_matrix["t1", cells_expressing_nt1] <- 100

  response_matrix["t2", cells_expressing_t2] <- 100
  response_matrix["t2", cells_expressing_nt1] <- 0

  response_matrix["t3", cells_expressing_t3] <- 0
  response_matrix["t3", cells_expressing_nt1] <- 100

  positive_control_pairs = data.frame(
    grna_target = c("t1", "t2", "t3"),
    response_id = c("t1", "t2", "t3")
  )
  discovery_pairs <- data.frame(
    grna_target = c("t1",        "t1",        "t2",         "t2"),
    response_id = c("response_4", "response_5",  "response_4", "response_6")
  )

  scep_high_all_pass <- import_data(
    grna_matrix = grna_matrix,
    response_matrix = response_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "high"
  ) |>
    set_analysis_parameters(
      positive_control_pairs = positive_control_pairs,
      discovery_pairs = discovery_pairs,
      control_group = "complement"
    ) |>
    assign_grnas(method = "thresholding", threshold = 40)
}
## run this
# bug_weird_error()  # throws error




# 3. problem in `run_qc()` if no NTs are considered expressed
# can't replicate

# crash_no_NTs_in_run_qc <- function() {
#   test_data_list <- make_mock_base_data_for_testing_run_qc(num_cells=24)
#   num_cells <- ncol(test_data_list$response_matrix)
#   num_grnas <- nrow(test_data_list$grna_matrix)
#   num_targets <- sum(test_data_list$grna_target_data_frame$grna_target != "non-targeting")
#   unique_targets <- unique(test_data_list$grna_target_data_frame$grna_target[1:num_targets])
#   nt_guides <- with(test_data_list$grna_target_data_frame, grna_id[grna_target == "non-targeting"])
#
#   set.seed(123)
#   grna_matrix <- test_data_list$grna_matrix
#   grna_matrix[] <- rpois(prod(dim(grna_matrix)), .5)
#   grna_matrix[1,] <- 10
#   # grna_matrix[9:11,] <- 0 # no expression for NTs
#
#   scep_high <- import_data(
#     grna_matrix = grna_matrix,
#     response_matrix = test_data_list$response_matrix,
#     grna_target_data_frame = test_data_list$grna_target_data_frame,
#     moi = "high"
#   ) |>
#     set_analysis_parameters(
#       positive_control_pairs = test_data_list$positive_control_pairs,
#       discovery_pairs = test_data_list$discovery_pairs
#     ) |>
#     assign_grnas(method="thresholding", threshold=10) |>
#     run_qc()
# }
