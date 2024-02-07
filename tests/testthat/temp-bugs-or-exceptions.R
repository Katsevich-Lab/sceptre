# POSSIBLE BUGS / UNHANDLED EXCEPTIONS?


# 1. problem with thresholding in `assign_grnas()`

# 2. seems really quick to throw low rank error?

bug_threshold_assign_grnas <- function(threshold) {

  test_data_list <- make_mock_base_data_for_testing_run_qc(24)
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

bug_threshold_assign_grnas(threshold=5)
bug_threshold_assign_grnas(threshold=10)


# 2. problem in `run_qc()` if no NTs are considered expressed
# can't replicate

# crash_no_NTs_in_run_qc <- function() {
#   test_data_list <- make_mock_base_data_for_testing_run_qc()
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
#     assign_grnas() |>
#     run_qc()
# }

