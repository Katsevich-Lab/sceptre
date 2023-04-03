#######
# SETUP
#######
data(response_matrix_lowmoi)
data(grna_matrix_lowmoi)
data(covariate_data_frame_lowmoi)
data(grna_group_data_frame_lowmoi)
formula_object <- formula(~log(response_n_umis) + log(response_n_nonzero) + bio_rep + p_mito)

grp_size <- sample(seq(2, 4), size = 1)
n_pairs <- sample(x = seq(10, 7000), size = 1)
response_grna_group_pairs <- generate_all_pairs(response_matrix_lowmoi,
                                                grna_group_data_frame_lowmoi) |> dplyr::sample_n(n_pairs)
n_nonzero_trt_thresh <- sample(x = seq(0L, 20), size = 1)
n_nonzero_cntrl_thresh <- sample(x = seq(0L, 50), size = 1)
grna_assignments <- assign_grnas_to_cells_lowmoi_v2(grna_matrix_lowmoi, grna_group_data_frame_lowmoi, TRUE, NULL)

calibration_result <- run_sceptre_lowmoi(response_matrix = response_matrix_lowmoi,
                                         grna_matrix = grna_matrix_lowmoi,
                                         covariate_data_frame = covariate_data_frame_lowmoi,
                                         grna_group_data_frame = grna_group_data_frame_lowmoi,
                                         formula_object = formula_object,
                                         response_grna_group_pairs = response_grna_group_pairs,
                                         calibration_check = TRUE,
                                         n_nonzero_trt_thresh = n_nonzero_trt_thresh,
                                         n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
                                         calibration_group_size = grp_size)

discovery_result <- run_sceptre_lowmoi(response_matrix = response_matrix_lowmoi,
                                       grna_matrix = grna_matrix_lowmoi,
                                       covariate_data_frame = covariate_data_frame_lowmoi,
                                       grna_group_data_frame = grna_group_data_frame_lowmoi,
                                       formula_object = formula_object,
                                       response_grna_group_pairs = response_grna_group_pairs,
                                       calibration_check = FALSE,
                                       n_nonzero_trt_thresh = n_nonzero_trt_thresh,
                                       n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,)
######
# TEST
######
test_that("calibration all pass qc", {
  expect_gte(min(calibration_result$n_nonzero_trt), n_nonzero_trt_thresh)
  expect_gte(min(calibration_result$n_nonzero_cntrl), n_nonzero_cntrl_thresh)
})


test_that("undercover group size", {
  n_nt <- strsplit(x = as.character(calibration_result$grna_group[1]),
                   split = "&", fixed = TRUE)[[1]] |> length()
  expect_equal(n_nt, grp_size)
})


test_that("calibration effective sample size", {
  samp <- sample(x = seq(1, nrow(calibration_result)), size = 10)
  for (idx in samp) {
    row <- calibration_result[idx,]
    response_id <- as.character(row$response_id)
    nt_grnas <- strsplit(x = as.character(row$grna_group), split = "&", fixed = TRUE)[[1]]
    y <- response_matrix_lowmoi[response_id,]
    y_nt <- y[grna_assignments$all_nt_idxs]
    n_trt <- sum(y_nt[unlist(grna_assignments$indiv_nt_grna_idxs[nt_grnas])] >= 1)
    n_cntrl <- sum(y_nt >= 1) - n_trt
    
    expect_equal(n_trt, row$n_nonzero_trt)
    expect_equal(n_cntrl, row$n_nonzero_cntrl)
  }
})


test_that("calibration discovery sample size agreement", {
  expect_equal(nrow(na.omit(discovery_result)), nrow(calibration_result))
})


test_that("discovery all pass qc", {
  discovery_result_pass_qc <- na.omit(discovery_result)
  expect_true(all(discovery_result_pass_qc$n_nonzero_trt >= n_nonzero_trt_thresh))
  expect_true(all(discovery_result_pass_qc$n_nonzero_cntrl >= n_nonzero_cntrl_thresh))  
})


test_that("discovery effective sample size", {
  samp <- sample(x = seq(1, nrow(calibration_result)), size = 10)
  for (idx in samp) {
    row <- discovery_result[idx,]
    response_id <- as.character(row$response_id)
    grna_group <- as.character(row$grna_group)
    y <- response_matrix_lowmoi[response_id,]
    n_cntrl <- sum(y[grna_assignments$all_nt_idxs] >= 1)
    n_trt <- sum(y[grna_assignments$grna_group_idxs[[grna_group]]] >= 1)
    expect_equal(n_trt, row$n_nonzero_trt)
    expect_equal(n_cntrl, row$n_nonzero_cntrl)
  }
})
