test_that("construct_cis_pairs loads packaged default response positions", {
  data(highmoi_example_data)
  data(grna_target_data_frame_highmoi)

  sceptre_object <- import_data(
    response_matrix = highmoi_example_data$response_matrix,
    grna_matrix = highmoi_example_data$grna_matrix,
    grna_target_data_frame = grna_target_data_frame_highmoi,
    moi = "high",
    extra_covariates = highmoi_example_data$extra_covariates,
    response_names = highmoi_example_data$gene_names
  )
  positive_control_pairs <- construct_positive_control_pairs(sceptre_object)

  expect_no_error({
    cis_pairs <- construct_cis_pairs(
      sceptre_object,
      positive_control_pairs = positive_control_pairs,
      distance_threshold = 5e6
    )
  })
  expect_named(cis_pairs, c("grna_target", "response_id"))
  expect_gt(nrow(cis_pairs), 0)
})
