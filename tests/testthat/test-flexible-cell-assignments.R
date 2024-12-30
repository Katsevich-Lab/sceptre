# # TODO: Test whether matches current package
#
# library(sceptredata)
# # 1. create the sceptre object
# data("lowmoi_example_data")
# sceptre_object <- import_data(
#   response_matrix = lowmoi_example_data$response_matrix,
#   grna_matrix = lowmoi_example_data$grna_matrix,
#   extra_covariates = lowmoi_example_data$extra_covariates,
#   grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#   moi = "low"
# )
#
# # 2. set the analysis parameters
# positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
# discovery_pairs <- construct_trans_pairs(
#   sceptre_object = sceptre_object,
#   positive_control_pairs = positive_control_pairs,
#   pairs_to_exclude = "pc_pairs"
# ) |>
#   dplyr::filter(response_id == "PCBP3")
#
# sceptre_object <- set_analysis_parameters(
#   sceptre_object = sceptre_object,
#   discovery_pairs = discovery_pairs,
#   positive_control_pairs = positive_control_pairs,
#   grna_integration_strategy = "union",
#   treatment_group = "inclusive"
# )
#
# # 3. assign grnas
# sceptre_object <- sceptre_object |> assign_grnas(method = "mixture")
#
# # 4. run qc
# sceptre_object <- sceptre_object |>
#   run_qc(p_mito_threshold = 0.075,
#          remove_cells_w_zero_or_twoplus_grnas = TRUE)
#
# # 5. run the calibration check
# sceptre_object <- run_calibration_check(sceptre_object, parallel = TRUE, n_processors = 2)
#
# # 6. run power check
# sceptre_object <- run_power_check(sceptre_object, parallel = TRUE, n_processors = 2)
#
# # 7. run discovery analysis
# sceptre_object <- run_discovery_analysis(sceptre_object, parallel = TRUE, n_processors = 2)
#
#
