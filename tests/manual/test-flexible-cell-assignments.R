# the goal of this script is to ensure that the new version of sceptre, when called
# with default arguments, produces the same results as the existing version of sceptre

# this script needs to be run after restarting the R session, to restore the
# existing version of sceptre prior to running the new version

run_lowmoi_example <- function(){
  library(sceptre)
#  library(sceptredata)
  # 1. create the sceptre object
  data("lowmoi_example_data")
  sceptre_object <- import_data(
    response_matrix = lowmoi_example_data$response_matrix,
    grna_matrix = lowmoi_example_data$grna_matrix,
    extra_covariates = lowmoi_example_data$extra_covariates,
    grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
    moi = "low"
  )

  # 2. set the analysis parameters
  positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
  discovery_pairs <- construct_trans_pairs(
    sceptre_object = sceptre_object,
    positive_control_pairs = positive_control_pairs,
    pairs_to_exclude = "pc_pairs"
  ) # |>
#    dplyr::filter(response_id == "PCBP3")

  sceptre_object <- set_analysis_parameters(
    sceptre_object = sceptre_object,
    discovery_pairs = discovery_pairs,
    positive_control_pairs = positive_control_pairs,
    grna_integration_strategy = "union"
  )

  # 3. assign grnas
  sceptre_object <- sceptre_object |>
    assign_grnas(method = "thresholding")

  # 4. run qc
  sceptre_object <- sceptre_object |>
    run_qc(p_mito_threshold = 0.075)

  # 5. run the calibration check
  sceptre_object <- run_calibration_check(sceptre_object, parallel = TRUE, n_processors = 2)

  # 6. run power check
  sceptre_object <- run_power_check(sceptre_object, parallel = TRUE, n_processors = 2)

  # 7. run discovery analysis
  sceptre_object <- run_discovery_analysis(sceptre_object, parallel = TRUE, n_processors = 2)

  return(sceptre_object)
}

run_highmoi_example <- function(){
  library(sceptre)
#  library(sceptredata)

  # 1. create the sceptre object from cellranger output
  directories <- paste0(
    system.file("extdata", package = "sceptre"),
    "/highmoi_example/gem_group_", c(1, 2)
  )
  data(grna_target_data_frame_highmoi)
  sceptre_object <- import_data_from_cellranger(
    directories = directories,
    moi = "high",
    grna_target_data_frame = grna_target_data_frame_highmoi
  )

  # 2. set the analysis parameters
  positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
  discovery_pairs <- construct_cis_pairs(sceptre_object,
                                         positive_control_pairs = positive_control_pairs,
                                         distance_threshold = 5e6
  ) # |>
#    dplyr::slice_head(n = 100)

  sceptre_object <- set_analysis_parameters(
    sceptre_object = sceptre_object,
    discovery_pairs = discovery_pairs,
    positive_control_pairs = positive_control_pairs,
    side = "left"
  )

  # 3. assign grnas
  sceptre_object <- sceptre_object |> assign_grnas(parallel = TRUE, n_processors = 2)

  # 4. run qc
  plot_covariates(sceptre_object, p_mito_threshold = 0.075)
  sceptre_object <- sceptre_object |> run_qc(p_mito_threshold = 0.075)

  # 5. run the calibration check
  sceptre_object <- run_calibration_check(sceptre_object, parallel = TRUE, n_processors = 2)

  # 6. run the power check
  sceptre_object <- run_power_check(sceptre_object, parallel = TRUE, n_processors = 2)

  # 7. run discovery analysis
  sceptre_object <- run_discovery_analysis(sceptre_object, parallel = TRUE, n_processors = 2)

  return(sceptre_object)
}

# run on existing package
sceptre_object_old_lowmoi <- run_lowmoi_example()
sceptre_object_old_highmoi <- run_highmoi_example()

# load new package
devtools::load_all()

# run on new package
sceptre_object_new_lowmoi <- run_lowmoi_example()
sceptre_object_new_highmoi <- run_highmoi_example()

# check that the low-MOI results match
expect_equal(sceptre_object_new_lowmoi@discovery_result |> dplyr::select(-n_trt, -n_cntrl),
             sceptre_object_old_lowmoi@discovery_result)

expect_equal(sceptre_object_new_lowmoi@calibration_result |> dplyr::select(-n_trt, -n_cntrl),
             sceptre_object_old_lowmoi@calibration_result)

expect_equal(sceptre_object_new_lowmoi@power_result |> dplyr::select(-n_trt, -n_cntrl),
             sceptre_object_old_lowmoi@power_result)

# check that the high-MOI results match
expect_equal(sceptre_object_new_highmoi@discovery_result |> dplyr::select(-n_trt, -n_cntrl),
             sceptre_object_old_highmoi@discovery_result)

expect_equal(sceptre_object_new_highmoi@calibration_result |> dplyr::select(-n_trt, -n_cntrl),
             sceptre_object_old_highmoi@calibration_result)

expect_equal(sceptre_object_new_highmoi@power_result |> dplyr::select(-n_trt, -n_cntrl),
             sceptre_object_old_highmoi@power_result)
