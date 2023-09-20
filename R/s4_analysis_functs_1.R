#' Carry out an analysis using `sceptre`
#'
#' @param response_matrix TBD
#' @param grna_matrix TBD
#' @param grna_target_data_frame TBD
#' @param moi TBD
#' @param extra_covariates TBD
#' @param response_names TBD
#'
#' @export
#' @examples
#' #################
#' # Low MOI example
#' #################
#' # 0. obtain the data required for a single-cell screen analysis
#' data(response_matrix_lowmoi) # response-by-cell expression matrix
#' data(grna_matrix_lowmoi) # gRNA-by-cell expression matrix
#' data(extra_covariates_lowmoi) # cell-by-covariate data frame
#' data(grna_target_data_frame_lowmoi) # gRNA group information
#'
#' # 1. create the sceptre object
#' sceptre_object <- import_data(
#' response_matrix = response_matrix_lowmoi,
#' grna_matrix = grna_matrix_lowmoi,
#' extra_covariates = extra_covariates_lowmoi,
#' grna_target_data_frame = grna_target_data_frame_lowmoi,
#' moi = "low")
#' print(sceptre_object)
#'
#' # 2. obtain the discovery and positive control pairs
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_trans_pairs(sceptre_object = sceptre_object,
#' positive_control_pairs = positive_control_pairs)
#'
#' # 3. set the analysis parameters
#' sceptre_object <- set_analysis_parameters(
#' sceptre_object = sceptre_object,
#' discovery_pairs = discovery_pairs,
#' positive_control_pairs = positive_control_pairs)
#' print(sceptre_object)
#'
#' # 4. optional: explicitly assign grnas, run QC
#' plot_grna_count_distributions(sceptre_object)
#' sceptre_object <- sceptre_object |> assign_grnas()
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' sceptre_object <- sceptre_object |> run_qc(p_mito_threshold = 0.1)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 5. run the calibration check
#' sceptre_object <- run_calibration_check(sceptre_object)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 6. run power check
#' sceptre_object <- run_power_check(sceptre_object)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 7. run discovery analysis
#' sceptre_object <- run_discovery_analysis(sceptre_object)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 8. obtain the results for downstream analysis
#' write_outputs_to_directory(sceptre_object = sceptre_object, "~/sceptre_outputs/")
#'
#' ##################
#' # High MOI example
#' ##################
#' # 1. create the sceptre object from CellRanger output
#' directories <- paste0(system.file("extdata", package = "sceptre"), "/highmoi_example/gem_group_", 1:2)
#' data(grna_target_data_frame_highmoi)
#' sceptre_object <- import_data_from_cellranger(directories = directories,
#' moi = "high",
#' grna_target_data_frame = grna_target_data_frame_highmoi)
#' print(sceptre_object)
#'
#' # 2. obtain the response-gRNA group pairs to analyze
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_cis_pairs(sceptre_object,
#' positive_control_pairs = positive_control_pairs,
#' distance_threshold = 5e6)
#'
#' # 3. set the analysis parameters
#' sceptre_object <- set_analysis_parameters(
#' sceptre_object = sceptre_object,
#' discovery_pairs = discovery_pairs,
#' positive_control_pairs = positive_control_pairs,
#' side = "left")
#' print(sceptre_object)
#'
#' # 4 (optional) manually assign grnas, run QC
#' plot_grna_count_distributions(sceptre_object)
#' sceptre_object <- sceptre_object |> assign_grnas()
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' sceptre_object <- sceptre_object |> run_qc(p_mito_threshold = 0.1)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 5. run the calibration check
#' sceptre_object <- run_calibration_check(sceptre_object)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 6. (optional) run the power check
#' sceptre_object <- run_power_check(sceptre_object)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 7. run discovery analysis
#' sceptre_object <- run_discovery_analysis(sceptre_object)
#' plot(sceptre_object)
#' print(sceptre_object)
#'
#' # 8. obtain results; write outputs to directory
#' write_outputs_to_directory(sceptre_object = sceptre_object, "~/sceptre_outputs/")
import_data <- function(response_matrix, grna_matrix, grna_target_data_frame, moi, extra_covariates = NULL, response_names = NA_character_) {
  # 1. perform initial check
  check_import_data_inputs(response_matrix, grna_matrix, grna_target_data_frame, moi, extra_covariates) |> invisible()

  # 2. compute the covariates
  covariate_data_frame <- auto_compute_cell_covariates(response_matrix = response_matrix,
                                                       grna_matrix = grna_matrix,
                                                       extra_covariates = extra_covariates,
                                                       response_names = if (identical(response_names, NA_character_)) rownames(response_matrix) else response_names)

  # 3. make the response matrix row accessible
  response_matrix <- set_matrix_accessibility(response_matrix, make_row_accessible = TRUE)

  # 4. update fields in output object and return
  sceptre_object <- methods::new("sceptre_object")
  sceptre_object@response_matrix <- response_matrix
  sceptre_object@grna_matrix <- grna_matrix
  sceptre_object@covariate_data_frame <- covariate_data_frame
  sceptre_object@grna_target_data_frame <- grna_target_data_frame |> dplyr::mutate(grna_id = as.character(grna_id), grna_target = as.character(grna_target))
  sceptre_object@response_names <- response_names
  sceptre_object@low_moi <- (moi == "low")
  if (!is.null(extra_covariates)) sceptre_object@user_specified_covariates <- colnames(extra_covariates)

  # 5. initialize flags
  sceptre_object@last_function_called <- "import_data"
  sceptre_object@functs_called <- c(import_data = TRUE, set_analysis_parameters = FALSE,
                                    assign_grnas = FALSE, run_qc = FALSE, run_calibration_check = FALSE,
                                    run_power_check = FALSE, run_discovery_analysis = FALSE)
  return(sceptre_object)
}


#' Set analysis parameters
#'
#' @param sceptre_object TBD
#' @param discovery_pairs TBD
#' @param positive_control_pairs TBD
#' @param formula_object TBD
#' @param side TBD
#' @param fit_parametric_curve TBD
#' @param control_group TBD
#' @param resampling_mechanism TBD
#' @param B1 TBD
#' @param B2 TBD
#' @param B3 TBD
#' @param multiple_testing_method TBD
#' @param multiple_testing_alpha TBD
#'
#' @export
set_analysis_parameters <- function(sceptre_object,
                                    discovery_pairs,
                                    positive_control_pairs = data.frame(),
                                    side = "both",
                                    grna_grouping_strategy = "union",
                                    formula_object = "default",
                                    fit_parametric_curve = TRUE,
                                    control_group = "default",
                                    resampling_mechanism = "default",
                                    B1 = 499L, B2 = 4999L, B3 = 24999L,
                                    multiple_testing_method = "BH",
                                    multiple_testing_alpha = 0.1) {
  # 0. verify that function called in correct order
  sceptre_object <- perform_status_check_and_update(sceptre_object, "set_analysis_parameters")

  # 1. handle default arguments
  if (!sceptre_object@low_moi) {
    control_group <- "complement"
    if (resampling_mechanism == "default") resampling_mechanism <- "crt"
  }
  if (sceptre_object@low_moi) {
    if (control_group == "default") control_group <- "nt_cells"
    if (resampling_mechanism == "default") resampling_mechanism <- "permutations"
  }
  if (identical(formula_object, "default")) {
    formula_object <- auto_construct_formula_object(cell_covariates = sceptre_object@covariate_data_frame,
                                                    include_grna_covariates = !sceptre_object@low_moi)
  }

  # 2. check inputs
  check_set_analysis_parameters(sceptre_object = sceptre_object, formula_object = formula_object,
                                response_grna_target_pairs_list = list(discovery_pairs = discovery_pairs,
                                                                       positive_control_pairs = positive_control_pairs),
                                control_group = control_group,
                                resampling_mechanism = resampling_mechanism,
                                side = side, low_moi = sceptre_object@low_moi,
                                grna_grouping_strategy = grna_grouping_strategy) |> invisible()

  # 3. determine whether to reset response precomputations
  reset_response_precomps <- !((length(sceptre_object@formula_object) >= 2) &&
                                 identical(sceptre_object@formula_object[[2L]], formula_object[[2L]]))

  # 4. update uncached fields of the sceptre object
  side_code <- which(side == c("left", "both", "right")) - 2L
  control_group_complement <- control_group == "complement"
  run_permutations <- resampling_mechanism == "permutations"
  sceptre_object@discovery_pairs <- discovery_pairs |> dplyr::mutate(grna_target = as.character(grna_target), response_id = as.character(response_id))
  sceptre_object@positive_control_pairs <- positive_control_pairs |> dplyr::mutate(grna_target = as.character(grna_target), response_id = as.character(response_id))
  sceptre_object@formula_object <- formula_object
  sceptre_object@side_code <- side_code
  sceptre_object@fit_parametric_curve <- fit_parametric_curve
  sceptre_object@control_group_complement <- control_group_complement
  sceptre_object@run_permutations <- run_permutations
  sceptre_object@B1 <- B1
  sceptre_object@B2 <- B2
  sceptre_object@B3 <- B3
  sceptre_object@multiple_testing_alpha <- multiple_testing_alpha
  sceptre_object@multiple_testing_method <- multiple_testing_method
  sceptre_object@grna_grouping_strategy <- grna_grouping_strategy
  sceptre_object@covariate_matrix  <- convert_covariate_df_to_design_matrix(covariate_data_frame = sceptre_object@covariate_data_frame,
                                                                            formula_object = formula_object)

  # 5. modify the grna target df and response grna target dfs
  sceptre_object <- update_dfs_based_on_grouping_strategy(sceptre_object)

  # 6. update cached fields
  if (reset_response_precomps) sceptre_object@response_precomputations <- list()

  # return
  return(sceptre_object)
}


#' Assign gRNAs to cells
#'
#' @param sceptre_object TBD
#' @param method TBD
#' @param hyperparameters TBD
#' @param print_progress TBD
#' @param parallel TBD
#'
#' @export
assign_grnas <- function(sceptre_object, method = "default", hyperparameters = "default", print_progress = TRUE, parallel = FALSE) {
  # 0. verify that function called in correct order
  sceptre_object <- perform_status_check_and_update(sceptre_object, "assign_grnas")

  # 1. handle the default arguments
  if (identical(method, "default")) {
    method <- if (sceptre_object@low_moi) "maximum" else "mixture"
  }
  hyperparameters_default <- if (method == "maximum") {
    list(umi_fraction_threshold = 0.8)
  } else if (method == "thresholding")  {
    list(threshold = 5)
  } else if (method == "mixture") {
    list(n_em_rep = 5L, pi_guess_range = c(1e-5, 0.1),
      g_pert_guess_range = log(c(10, 5000)), n_nonzero_cells_cutoff = 10L,
      backup_threshold = 5, probability_threshold = 0.8,
      formula_object = auto_construct_formula_object(cell_covariates = sceptre_object@covariate_data_frame,
                                                     include_grna_covariates = TRUE))
  } else if (method == "user_supplied") {
    list()
  }
  if (identical(hyperparameters, "default")) hyperparameters <- hyperparameters_default
  if (methods::is(hyperparameters, "list")) {
    hyperparam_names <- names(hyperparameters)
    for (hyperparam_name in hyperparam_names) hyperparameters_default[[hyperparam_name]] <- hyperparameters[[hyperparam_name]]
    hyperparameters <- hyperparameters_default
  }

  # 2. check inputs
  check_assign_grna_inputs(sceptre_object, method, hyperparameters) |> invisible()

  # 3. determine whether to reset response precomputations
  reset_response_precomps <- sceptre_object@low_moi &&
    (!identical(sceptre_object@grna_assignment_method, method) ||
     !identical(sceptre_object@grna_assignment_hyperparameters, hyperparameters))
  if (reset_response_precomps) sceptre_object@response_precomputations <- list()

  # 4. update uncached fields
  sceptre_object@grna_assignment_method <- method
  sceptre_object@grna_assignment_hyperparameters <- hyperparameters

  # 5. assign the grnas
  sceptre_object <- assign_grnas_to_cells(sceptre_object, print_progress, parallel)

  # return
  return(sceptre_object)
}


#' Run QC
#'
#' @param sceptre_object TBD
#' @param n_nonzero_trt_thresh TBD
#' @param n_nonzero_cntrl_thresh TBD
#' @param response_n_umis_range TBD
#' @param response_n_nonzero_range TBD
#' @param p_mito_threshold TBD
#' @param additional_cells_to_remove TBD
#'
#' @export
run_qc <- function(sceptre_object,
                   n_nonzero_trt_thresh = 7L,
                   n_nonzero_cntrl_thresh = 7L,
                   response_n_umis_range = c(0.01, 0.99),
                   response_n_nonzero_range = c(0.01, 0.99),
                   p_mito_threshold = 0.2,
                   additional_cells_to_remove = integer()) {
  # 1. verify that function called in correct order
  sceptre_object <- perform_status_check_and_update(sceptre_object, "run_qc")

  # 2. check inputs
  check_run_qc_inputs(n_nonzero_trt_thresh,
                      n_nonzero_cntrl_thresh,
                      response_n_umis_range,
                      response_n_nonzero_range) |> invisible()

  # 3. obtain previous cells_in_use for caching purposes
  current_cells_in_use <- sceptre_object@cells_in_use

  # 4. update uncached fields of the sceptre object
  sceptre_object@n_nonzero_trt_thresh <- n_nonzero_trt_thresh
  sceptre_object@n_nonzero_cntrl_thresh <- n_nonzero_cntrl_thresh

  # 5. determine the cells to retain after cellwise qc
  sceptre_object <- determine_cells_to_retain(sceptre_object, response_n_umis_range, response_n_nonzero_range,
                                              p_mito_threshold, additional_cells_to_remove)

  # 6. determine whether to reset response precomputation
  if (!identical(current_cells_in_use, sceptre_object@cells_in_use)) {
    sceptre_object@response_precomputations <- list()
  }

  # 7. update the grna assignments given the cellwise qc
  sceptre_object <- update_grna_assignments_given_qc(sceptre_object)

  # 8. compute (i) the NT M matrix, (ii), n nonzero total vector, (iii) n_nonzero_trt, and (iv) n_nonzero_cntrl vectors
  sceptre_object <- compute_pairwise_qc_information(sceptre_object)

  # 9. compute the number of discovery pairs and (if applicable) pc pairs passing qc
  sceptre_object <- compute_qc_metrics(sceptre_object)

  # return
  return(sceptre_object)
}
