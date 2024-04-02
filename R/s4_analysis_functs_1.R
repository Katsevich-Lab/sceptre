#' Set analysis parameters
#'
#' `set_analysis_parameters()` sets the analysis parameters that control how the statistical analysis is to be conducted. See \href{https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html}{Chapter 2 of the manual} for more detailed information about this function.
#'
#' @note Every argument to this function is optional, but typically, users want to specify `discovery_pairs` at minimum.
#'
#' @param sceptre_object a `sceptre_object`
#' @param discovery_pairs (optional) a data frame with columns `grna_target` and `response_id` specifying the discovery pairs to analyze
#' @param positive_control_pairs (optional) a data frame with columns `grna_target` and `response_id` specifying the positive control pairs to analyze
#' @param formula_object (optional) a formula object specifying how to adjust for the covariates in the model
#' @param side (optional; default `"both"`) the sidedness of the test, one of `"left"`, `"right"`, or `"both"`
#' @param grna_integration_strategy (optional; default `"union"`) a string specifying the gRNA integration strategy, either `"singleton"`, `"union"`, or `"bonferroni"`
#' @param resampling_approximation (optional; default `"skew_normal"`) a string indicating the resampling approximation to make to the null distribution of test statistics, either `"skew_normal"` or `"no_approximation"`
#' @param control_group (optional) a string specifying the control group to use in the differential expression analysis, either `"complement"` or `"nt_cells"`
#' @param resampling_mechanism (optional) a string specifying the resampling mechanism to use, either `"permutations"` or `"crt"`
#' @param multiple_testing_method (optional; default `"BH"`) a string specifying the multiple testing correction method to use; see `p.adjust.methods` for options
#' @param multiple_testing_alpha (optional; default `0.1`) a numeric specifying the nominal level of the multiple testing correction method
#'
#' @return an updated `sceptre_object` in which the analysis parameters have been set
#'
#' @export
#' @examples
#' library(sceptredata)
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#'
#' # set analysis parameters
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_cis_pairs(sceptre_object,
#'   positive_control_pairs = positive_control_pairs,
#'   distance_threshold = 5e6
#' )
#' sceptre_object <- sceptre_object |>
#'   set_analysis_parameters(
#'     discovery_pairs = discovery_pairs,
#'     positive_control_pairs = positive_control_pairs,
#'     side = "left"
#'   )
set_analysis_parameters <- function(sceptre_object,
                                    discovery_pairs = data.frame(grna_target = character(0),
                                                                 response_id = character(0)),
                                    positive_control_pairs = data.frame(grna_target = character(0),
                                                                        response_id = character(0)),
                                    side = "both",
                                    grna_integration_strategy = "union",
                                    formula_object = "default",
                                    resampling_approximation = "skew_normal",
                                    control_group = "default",
                                    resampling_mechanism = "default",
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
    formula_object <- auto_construct_formula_object(
      cell_covariates = sceptre_object@covariate_data_frame,
      include_grna_covariates = !sceptre_object@low_moi
    )
  }
  B1 <- 499L
  if (resampling_approximation == "skew_normal") {
    B2 <- 4999L
    B3 <- if (resampling_mechanism == "permutations") 24999L else 0L
  } else if (resampling_approximation == "no_approximation") {
    B2 <- 0L # no curve fitting; thus, B2 = 0L
    B3 <- 0L # to be updated in the run_qc step
  }

  # 2. check inputs
  check_set_analysis_parameters(
    sceptre_object = sceptre_object,
    formula_object = formula_object,
    response_grna_target_pairs_list = list(
      discovery_pairs = discovery_pairs,
      positive_control_pairs = positive_control_pairs
    ),
    control_group = control_group,
    resampling_mechanism = resampling_mechanism,
    side = side, low_moi = sceptre_object@low_moi,
    grna_integration_strategy = grna_integration_strategy,
    resampling_approximation = resampling_approximation
  ) |> invisible()

  # 3. determine whether to reset response precomputations
  reset_response_precomps <- !((length(sceptre_object@formula_object) >= 2) &&
    identical(sceptre_object@formula_object[[2L]], formula_object[[2L]]))

  # 4. update uncached fields of the sceptre object
  side_code <- which(side == c("left", "both", "right")) - 2L
  control_group_complement <- control_group == "complement"
  run_permutations <- resampling_mechanism == "permutations"
  sceptre_object@discovery_pairs <- discovery_pairs |> dplyr::mutate(grna_target = as.character(grna_target), response_id = as.character(response_id))
  sceptre_object@positive_control_pairs <- positive_control_pairs |> dplyr::mutate(grna_target = as.character(grna_target), response_id = as.character(response_id))
  sceptre_object@n_discovery_pairs <- nrow(sceptre_object@discovery_pairs)
  sceptre_object@n_positive_control_pairs <- nrow(sceptre_object@positive_control_pairs)
  sceptre_object@formula_object <- formula_object
  sceptre_object@side_code <- side_code
  sceptre_object@resampling_approximation <- resampling_approximation
  sceptre_object@control_group_complement <- control_group_complement
  sceptre_object@run_permutations <- run_permutations
  sceptre_object@B1 <- B1
  sceptre_object@B2 <- B2
  sceptre_object@B3 <- B3
  sceptre_object@multiple_testing_alpha <- multiple_testing_alpha
  sceptre_object@multiple_testing_method <- multiple_testing_method
  sceptre_object@grna_integration_strategy <- grna_integration_strategy
  sceptre_object@covariate_matrix <- convert_covariate_df_to_design_matrix(
    covariate_data_frame = sceptre_object@covariate_data_frame,
    formula_object = formula_object
  )

  # 5. modify the grna target df and response grna target dfs
  sceptre_object <- update_dfs_based_on_grouping_strategy(sceptre_object)

  # 6. update cached fields
  if (reset_response_precomps) sceptre_object@response_precomputations <- list()

  # return
  return(sceptre_object)
}


#' Assign gRNAs to cells
#'
#' `assign_grnas()` performs the gRNA-to-cell assignments. `sceptre` provides three gRNA-to-cell assignment strategies: the mixture method, the thresholding method, and the maximum method. The mixture method involves assigning gRNAs to cells using a principled mixture model. Next, the thresholding method assigns a gRNA to a cell if the UMI count of the gRNA in the cell is greater than or equal to some integer threshold. Finally, the maximum method assigns the gRNA that accounts for the greatest number of UMIs in a given cell to that cell. The maximum method is available only in low MOI. See \href{https://timothy-barry.github.io/sceptre-book/assign-grnas.html}{Chapter 3 of the manual} for more detailed information about `assign_grnas()`.
#'
#' @note See the manual for information about the method-specific additional arguments.
#' @param sceptre_object a `sceptre_object`
#' @param method (optional) a string indicating the method to use to assign the gRNAs to cells, one of `"mixture"`, `"thresholding"`, or `"maximum"`
#' @param print_progress (optional; default `TRUE`) a logical indicating whether to print progress updates
#' @param parallel (optional; default `FALSE`) a logical indicating whether to run the function in parallel
#' @param n_processors (optional; default "auto") an integer specifying the number of processors to use if `parallel` is set to `TRUE`. The default, `"auto"`, automatically detects the number of processors available on the machine.
#' @param log_dir (optional; default `tempdir()`) a string indicating the directory in which to write the log files (ignored if `parallel = FALSE`)
#' @param ... optional method-specific additional arguments
#'
#' @return an updated `sceptre_object` in which the gRNA assignments have been carried out
#' @export
#' @examples
#' library(sceptredata)
#' data("lowmoi_example_data")
#' # 1. import data, set default analysis parameters
#' sceptre_object <- import_data(
#'   response_matrix = lowmoi_example_data$response_matrix,
#'   grna_matrix = lowmoi_example_data$grna_matrix,
#'   extra_covariates = lowmoi_example_data$extra_covariates,
#'   grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
#'   moi = "low"
#' ) |> set_analysis_parameters()
#'
#' # 2. assign gRNAs (three different methods)
#' sceptre_object <- sceptre_object |> assign_grnas(method = "thresholding")
#' sceptre_object <- sceptre_object |> assign_grnas(method = "maximum")
#' sceptre_object <- sceptre_object |> assign_grnas(
#'   method = "mixture", parallel = TRUE, n_processors = 2
#' )
assign_grnas <- function(sceptre_object, method = "default", print_progress = TRUE, parallel = FALSE,
                         n_processors = "auto", log_dir = tempdir(), ...) {
  # 0. verify that function called in correct order
  sceptre_object <- perform_status_check_and_update(sceptre_object, "assign_grnas")

  # 1. handle the default arguments
  if (identical(method, "default")) {
    method <- if (sceptre_object@low_moi) "maximum" else "mixture"
  }
  hyperparameters_default <- if (method == "maximum") {
    list(
      umi_fraction_threshold = 0.8,
      min_grna_n_umis_threshold = 5L
    )
  } else if (method == "thresholding") {
    list(threshold = 5)
  } else if (method == "mixture") {
    list(
      n_em_rep = 5L, pi_guess_range = c(1e-5, 0.1),
      g_pert_guess_range = log(c(10, 5000)), n_nonzero_cells_cutoff = 10L,
      backup_threshold = 5L, probability_threshold = 0.8,
      formula_object = auto_construct_formula_object(
        cell_covariates = sceptre_object@covariate_data_frame,
        include_grna_covariates = TRUE
      )
    )
  }
  hyperparameters <- list(...)
  if (length(hyperparameters) == 0L) hyperparameters <- hyperparameters_default
  for (hyperparam_name in names(hyperparameters)) hyperparameters_default[[hyperparam_name]] <- hyperparameters[[hyperparam_name]]
  hyperparameters <- hyperparameters_default

  # 2. check inputs
  check_assign_grna_inputs(sceptre_object, method, hyperparameters, n_processors) |> invisible()

  # 3. determine whether to reset response precomputations
  reset_response_precomps <- sceptre_object@low_moi &&
    (!identical(sceptre_object@grna_assignment_method, method) ||
      !identical(sceptre_object@grna_assignment_hyperparameters, hyperparameters))
  if (reset_response_precomps) sceptre_object@response_precomputations <- list()

  # 4. update uncached fields
  sceptre_object@grna_assignment_method <- method
  sceptre_object@grna_assignment_hyperparameters <- hyperparameters

  # 5. assign the grnas
  sceptre_object <- assign_grnas_to_cells(sceptre_object, print_progress, parallel, n_processors, log_dir)

  # return
  return(sceptre_object)
}


#' Run QC
#'
#' `run_qc()` runs cellwise and pairwise QC on the data. Cellwise QC involves filtering cells on the covariates `response_n_nonzero`, `response_n_umis`, and `response_p_mito`. In low-MOI we additionally remove cells that contain zero or multiple gRNAs. Next, pairwise QC involves filtering out target-response pairs whose data are too sparse to be analyzed reliably. In this context we define the “number of nonzero treatment cells” (resp., the “number of nonzero control cells”) as the number of cells in the treatment group (resp., control group) that contain nonzero expression of the response. (We sometimes use the shorthand `n_nonzero_trt` and `n_nonzero_cntrl` to refer to the number of nonzero treatment cells and control cells, respectively.) Pairwise QC involves filtering target-response pairs on `n_nonzero_trt` and `n_nonzero_cntrl`. See \href{https://timothy-barry.github.io/sceptre-book/run-qc.html}{Chapter 4 of the manual} for more detailed information about this function.
#'
#' @param sceptre_object a `sceptre_object`
#' @param n_nonzero_trt_thresh (optional; default `7L`) an integer specifying the number of nonzero *treatment* cells a pair must contain for it to be retained
#' @param n_nonzero_cntrl_thresh (optional; default `7L`) an integer specifying the number of nonzero *control* cells a pair must contain for it to be retained
#' @param response_n_umis_range (optional; default `c(0.01, 0.99)`) a length-two vector of percentiles specifying the location at which to clip the left and right tails of the `response_n_umis` distribution
#' @param response_n_nonzero_range (optional; default `c(0.01, 0.99)`) a length-two vector of percentiles specifying the location at which to clip the left and right tails of the `response_n_nonzero` distribution
#' @param p_mito_threshold (optional; default `0.2`) a numeric value specifying the location at which to clip the right tail of the `response_p_mito` distribution
#' @param additional_cells_to_remove (optional) a vector of integer indices specifying additional cells to remove
#'
#' @return an updated `sceptre_object` in which cellwise and pairwise QC have been applied
#' @export
#' @examples
#' library(sceptredata)
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#'
#' # set analysis parameters, assign grnas, run qc
#' positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
#' discovery_pairs <- construct_cis_pairs(sceptre_object,
#'   positive_control_pairs = positive_control_pairs,
#'   distance_threshold = 5e6
#' )
#' sceptre_object <- sceptre_object |>
#'   set_analysis_parameters(
#'     discovery_pairs = discovery_pairs,
#'     positive_control_pairs = positive_control_pairs,
#'     side = "left"
#'   ) |>
#'   assign_grnas(method = "thresholding") |>
#'   run_qc()
run_qc <- function(sceptre_object,
                   n_nonzero_trt_thresh = 7L,
                   n_nonzero_cntrl_thresh = 7L,
                   response_n_umis_range = c(0.01, 0.99),
                   response_n_nonzero_range = c(0.01, 0.99),
                   p_mito_threshold = 0.2,
                   additional_cells_to_remove = integer()) {
  run_qc_pt_1(
    sceptre_object,
    n_nonzero_trt_thresh,
    n_nonzero_cntrl_thresh,
    response_n_umis_range,
    response_n_nonzero_range,
    p_mito_threshold,
    additional_cells_to_remove
  ) |>
    run_qc_pt_2()
}


run_qc_pt_1 <- function(sceptre_object,
                        n_nonzero_trt_thresh = 7L,
                        n_nonzero_cntrl_thresh = 7L,
                        response_n_umis_range = c(0.01, 0.99),
                        response_n_nonzero_range = c(0.01, 0.99),
                        p_mito_threshold = 0.2,
                        additional_cells_to_remove = integer()) {
  # cellwise start
  # 1. verify that function called in correct order
  sceptre_object <- perform_status_check_and_update(sceptre_object, "run_qc")

  # 2. check inputs
  check_run_qc_inputs(
    n_nonzero_trt_thresh,
    n_nonzero_cntrl_thresh,
    response_n_umis_range,
    response_n_nonzero_range,
    sceptre_object@initial_grna_assignment_list
  ) |> invisible()

  # 2.5 update fields of sceptre_object
  sceptre_object@cellwise_qc_thresholds <- list(
    response_n_umis_range = response_n_umis_range,
    response_n_nonzero_range = response_n_nonzero_range,
    p_mito_threshold = p_mito_threshold
  )

  # 3. obtain previous cells_in_use for caching purposes
  current_cells_in_use <- sceptre_object@cells_in_use

  # 4. update uncached fields of the sceptre object
  sceptre_object@n_nonzero_trt_thresh <- as.integer(n_nonzero_trt_thresh)
  sceptre_object@n_nonzero_cntrl_thresh <- as.integer(n_nonzero_cntrl_thresh)

  # 5. determine the cells to retain after cellwise qc
  sceptre_object <- determine_cells_to_retain(
    sceptre_object, response_n_umis_range, response_n_nonzero_range,
    p_mito_threshold, additional_cells_to_remove
  )

  # 6. determine whether to reset response precomputation
  if (!identical(current_cells_in_use, sceptre_object@cells_in_use)) {
    sceptre_object@response_precomputations <- list()
  }

  # 7. update the grna assignments given the cellwise qc
  sceptre_object <- update_grna_assignments_given_qc(sceptre_object)
  return(sceptre_object)
}


run_qc_pt_2 <- function(sceptre_object) {
  # pairwise start
  # 8. compute (i) the NT M matrix, (ii), n nonzero total vector, (iii) n_nonzero_trt, and (iv) n_nonzero_cntrl vectors
  sceptre_object <- compute_pairwise_qc_information(sceptre_object)

  # 9. compute the number of discovery pairs and (if applicable) pc pairs passing qc
  sceptre_object <- compute_qc_metrics(sceptre_object)

  # update B3, the number of resamples to draw, if resampling_approximation is no_approximation
  if (sceptre_object@resampling_approximation == "no_approximation") {
    mult_fact <- if (sceptre_object@side_code == 0L) 10 else 5
    sceptre_object@B3 <- ceiling(mult_fact * max(
      sceptre_object@n_ok_discovery_pairs,
      sceptre_object@n_ok_positive_control_pairs
    ) / sceptre_object@multiple_testing_alpha) |>
      as.integer()
  }

  # return
  return(sceptre_object)
}
