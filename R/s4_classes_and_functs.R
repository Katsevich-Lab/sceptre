#' Carry out an analysis using `sceptre`
#'
#' @param response_matrix (required) a matrix of raw expression counts. The responses (e.g., genes or proteins) should be in the rows, and the cells should be in the columns. The row names should be the unique IDs of the responses. The matrix can be a standard (dense) matrix or a sparse matrix of class \code{dgCMatrix}, \code{dgRMatrix}, or \code{dgTMatrix}.
#' @param grna_matrix (required) a matrix of gRNA expression counts. The gRNAs should be in the rows, and the cells should be in the columns. The row names should be the unique IDs of the gRNAs. The matrix can be a standard (dense) matrix or a sparse matrix of class \code{dgCMatrix}, \code{dgRMatrix}, or \code{dgTMatrix}. (See "Notes" for details about passing a gRNA matrix containing user-specified gRNA-to-cell assignments.)
#' @param covariate_data_frame (required) a data frame containing the cell-specific covariates (e.g., sequencing batch, total UMI count, etc.), with cells in the rows and covariates in the columns.
#' @param grna_group_data_frame (required) a data frame specifying the "group" to which each individual gRNA belongs. This data frame should contain columns \code{grna_id} and \code{grna_group}, with \code{grna_group} specifying the group membership of a given gRNA in the \code{grna_id} column. Typically, gRNAs that target the same site are grouped into the same gRNA group. Non-targeting gRNAs (if present) should be assigned a gRNA group label of "non-targeting".
#' @param moi the MOI of the dataset, either "low" or "high".
#' @param formula_object (required) an R formula object specifying how to adjust for the covariates in the \code{covariate_data_frame}.
#' @param response_grna_group_pairs (required) a data frame specifying the set of response-gRNA group pairs to test for association. The data frame should contain columns \code{response_id} and \code{grna_group}.
#' @param calibration_check (required) a logical (i.e., \code{TRUE}/\code{FALSE}) value indicating whether to run a calibration check (\code{TRUE}) or a discovery analysis (\code{FALSE}). See the "Details" section for a more complete description of the calibration and discovery analyses.
#' @param side (optional; moderate importance; default "both") the sidedness of the test, either two-sided ("both"), left-tailed ("left"), or right-tailed ("right").
#' @param control_group (optional; moderate importance; default depends on MOI) the set of cells against which the cells that received a given targeting gRNA group are compared, either "complement" (for the complement set) or "nt" (for the cells containing a non-targeting gRNA). The only valid option for high MOI data is "complement", as few (if any) cells contain exclusively non-targeting gRNAs in high MOI. The default for low MOI data is "nt".
#' @param resampling_mechanism (optional; moderate importance; default depends on MOI) the resampling mechanism used to carry out inference, either "crt" (for the conditional randomization test) or "permutations" (for the permutation test). The default is to use "crt" in high MOI and "permutations" in low MOI.
#' @param n_nonzero_trt_thresh (optional; moderate importance; default 7) For a given response-gRNA group pair, \code{n_nonzero_trt} is the number of "treatment cells" (i.e., cells that have received the given targeting gRNA group) with nonzero expression. \code{sceptre} filters for response-gRNA group pairs with an \code{n_nonzero_trt} value equal to or greater than \code{n_nonzero_trt_thresh}.
#' @param n_nonzero_cntrl_thresh (optional; moderate importance; default 7) For a given response-gRNA group pair, \code{n_nonzero_cntrl} is the number of "control cells" (i.e., cells against which the treatment cells are compared) with nonzero expression. \code{sceptre} filters for pairs with a \code{n_nonzero_control} value equal to or greater than \code{n_nonzero_control_thresh}.
#' @param grna_assign_threshold (optional; lower importance; default 5) the threshold used to assign gRNAs to cells on high MOI data. A given gRNA is assigned to a cell if the UMI count of the gRNA within that cell is greater than or equal to the threshold. (In low MOI, this argument is ignored, as gRNAs are assigned to cells via a maximum as opposed to thresholding operation.)
#' @param calibration_group_size (optional; lower importance; default \code{NULL}) the number of NT gRNAs to assign to each negative control gRNA group in the calibration check. By default, \code{calibration_group_size} is set equal to the median group size of the targeting gRNA groups.
#' @param n_calibration_pairs (optional; lower importance; default \code{NULL}) the number of negative control pairs to analyze in the calibration check. By default, this argument is set equal to the number of discovery pairs (i.e., pairs specified in the data frame \code{response_grna_group_pairs}) that passes pairwise QC.
#' @param fit_parametric_curve (optional; moderate importance; default \code{TRUE}) a logical (i.e., \code{TRUE}/\code{FALSE}) indicating whether to fit a parametric curve to the distribution of null test statistics and to compute a p-value using this fitted parametric curve. Setting this argument to \code{TRUE} (the default) enables \code{sceptre} to return very small p-values.
#' @param B1 (optional; lower importance; default 499) the number of null test statistics to compute in the first stage of the permutation test.
#' @param B2 (optional; lower importance; default 4999) the number of null test statistics to compute the second stage of the permutation test.
#' @param B3 (optional; lower importance; default 24999) the number of null test statistics to compute the third stage of the permutation test.
#' @param print_progress (optional; lower importance; default TRUE) a logical (i.e., \code{TRUE}/\code{FALSE} value) indicating whether to print the progress of the function.
#' @param output_amount (optional; lower importance; default 1) an integer specifying the amount of information to return as part of the output, either 1 (least), 2 (intermediate), or 3 (most). (See "Value".)
#'
#' @return an object of class `sceptre_object`
#' @export
#'
#' @examples
#' #################
#' # Low MOI example
#' #################
#' # 0. obtain the data required for a single-cell screen analysis
#' data(response_matrix_lowmoi) # response-by-cell expression matrix
#' data(grna_matrix_lowmoi) # gRNA-by-cell expression matrix
#' data(extra_covariates_lowmoi) # cell-by-covariate data frame
#' data(grna_group_data_frame_lowmoi) # gRNA group information
#'
#' # 1. create the sceptre object
#' sceptre_object <- import_data(
#' response_matrix = response_matrix_lowmoi,
#' grna_matrix = grna_matrix_lowmoi,
#' extra_covariates = extra_covariates_lowmoi,
#' grna_group_data_frame = grna_group_data_frame_lowmoi,
#' moi = "low")
#'
#' # 2. prepare the analysis
#' sceptre_object <- set_analysis_parameters(
#' sceptre_object = sceptre_object,
#' discovery_pairs = "all")
#'
#' # 3. run the calibration check
#' sceptre_object <- run_calibration_check(sceptre_object)
#' plot(sceptre_object)
#'
#' # 4. run discovery analysis
#' sceptre_object <- run_discovery_analysis(sceptre_object)
#' plot(sceptre_object)
#'
#' # 5. obtain the results for downstream analysis
#' discovery_result <- get_result(sceptre_object, "discovery")
#'
#' ##################
#' # High MOI example
#' ##################
#' # 0. obtain the data required for a single-cell screen analysis
#' data(response_matrix_highmoi_experimental)
#' data(grna_matrix_highmoi_experimental)
#' data(extra_covariates_highmoi_experimental)
#' data(grna_group_data_frame_highmoi_experimental)
#'
#' # 1. create the sceptre object
#' sceptre_object <- import_data(
#' response_matrix = response_matrix_highmoi_experimental,
#' grna_matrix = grna_matrix_highmoi_experimental,
#' grna_group_data_frame = grna_group_data_frame_highmoi_experimental,
#' moi = "high",
#' extra_covariates = extra_covariates_highmoi_experimental)
#'
#' # 2. set the discovery pairs and positive control pairs
#' data(discovery_pairs_highmoi_experimental)
#' data(pc_pairs_highmoi_experimental)
#'
#' # 3. prepare the analysis
#' sceptre_object <- set_analysis_parameters(
#' sceptre_object = sceptre_object,
#' discovery_pairs = discovery_pairs_highmoi_experimental,
#' positive_control_pairs = pc_pairs_highmoi_experimental,
#' resampling_mechanism = "permutations",
#' side = "left")
#'
#' # 4. run the calibration check; plot the result
#' sceptre_object <- run_calibration_check(sceptre_object)
#' plot(sceptre_object)
#'
#' # 5. (optional) run the power check
#' sceptre_object <- run_power_check(sceptre_object)
#'
#' # 6. run discovery analysis
#' sceptre_object <- run_discovery_analysis(sceptre_object)
#' plot(sceptre_object)
#'
#' # 7. obtain the results for downstream analysis
#' discovery_result <- get_result(sceptre_object, "discovery")
import_data <- function(response_matrix, grna_matrix, grna_group_data_frame, moi, extra_covariates = NULL) {
  # 1. perform initial check
  check_import_data_inputs(response_matrix, grna_matrix,
                                     grna_group_data_frame, moi, extra_covariates) |> invisible()

  # 2. compute the covariates
  covariate_data_frame <- auto_compute_cell_covariates(response_matrix, grna_matrix, extra_covariates)

  # 3. make the response matrix row accessible
  response_matrix <- set_matrix_accessibility(response_matrix, make_row_accessible = TRUE)

  # 4. update fields in output object and return
  sceptre_object <- new("sceptre_object")
  sceptre_object@response_matrix <- response_matrix
  sceptre_object@grna_matrix <- grna_matrix
  sceptre_object@covariate_data_frame <- covariate_data_frame
  sceptre_object@grna_group_data_frame <- grna_group_data_frame
  sceptre_object@low_moi <- (moi == "low")
  if (!is.null(extra_covariates)) sceptre_object@user_specified_covariates <- colnames(extra_covariates)
  sceptre_object@last_function_called <- "import_data"
  return(sceptre_object)
}


# step 2: set analysis parameters
set_analysis_parameters <- function(sceptre_object,
                                    discovery_pairs,
                                    positive_control_pairs = data.frame(),
                                    formula_object = "default",
                                    side = "both",
                                    fit_parametric_curve = TRUE,
                                    control_group = "default",
                                    resampling_mechanism = "default",
                                    B1 = 499L, B2 = 4999L, B3 = 24999L,
                                    multiple_testing_method = "BH",
                                    multiple_testing_alpha = 0.1) {
  # 0. verify that function called in correct order
  check_function_call(sceptre_object, "set_analysis_parameters")

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
  if (identical(discovery_pairs, "all")) {
    discovery_pairs <- generate_all_pairs(sceptre_object@response_matrix,
                                          sceptre_object@grna_group_data_frame)
  }

  # 2. check inputs
  check_set_analysis_parameters(response_matrix = sceptre_object@response_matrix,
                                grna_matrix = sceptre_object@grna_matrix,
                                covariate_data_frame = sceptre_object@covariate_data_frame,
                                grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                formula_object = formula_object,
                                response_grna_group_pairs_list = list(discovery_pairs, positive_control_pairs),
                                control_group = control_group,
                                resampling_mechanism = resampling_mechanism,
                                side = side, low_moi = sceptre_object@low_moi) |> invisible()

  # 3. reset results
  sceptre_object <- reset_results(sceptre_object)

  # 4. determine whether to reset response precomputations
  reset_response_precomps <- !((length(sceptre_object@formula_object) >= 2) &&
                                 identical(sceptre_object@formula_object[[2L]], formula_object[[2L]]))

  # 5. update uncached fields of the sceptre object
  side_code <- which(side == c("left", "both", "right")) - 2L
  control_group_complement <- control_group == "complement"
  run_permutations <- resampling_mechanism == "permutations"
  sceptre_object@discovery_pairs <- discovery_pairs
  sceptre_object@positive_control_pairs <- positive_control_pairs
  sceptre_object@formula_object <- formula_object
  sceptre_object@side_code <- side_code
  sceptre_object@fit_parametric_curve <- fit_parametric_curve
  sceptre_object@control_group_complement <- control_group_complement
  sceptre_object@run_permutations <- run_permutations
  sceptre_object@B1 <- B1
  sceptre_object@B2 <- B2
  sceptre_object@B3 <- B3
  sceptre_object@last_function_called <- "set_analysis_parameters"
  sceptre_object@multiple_testing_alpha <- multiple_testing_alpha
  sceptre_object@multiple_testing_method <- multiple_testing_method
  sceptre_object@covariate_matrix  <- convert_covariate_df_to_design_matrix(covariate_data_frame = sceptre_object@covariate_data_frame,
                                                                            formula_object = formula_object)

  # 6. update cached fields
  if (reset_response_precomps) sceptre_object@response_precomputations <- list()

  # return
  return(sceptre_object)
}


# step 3: assign grnas to cells
assign_grnas <- function(sceptre_object, assignment_method = "default", hyperparameters = "default") {
  # 0. verify that function called in correct order
  check_function_call(sceptre_object, "assign_grnas")

  # 1. handle the default arguments
  if (identical(assignment_method, "default")) {
    assignment_method <- if (sceptre_object@low_moi) "maximum" else "thresholding"
  }
  if (identical(hyperparameters, "default")) {
    hyperparameters <- if (assignment_method == "maximum") {
      list(umi_fraction_threshold = 0.5)
    } else if (assignment_method == "thresholding")  {
      list(threshold = 5L)
    } else if (assignment_method == "mixture") {
      list()
    } else {
      list()
    }
  }

  # 2. check inputs
  check_assign_grna_inputs(sceptre_object, assignment_method, hyperparameters) |> invisible()

  # 3. reset results
  sceptre_object <- reset_results(sceptre_object)

  # 4. determine whether to reset response precomputations
  reset_response_precomps <- !(identical(sceptre_object@grna_assignment_method, assignment_method) &&
                               identical(sceptre_object@grna_assignment_hyperparameters, hyperparameters))
  if (reset_response_precomps) sceptre_object@response_precomputations <- list()

  # 5. update uncached fields
  sceptre_object@grna_assignment_method <- assignment_method
  sceptre_object@grna_assignment_hyperparameters <- hyperparameters
  sceptre_object@last_function_called <- "assign_grnas"

  # 6. assign the grnas
  sceptre_object <- assign_grnas_to_cells(sceptre_object)

  # return
  return(sceptre_object)
}

# step 4: run cellwise and pairwise qc
run_qc <- function(sceptre_object,
                   n_nonzero_trt_thresh = 7L,
                   n_nonzero_cntrl_thresh = 7L,
                   response_n_umis_range = c(0.01, 0.99),
                   additional_cells_to_remove = integer()) {
  # 0. verify that function called in correct order
  check_function_call(sceptre_object, "run_qc")

  # 2. check inputs
  check_run_qc_inputs(n_nonzero_trt_thresh,
                      n_nonzero_cntrl_thresh,
                      response_n_umis_range) |> invisible()

  # 3. reset results
  sceptre_object <- reset_results(sceptre_object)

  # 4. obtain previous cells_in_use for caching purposes
  current_cells_in_use <- sceptre_object@cells_in_use

  # 5. update uncached fields of the sceptre object
  sceptre_object@n_nonzero_trt_thresh <- n_nonzero_trt_thresh
  sceptre_object@n_nonzero_cntrl_thresh <- n_nonzero_cntrl_thresh
  sceptre_object@last_function_called <- "run_qc"

  # 6. determine the cells to retain after cellwise qc
  sceptre_object <- determine_cells_to_retain(sceptre_object, response_n_umis_range, additional_cells_to_remove)

  # 7. determine whether to reset response precomputation
  if (!identical(current_cells_in_use, sceptre_object@cells_in_use)) {
    sceptre_object@response_precomputations <- list()
  }

  # 8. update the grna assignments given the cellwise qc
  sceptre_object <- update_grna_assignments_given_qc(sceptre_object)

  # 9. compute (i) the NT M matrix, (ii), n nonzero total vector, (iii) n_nonzero_trt, and (iv) n_nonzero_cntrl vectors
  sceptre_object <- compute_pairwise_qc_information(sceptre_object)

  # 10. compute the number of discovery pairs and (if applicable) pc pairs passing qc
  sceptre_object <- compute_qc_metrics(sceptre_object)

  # return
  return(sceptre_object)
}
