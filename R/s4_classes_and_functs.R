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
#' sceptre_object <- create_sceptre_object(
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
#' # 3. assign gRNAs as run cellwise and pairwise QC
#' sceptre_object <- sceptre_object |> assign_grnas()
#' sceptre_object <- sceptre_object |> run_qc()
#'
#' # 4. run the calibration check
#' sceptre_object <- run_calibration_check(sceptre_object)
#' plot(sceptre_object)
#'
#' # 5. run discovery analysis
#' sceptre_object <- run_discovery_analysis(sceptre_object)
#' plot(sceptre_object)
#'
#' # 6. obtain the results for downstream analysis
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
#' sceptre_object <- create_sceptre_object(
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
#' # 4. assign gRNAs as run cellwise and pairwise QC
#' sceptre_object <- sceptre_object |> assign_grnas()
#' sceptre_object <- sceptre_object |> run_qc()
#'
#' # 5. run the calibration check; plot the result
#' sceptre_object <- run_calibration_check(sceptre_object)
#' plot(sceptre_object)
#'
#' # 6. (optional) run the power check
#' sceptre_object <- run_power_check(sceptre_object)
#'
#' # 7. run discovery analysis
#' sceptre_object <- run_discovery_analysis(sceptre_object)
#' plot(sceptre_object)
#'
#' # 8. obtain the results for downstream analysis
#' discovery_result <- get_result(sceptre_object, "discovery")
create_sceptre_object <- function(response_matrix, grna_matrix, grna_group_data_frame, moi, extra_covariates = NULL) {
  # 1. perform initial check
  check_create_sceptre_object_inputs(response_matrix, grna_matrix,
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
  sceptre_object@last_function_called <- "create_sceptre_object"
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
                                    n_nonzero_trt_thresh = 7L,
                                    n_nonzero_cntrl_thresh = 7L,
                                    B1 = 499L, B2 = 4999L, B3 = 24999L) {
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
    formula_object <- auto_construct_formula_object(sceptre_object@covariate_data_frame,
                                                    sceptre_object@low_moi)
  }
  if (identical(discovery_pairs, "all")) {
    discovery_pairs <- generate_all_pairs(sceptre_object@response_matrix,
                                          sceptre_object@grna_group_data_frame)
  }

  # 2. check inputs
  check_prepare_analysis_inputs(response_matrix = sceptre_object@response_matrix,
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
  sceptre_object@n_nonzero_trt_thresh <- n_nonzero_trt_thresh
  sceptre_object@n_nonzero_cntrl_thresh <- n_nonzero_cntrl_thresh
  sceptre_object@B1 <- B1
  sceptre_object@B2 <- B2
  sceptre_object@B3 <- B3
  sceptre_object@last_function_called <- "set_analysis_parameters"
  sceptre_object <- convert_covariate_df_to_design_matrix(sceptre_object)

  # 6. update cached fields
  if (reset_response_precomps) sceptre_object@response_precomputations <- list()

  # return
  return(sceptre_object)

  # compute (i) the NT M matrix, (ii), n nonzero total vector, (iii) n_nonzero_trt, and (iv) n_nonzero_cntrl vectors
  # sceptre_object <- compute_pairwise_qc_information(sceptre_object)
  # compute the number of discovery pairs and (if applicable) pc pairs passing qc
  # sceptre_object <- compute_qc_metrics(sceptre_object)
  #if (update_cached_objects$discard_negative_control_pairs) {
  #  cat("reseting negative control pairs\n")
  #  sceptre_object@negative_control_pairs <- data.frame()
  #}
  #if (update_cached_objects$discard_grna_assignments) {
  #  cat("reseting grna assignments\n")
  #  sceptre_object <- assign_grnas_to_cells(sceptre_object)
  #}
}


# step 3: assign grnas to cells
assign_grnas <- function(sceptre_object, assignment_method = "default", hyperparameters = "default") {
  # 0. verify that function called in correct order
  check_function_call(sceptre_object, "assign_grnas")

  # 1. handle the default arguments
  if (assignment_method == "default") {
    assignment_method <- if (sceptre_object@low_moi) "maximum" else "thresholding"
  }
  if (hyperparameters == "default") {
    hyperparameters <- if (sceptre_object@low_moi) {
      list(umi_fraction_threshold = 0.8)
    } else {
      list(threshold = 5L)
    }
  }

  # 2. check inputs
  check_assign_grna_inputs(sceptre_object, assignment_method, hyperparameters) |> invisible()

  # 3. reset results
  sceptre_object <- reset_results(sceptre_object)

  # 4. determine whether to update cached fields (perhaps add later)

  # 5. update uncached fields
  sceptre_object@grna_assignment_method <- assignment_method
  sceptre_object@grna_assignment_hyperparameters <- hyperparameters
  sceptre_object@last_function_called <- "assign_grnas"

  # 6. assign the grnas
  sceptre_object <- assign_grnas_to_cells(sceptre_object)

  # return
  return(sceptre_object)
}


# step 4: run qc



# step 5: calibration check

# step 6: positive power check

# step 7: discovery analysis















###########################################
# OLDER STUFF
###########################################
prepare_analysis <- function(sceptre_object,
                             discovery_pairs = data.frame(),
                             positive_control_pairs = data.frame(),
                             formula_object = "default",
                             side = "both",
                             fit_parametric_curve = TRUE,
                             control_group = "default",
                             resampling_mechanism = "default",
                             n_nonzero_trt_thresh = 7L,
                             n_nonzero_cntrl_thresh = 7L,
                             grna_assign_threshold = 5L,
                             B1 = 499L, B2 = 4999L, B3 = 24999L,
                             calibration_group_size = NA_integer_,
                             n_calibration_pairs = NA_integer_) {
  # 1. handle default arguments
  updated_args <- handle_default_arguments(sceptre_object, control_group, resampling_mechanism,
                                           formula_object, discovery_pairs)
  discovery_pairs <- updated_args$discovery_pairs
  control_group <- updated_args$control_group
  resampling_mechanism <- updated_args$resampling_mechanism
  formula_object <- updated_args$formula_object

  # 2. check inputs
  check_prepare_analysis_inputs(response_matrix = sceptre_object@response_matrix,
                                grna_matrix = sceptre_object@grna_matrix,
                                covariate_data_frame = sceptre_object@covariate_data_frame,
                                grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                formula_object = formula_object,
                                response_grna_group_pairs_list = list(discovery_pairs, positive_control_pairs),
                                control_group = control_group,
                                resampling_mechanism = resampling_mechanism,
                                side = side, low_moi = sceptre_object@low_moi) |> invisible()

  # 3. determine whether to update cached objects
  update_cached_objects <- check_whether_to_update_cached_objects(sceptre_object, formula_object, grna_assign_threshold,
                                                                  n_calibration_pairs, calibration_group_size, n_nonzero_trt_thresh,
                                                                  n_nonzero_cntrl_thresh, discovery_pairs, control_group)
  # 4. update non-cached fields of the sceptre object
  sceptre_object <- update_uncached_fields(sceptre_object, control_group, resampling_mechanism,
                                           side, fit_parametric_curve, B1, B2, B3, n_nonzero_trt_thresh,
                                           n_nonzero_cntrl_thresh, discovery_pairs, positive_control_pairs,
                                           grna_assign_threshold, formula_object, calibration_group_size,
                                           n_calibration_pairs)

  # update cached fields (three fields cached: response precomputations, grna assignments, negative control pairs)
  if (update_cached_objects$discard_response_precomputations) {
    cat("reseting response precomputations\n")
    sceptre_object@response_precomputations <- list()
  }
  if (update_cached_objects$discard_negative_control_pairs) {
    cat("reseting negative control pairs\n")
    sceptre_object@negative_control_pairs <- data.frame()
  }
  if (update_cached_objects$discard_grna_assignments) {
    cat("reseting grna assignments\n")
    sceptre_object <- assign_grnas_to_cells(sceptre_object)
  }

  # compute (i) the NT M matrix, (ii), n nonzero total vector, (iii) n_nonzero_trt, and (iv) n_nonzero_cntrl vectors
  sceptre_object <- compute_pairwise_qc_information(sceptre_object)

  # compute the number of discovery pairs and (if applicable) pc pairs passing qc
  sceptre_object <- compute_qc_metrics(sceptre_object)

  # create the covariate matrix, check it for correctness, and insert it into object
  sceptre_object <- convert_covariate_df_to_design_matrix(sceptre_object)

  # return
  return(sceptre_object)
}


run_calibration_check <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE) {
  # 1. do argument check
  check_calibration_check_inputs(analysis_prepared = sceptre_object@analysis_prepared,
                                 grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                 control_group_complement = sceptre_object@control_group_complement) |> invisible()

  # 2. get grna assignments
  grna_assignments <- sceptre_object@grna_assignments

  # 3. construct the negative control pairs
  if (is.na(sceptre_object@calibration_group_size)) {
    calibration_group_size <- compute_calibration_group_size(sceptre_object@grna_group_data_frame)
  }  else {
    calibration_group_size <- sceptre_object@calibration_group_size
  }

  cat("Constructing negative control pairs.")
  if (nrow(sceptre_object@negative_control_pairs) == 0L) {
    response_grna_group_pairs <- construct_negative_control_pairs(n_calibration_pairs = sceptre_object@n_calibration_pairs,
                                                                  calibration_group_size = calibration_group_size,
                                                                  grna_assignments = grna_assignments,
                                                                  response_matrix = sceptre_object@response_matrix,
                                                                  grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                                                  low_moi = sceptre_object@low_moi,
                                                                  n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                                                  n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                                                  n_nonzero_m = sceptre_object@M_matrix,
                                                                  n_nonzero_tot = sceptre_object@n_nonzero_tot_vector)
  } else {
    response_grna_group_pairs <- sceptre_object@negative_control_pairs
  }
  cat(crayon::green(' \u2713\n'))

  # 4. generate the set of synthetic indicator idxs (if running permutations)
  if (sceptre_object@run_permutations) {
    cat("Generating permutation resamples.")
    n_cells <- nrow(sceptre_object@covariate_matrix)
    synthetic_idxs <- get_synthetic_permutation_idxs(grna_assignments = grna_assignments,
                                                     B = sceptre_object@B1 + sceptre_object@B2 + sceptre_object@B3,
                                                     calibration_check = TRUE,
                                                     control_group_complement = sceptre_object@control_group_complement,
                                                     calibration_group_size = calibration_group_size,
                                                     n_cells = n_cells)
    cat(crayon::green(' \u2713\n'))
  }
  gc() |> invisible()

  # 5. run the method
  if (sceptre_object@run_permutations) {
    out <- run_perm_test_in_memory(response_matrix = sceptre_object@response_matrix,
                                   grna_assignments = grna_assignments,
                                   covariate_matrix = sceptre_object@covariate_matrix,
                                   response_grna_group_pairs = response_grna_group_pairs,
                                   synthetic_idxs = synthetic_idxs,
                                   output_amount = output_amount,
                                   fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                   B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                   B3 = sceptre_object@B3, calibration_check = TRUE,
                                   control_group_complement = sceptre_object@control_group_complement,
                                   n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                   n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                   side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                   response_precomputations = sceptre_object@response_precomputations,
                                   print_progress = print_progress, parallel = parallel)
  } else {
    stop("CRT not yet implemented.")
    ret <- run_crt_in_memory_v2(response_matrix = sceptre_object@response_matrix,
                                grna_assignments = grna_assignments,
                                covariate_matrix = sceptre_object@covariate_matrix,
                                response_grna_group_pairs = response_grna_group_pairs,
                                output_amount = output_amount,
                                fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                B3 = sceptre_object@B3, calibration_check = TRUE,
                                control_group_complement = sceptre_object@control_group_complement,
                                n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                print_progress = print_progress)
  }

  # update fields of sceptre object
  sceptre_object@calibration_result <- out$ret_pass_qc
  sceptre_object@negative_control_pairs <- response_grna_group_pairs
  sceptre_object@response_precomputations <- out$response_precomputations
  sceptre_object@calibration_check_run <- TRUE
  sceptre_object@last_function_called <- "run_calibration_check"
  return(sceptre_object)
}


run_power_check <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE) {
  # 1. do argument check
  check_discovery_analysis_inputs(response_grna_group_pairs = sceptre_object@positive_control_pairs_with_info,
                                  control_group_complement = sceptre_object@control_group_complement,
                                  grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                  calibration_check_run = sceptre_object@calibration_check_run,
                                  pc_analysis = TRUE) |> invisible()

  # 2. get grna assignments
  grna_assignments <- sceptre_object@grna_assignments

  # 4. generate the set of synthetic indicator idxs (if running permutations)
  if (sceptre_object@run_permutations) {
    cat("Generating permutation resamples.")
    n_cells <- nrow(sceptre_object@covariate_matrix)
    synthetic_idxs <- get_synthetic_permutation_idxs(grna_assignments = grna_assignments,
                                                     B = sceptre_object@B1 + sceptre_object@B2 + sceptre_object@B3,
                                                     calibration_check = FALSE,
                                                     control_group_complement = sceptre_object@control_group_complement,
                                                     calibration_group_size = NULL,
                                                     n_cells = n_cells)
    cat(crayon::green(' \u2713\n'))
  }
  gc() |> invisible()

  # 5. run the method
  if (sceptre_object@run_permutations) {
    out <- run_perm_test_in_memory(response_matrix = sceptre_object@response_matrix,
                                   grna_assignments = grna_assignments,
                                   covariate_matrix = sceptre_object@covariate_matrix,
                                   response_grna_group_pairs = sceptre_object@positive_control_pairs_with_info |> dplyr::filter(pass_qc),
                                   synthetic_idxs = synthetic_idxs,
                                   output_amount = output_amount,
                                   fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                   B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                   B3 = sceptre_object@B3, calibration_check = FALSE,
                                   control_group_complement = sceptre_object@control_group_complement,
                                   n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                   n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                   side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                   response_precomputations = sceptre_object@response_precomputations,
                                   print_progress = print_progress, parallel = parallel)
  } else {
    stop("CRT not yet implemented.")
    ret <- run_crt_in_memory_v2(response_matrix = sceptre_object@response_matrix,
                                grna_assignments = grna_assignments,
                                covariate_matrix = sceptre_object@covariate_matrix,
                                response_grna_group_pairs = sceptre_object@positive_control_pairs_with_info,
                                output_amount = output_amount,
                                fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                B3 = sceptre_object@B3, calibration_check = FALSE,
                                control_group_complement = sceptre_object@control_group_complement,
                                n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                print_progress = print_progress)
  }
  # 6. create the final results data frame by combining the pairs that pass qc with those that do not
  final_result <- data.table::rbindlist(list(out$ret_pass_qc,
                                             sceptre_object@positive_control_pairs_with_info |> dplyr::filter(!pass_qc)), fill = TRUE)
  data.table::setorderv(final_result, cols = c("p_value", "response_id"), na.last = TRUE)

  # update fields of sceptre object
  sceptre_object@power_result <- out$ret
  sceptre_object@response_precomputations <- out$response_precomputations
  sceptre_object@power_check_run <- TRUE
  sceptre_object@last_function_called <- "run_power_check"
  return(sceptre_object)
}


run_discovery_analysis <- function(sceptre_object, output_amount = 1, print_progress = TRUE, parallel = FALSE) {
  # 1. do argument check
  check_discovery_analysis_inputs(response_grna_group_pairs = sceptre_object@discovery_pairs_with_info,
                                  control_group_complement = sceptre_object@control_group_complement,
                                  grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                  calibration_check_run = sceptre_object@calibration_check_run,
                                  pc_analysis = FALSE) |> invisible()

  # 2. get grna assignments
  grna_assignments <- sceptre_object@grna_assignments

  # 4. generate the set of synthetic indicator idxs (if running permutations)
  if (sceptre_object@run_permutations) {
    cat("Generating permutation resamples.")
    n_cells <- nrow(sceptre_object@covariate_matrix)
    synthetic_idxs <- get_synthetic_permutation_idxs(grna_assignments = grna_assignments,
                                                     B = sceptre_object@B1 + sceptre_object@B2 + sceptre_object@B3,
                                                     calibration_check = FALSE,
                                                     control_group_complement = sceptre_object@control_group_complement,
                                                     calibration_group_size = NULL,
                                                     n_cells = n_cells)
    cat(crayon::green(' \u2713\n'))
  }
  gc() |> invisible()

  # 5. run the method on the pairs passing qc
  if (sceptre_object@run_permutations) {
    out <- run_perm_test_in_memory(response_matrix = sceptre_object@response_matrix,
                                   grna_assignments = grna_assignments,
                                   covariate_matrix = sceptre_object@covariate_matrix,
                                   response_grna_group_pairs = sceptre_object@discovery_pairs_with_info |> dplyr::filter(pass_qc),
                                   synthetic_idxs = synthetic_idxs,
                                   output_amount = output_amount,
                                   fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                   B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                   B3 = sceptre_object@B3, calibration_check = FALSE,
                                   control_group_complement = sceptre_object@control_group_complement,
                                   n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                   n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                   side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                   response_precomputations = sceptre_object@response_precomputations,
                                   print_progress = print_progress, parallel = parallel)
  } else {
    stop("CRT not yet implemented.")
    ret <- run_crt_in_memory_v2(response_matrix = sceptre_object@response_matrix,
                                grna_assignments = grna_assignments,
                                covariate_matrix = sceptre_object@covariate_matrix,
                                response_grna_group_pairs = sceptre_object@discovery_pairs_with_info,
                                output_amount = output_amount,
                                fit_parametric_curve = sceptre_object@fit_parametric_curve,
                                B1 = sceptre_object@B1, B2 = sceptre_object@B2,
                                B3 = sceptre_object@B3, calibration_check = FALSE,
                                control_group_complement = sceptre_object@control_group_complement,
                                n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                side_code = sceptre_object@side_code, low_moi = sceptre_object@low_moi,
                                print_progress = print_progress)
  }

  # 6. create the final results data frame by combining the pairs that pass qc with those that do not
  final_result <- data.table::rbindlist(list(out$ret_pass_qc,
                                             sceptre_object@discovery_pairs_with_info |> dplyr::filter(!pass_qc)), fill = TRUE)
  data.table::setorderv(final_result, cols = c("p_value", "response_id"), na.last = TRUE)

  # update fields of sceptre object
  sceptre_object@discovery_result <- final_result
  sceptre_object@response_precomputations <- out$response_precomputations
  sceptre_object@discovery_analysis_run <- TRUE
  sceptre_object@last_function_called <- "run_discovery_analysis"
  return(sceptre_object)
}


get_result <- function(sceptre_object, analysis_type, alpha = 0.1, multiple_testing_correction = "BH") {
  if (!(analysis_type %in% c("calibration", "power", "discovery"))) {
    stop("`analysis_type` must be one of `calibration`, `power`, or `discovery`.")
  }
  if (analysis_type == "calibration" && !sceptre_object@calibration_check_run) {
    stop("Calibration check has not yet been run.")
  }
  if (analysis_type == "power" && !sceptre_object@power_check_run) {
    stop("Power check has not yet been run.")
  }
  if (analysis_type == "discovery" && !sceptre_object@discovery_analysis_run) {
    stop("Discovery analysis not yet run.")
  }
  field_to_extract <- switch(EXPR = analysis_type,
                             calibration = "calibration_result",
                             power = "power_result",
                             discovery = "discovery_result")
  out <- slot(sceptre_object, field_to_extract)
  # if calibration check or discovery analysis, apply multiplicity correction
  if (analysis_type %in% c("calibration", "discovery")) {
    out <- out |>
      dplyr::mutate(reject = stats::p.adjust(p_value, method = multiple_testing_correction) < alpha) |>
      dplyr::select(-pass_qc)
  }
  return(out)
}
