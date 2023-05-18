# the "master" run sceptre function


#' Run sceptre
#'
#' @param response_matrix (required) a matrix of raw expression counts. The responses (e.g., genes or proteins) should be in the rows, and the cells should be in the columns. The row names should be the unique IDs of the responses. The matrix can be a standard (dense) matrix or a sparse matrix of class \code{dgCMatrix}, \code{dgRMatrix}, or \code{dgTMatrix}.
#' @param grna_matrix (required) a matrix of gRNA expression counts. The gRNAs should be in the rows, and the cells should be in the columns. The row names should be the unique IDs of the gRNAs. The matrix can be a standard (dense) matrix or a sparse matrix of class \code{dgCMatrix}, \code{dgRMatrix}, or \code{dgTMatrix}. (See note about passing logical gRNA matrices.)
#' @param covariate_data_frame (required) a data frame containing the cell-specific covariates (e.g., sequencing batch, total UMI count, etc.), with cells in the rows and covariates in the columns. The column names should be the names of the covariates.
#' @param grna_group_data_frame (required) a data frame specifying the group to which each individual gRNA belongs. This data frame should contain columns \code{grna_id} and \code{grna_group}, with \code{grna_group} specifying the group membership of a given gRNA in the \code{grna_id} column. Typically, gRNAs that target the same site are grouped into the same gRNA group. IMPORTANT: Non-targeting gRNAs should be assigned a gRNA group label of "non-targeting".
#' @param formula_object (required) an R formula object specifying how \code{sceptre} is to adjust for the covariates (in the \code{covariate_data_frame}).
#' @param response_grna_group_pairs (required) a data frame specifying the set of response-gRNA group pairs to test for association. The data frame should contain columns \code{response_id} and \code{grna_group}.
#' @param calibration_check (required) a logical (i.e., \code{TRUE}/\code{FALSE}) value indicating whether to run a calibration check (\code{TRUE}) or a discovery analysis (\code{FALSE}). See the "Details" section for a more complete description of the calibration and discovery analyses.
#' @param moi the MOI of the experiment, either "low" or "high"
#' @param control_group the set of cells against which the cells that received a targeting peturbation are compared, either "complement" (for the complement set) or "nt" (for the NT cells)
#' @param n_nonzero_trt_thresh (optional; moderate importance; default 7) For a given response-gRNA group pair, \code{n_nonzero_trt} is the number of "treatment cells" (i.e., cells that have received a targeting gRNA) with nonzero expression. \code{sceptre} filters for response-gRNA group pairs with an \code{n_nonzero_trt} value equal to or greater than \code{n_nonzero_trt_thresh}.
#' @param n_nonzero_cntrl_thresh (optional; moderate importance; default 7) \code{n_nonzero_cntrl} is the number of "control cells" (i.e., cells that have received an NT gRNA) with nonzero expression. \code{sceptre} filters for pairs with a \code{n_nonzero_control} value equal to or greater than \code{n_nonzero_control_thresh}.
#' @param return_debugging_metrics (optional; lower importance; default \code{FALSE}) return metrics that facilitate debugging calibration issues? (See "Note" for details.)
#' @param return_resampling_dist (optional; lower importance; default \code{FALSE}) return the resampling distribution of the null test statistics?
#' @param calibration_group_size (optional; lower importance; default \code{NULL}) the number of NT gRNAs to assign to each negative control gRNA group in the calibration check analysis. By default, \code{sceptre} sets \code{calibration_group_size} to the median group size of the targeting gRNAs.
#' @param n_calibration_pairs (optional; lower importance; default \code{NULL}) the number of negative control pairs to analyze in the calibration check. By default, the number of negative control pairs that \code{sceptre} analyzes is equal to the number of discovery pairs (i.e., pairs specified in the data frame \code{response_grna_group_pairs}) that passes pairwise QC.
#' @param fit_skew_normal (optional; lower importance; default \code{TRUE}) fit a skew-normal distribution to the null distribution of test statistics and use this fitted distribution to compute a more accurate p-value?
#' @param B1 (optional; lower importance; default 499) the number of null test statistics to compute in the first round of the permutation test.
#' @param B2 (optional; lower importance; default 4999) the number of null test statistics to compute the second round of the permutation test.
#' @param B3 (optional; lower importance; default 24999) the number of null test statistics to compute the third round of the permutation test.
#' @param regression_method (optional; lower importance; default "poisson_glm") The regression method to use to regress out the covariates, either "poisson_glm" for a Poisson GLM or "nb_glm" for a negative binomial GLM.
#' @param print_progress (optional; default TRUE) print the progress of the function?
#'
#' @return a data frame containing the following columns: \code{response_id}, \code{grna_group}, \code{n_nonzero_trt}, \code{n_nonzero_cntrl}, and \code{p_value}.
#' \itemize{
#' \item{\code{response_id}}: the response ID
#' \item{\code{grna_group}}: the gRNA group. When running a calibration check, \code{grna_group} is a negative control gRNA group constructed by grouping together individual NT gRNAs. When running a discovery analysis, by contrast, \code{grna_group} is a targeting gRNA group.
#' \item{\code{n_nonzero_trt}}: the number of treatment cells (i.e., cells that have received a targeting gRNA) with nonzero expression
#' \item{\code{n_nonzero_cntrl}}: the number of control cells (i.e., cells that have received a non-targeting gRNA) with nonzero expression
#' \item{\code{p_value}}: the \code{sceptre} p-value
#' \item{\code{log_2_fold_change}}: the \code{sceptre}-estimated log (base 2) fold change in expression
#' }
#' Rows are ordered according to \code{p_value}. Response-gRNA group pairs that fail to pass pairwise QC are assigned a p-value and log fold change estimate of \code{NA}.
#'
#' @export
#'
#' @examples
#' library(Matrix)
#'
#' ####################
#' # EXAMPLE 1: LOW MOI
#' ####################
#'
#' # 0. load the data associated with the experiment
#' data(response_matrix_lowmoi) # response-by-cell expression matrix
#' data(grna_matrix_lowmoi) # gRNA-by-cell expression matrix
#' data(covariate_data_frame_lowmoi) # cell-by-covariate data frame
#' data(grna_group_data_frame_lowmoi) # gRNA group information
#'
#' # 1. obtain the set of pairs to analyze
#' response_grna_group_pairs <- generate_all_pairs(response_matrix_lowmoi,
#' grna_group_data_frame_lowmoi)
#'
#' # 2. set the formula object
#' formula_object <- formula(~log(response_n_umis) + log(response_n_nonzero) +
#' bio_rep + p_mito)
#'
#' # 3. run the calibration check analysis (NOTE: calibration_check set to TRUE)
#' calibration_result <- run_sceptre(response_matrix = response_matrix_lowmoi,
#' grna_matrix = grna_matrix_lowmoi,
#' covariate_data_frame = covariate_data_frame_lowmoi,
#' grna_group_data_frame = grna_group_data_frame_lowmoi,
#' formula_object = formula_object,
#' response_grna_group_pairs = response_grna_group_pairs,
#' calibration_check = TRUE,
#' moi = "low",
#' control_group = "complement")
run_sceptre <- function(response_matrix, grna_matrix,
                        covariate_data_frame, grna_group_data_frame,
                        formula_object, response_grna_group_pairs,
                        calibration_check, moi, control_group,
                        n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L,
                        grna_assign_threshold = 5L, return_debugging_metrics = FALSE,
                        return_resampling_dist = FALSE, fit_skew_normal = TRUE,
                        calibration_group_size = NULL, n_calibration_pairs = NULL,
                        B1 = 499L, B2 = 4999L, B3 = 24999L,
                        regression_method = "poisson_glm", print_progress = TRUE) {
  ###############
  # PART 1: SETUP
  ###############
  cat("Running setup. ")
  # 1. check function input arguments
  check_inputs(response_matrix, grna_matrix, covariate_data_frame,
               grna_group_data_frame, formula_object, calibration_check,
               response_grna_group_pairs, regression_method, moi, control_group) |> invisible()

  # 2. order the pairs to analyze data frame
  response_grna_group_pairs <- order_pairs_to_analyze(response_grna_group_pairs)

  # 3. harmonize arguments (called for side-effects)
  harmonize_arguments(return_resampling_dist, fit_skew_normal, moi, control_group) |> invisible()

  # 4. make the response matrix row accessible
  response_matrix <- set_matrix_accessibility(response_matrix, make_row_accessible = TRUE)

  # 5. convert the cell covariate data frame into a design matrix
  covariate_matrix <- convert_covariate_df_to_design_matrix(covariate_data_frame, formula_object)
  rm(covariate_data_frame)

  # 6. assign gRNAs to cells
  grna_assignments <- assign_grnas_to_cells(grna_matrix, grna_group_data_frame, threshold, low_moi, control_group_complement, calibration_check, n_calibration_pairs)

  # 7. construct the negative control pairs
  # 8. generate the set of synthetic indicator idxs

  ####################
  # PART 2: RUN METHOD
  ####################

}
