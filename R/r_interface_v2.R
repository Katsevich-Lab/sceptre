#' Run \code{sceptre}
#'
#' Applies \code{sceptre} to analyze a single-cell CRISPR screen dataset.
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
#' @return A data frame containing the following columns: \code{response_id}, \code{grna_group}, \code{n_nonzero_trt}, \code{n_nonzero_cntrl}, and \code{p_value}.
#' \itemize{
#' \item{\code{response_id}}: the ID of the response
#' \item{\code{grna_group}}: the gRNA group against which the response is tested. When running a calibration check, \code{grna_group} is a negative control gRNA group constructed by grouping together individual NT gRNAs. When running a discovery analysis, by contrast, \code{grna_group} is a targeting gRNA group.
#' \item{\code{n_nonzero_trt}}: the number of treatment cells (i.e., cells that have received a targeting gRNA) with nonzero expression
#' \item{\code{n_nonzero_cntrl}}: the number of control cells (i.e., cells that have received a non-targeting gRNA) with nonzero expression
#' \item{\code{p_value}}: the \code{sceptre} p-value
#' \item{\code{log_2_fold_change}}: the \code{sceptre}-estimated log (base 2) fold change in expression, controlling for the covariates
#' }
#' Rows are ordered according to \code{p_value}. Response-gRNA group pairs that fail to pass pairwise QC are assigned a p-value and log fold change estimate of \code{NA}. When \code{output_amount} is set to 2, the following additional columns are returned:
#'
#' \itemize{
#'  \item{\code{stage}}: the stage at which the p-value was computed, one of 1, 2, or 3 (See Notes)
#'  \item{\code{z_orig}}: the original z-score computed on the raw data
#'  \item{\code{xi}, \code{omega}, \code{alpha}}: the fitted parameters of the skew-normal distribution. (These columns are set to \code{NA} if a parametric curve was not fitted to the null distribution of test statistics.)
#' }
#'
#' When \code{output_amount} is set to 3, the resampling distribution of the test statistics also is returned, stored in the columns \code{z_null_1}, \code{z_null_2}, \code{z_null_3}, ... .
#'
#' @details A typical \code{sceptre} analysis consists of three main steps. First, the user sets up the analysis, which entails specifying (i) the response expression matrix, (ii) the gRNA expression matrix, (iii) the cell covariate data frame, (iv) the gRNA group data frame, (v) the formula object, and (vi) the set of response-gRNA group pairs to test for association. Second, the user carries out the calibration check by passing the above objects to the function \code{run_sceptre}, setting \code{calibration_check} to \code{TRUE}. In carrying out the calibration check, \code{sceptre} automatically constructs a set of negative control response-gRNA group pairs to test for association; these pairs are "matched" to the discovery pairs in several respects. The user verifies calibration of the negative control p-values by plotting the calibration results via \code{plot_calibration_result}. Finally, the user runs the discovery analysis by again calling the function \code{run_sceptre}, this time setting \code{calibration_check} to \code{FALSE}. The user visualizes the discovery results via calls to \code{compare_calibration_and_discovery_results} and \code{make_volcano_plot} and obtains the final set of discoveries via a call to \code{obtain_discovery_set}.
#'
#' @note
#'
#' The notes are arranged roughly in order of most important to least important.
#'
#' \itemize{
#' \item{The function \code{compute_cell_covariates} can be used to help compute the \code{covariate_data_frame}, and the function \code{generate_all_pairs} can be used to help compute \code{response_grna_group_pairs}.}
#' \item{When constructing the \code{formula_object}, it is best practice to log-transform the integer count-based variables (e.g., total UMI count, number of expressed responses).}
#' \item{By default, \code{sceptre} assigns gRNAs to cells for users. Users instead may assign gRNAs to cells themselves, passing a logical (i.e., \code{TRUE}/\code{FALSE})  \code{grna_matrix} to \code{sceptre}. The logical matrix should be a standard (dense) R matrix or a sparse matrix of type \code{lgCMatrix}, \code{lgRMatrix}, \code{lgTMatrix}.} In low MOI each column should contain a single \code{TRUE} value, which indicates which gRNA has infected the given cell.
#' \item{\code{sceptre} tests a given pair for association in three stages. First, \code{sceptre} computes \code{B1} null statistics and calculates an initial empirical p-value \code{p1} using these null statistics. If \code{p1} is promising (i.e., \code{p1} < 0.02), then \code{sceptre} proceeds to stage 2; otherwise, \code{sceptre} returns \code{p1}. In stage 2, \code{sceptre} computes \code{B2} null statistics and fits a parametric curve (by default, a skew normal distribution) to these null statistics. If the fit of the parametric curve is good, then \code{sceptre} returns \code{p2}. Otherwise, \code{sceptre} proceeds to stage 3. In stage 3, \code{sceptre} computes \code{B3} null statistics and calculates an empirical p-value \code{p3} using these null statistics. \code{sceptre} then returns \code{p3}. When \code{resampling_mechanism} is set to "crt" (the default for high MOI data), the stage 3 empirical p-value is computed using the null test statistics from stage 2 (to save compute).}
#' }
#' @keywords internal
run_sceptre <- function(response_matrix, grna_matrix, covariate_data_frame, grna_group_data_frame,
                        moi, formula_object, calibration_check, response_grna_group_pairs = NULL,
                        control_group = "default", resampling_mechanism = "default", side = "both",
                        fit_parametric_curve = TRUE, n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L,
                        grna_assign_threshold = 5L, calibration_group_size = NULL, n_calibration_pairs = NULL,
                        B1 = 499L, B2 = 4999L, B3 = 24999L, print_progress = TRUE, output_amount = 1L) {
  ###############
  # PART 1: SETUP
  ###############
  cat("Running setup. ")
  # 1. check function input arguments
  check_inputs(response_matrix, grna_matrix, covariate_data_frame,
               grna_group_data_frame, formula_object, calibration_check,
               response_grna_group_pairs, moi, control_group, resampling_mechanism, side) |> invisible()

  # 2. harmonize arguments (called for side-effects)
  harmonize_arguments(fit_parametric_curve, moi, control_group, resampling_mechanism, side) |> invisible()

  # 3. make the response matrix row accessible
  response_matrix <- set_matrix_accessibility(response_matrix, make_row_accessible = TRUE)

  # 4. convert the cell covariate data frame into a design matrix
  covariate_matrix <- convert_covariate_df_to_design_matrix(covariate_data_frame, formula_object)
  # rm(covariate_data_frame)

  # 5. assign gRNAs to cells
  grna_assignments <- assign_grnas_to_cells(grna_matrix, grna_group_data_frame, grna_assign_threshold,
                                            low_moi, control_group_complement, calibration_check)
  # rm(grna_matrix)
  cat(crayon::green(' \u2713\n'))

  # 6. construct the negative control pairs
  if (calibration_check) {
    cat("Constructing negative control pairs.")
    response_grna_group_pairs <- construct_negative_control_pairs(n_calibration_pairs, calibration_group_size,
                                                                  grna_assignments, response_matrix, n_nonzero_trt_thresh,
                                                                  n_nonzero_cntrl_thresh, grna_group_data_frame,
                                                                  response_grna_group_pairs, control_group_complement, low_moi)
    cat(crayon::green(' \u2713\n'))
  }

  # 7. generate the set of synthetic indicator idxs
  if (run_permutations) {
    cat("Generating permutation resamples.")
    n_cells <- nrow(covariate_matrix)
    synthetic_idxs <- get_synthetic_permutation_idxs(grna_assignments, B1 + B2 + B3, calibration_check,
                                                     control_group_complement, calibration_group_size, n_cells)
    cat(crayon::green(' \u2713\n'))
  }
  gc() |> invisible()

  ####################
  # PART 2: RUN METHOD
  ####################
  if (run_permutations) {
    ret <- run_perm_test_in_memory(response_matrix, grna_assignments, covariate_matrix, response_grna_group_pairs,
                                   synthetic_idxs, output_amount, fit_parametric_curve, B1, B2, B3, calibration_check,
                                   control_group_complement, n_nonzero_trt_thresh, n_nonzero_cntrl_thresh,
                                   side_code, low_moi, print_progress)
  } else {
    ret <- run_crt_in_memory_v2(response_matrix, grna_assignments, covariate_matrix, response_grna_group_pairs,
                                output_amount, fit_parametric_curve, B1, B2, B3, calibration_check, control_group_complement,
                                n_nonzero_trt_thresh, n_nonzero_cntrl_thresh, side_code, low_moi, print_progress)
  }
  return(ret)
}

#' Run \code{sceptre} (low MOI)
#'
#' @inherit run_sceptre
#' @examples
#' \dontrun{
#' library(Matrix)
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
#' calibration_result <- run_sceptre_lowmoi(response_matrix = response_matrix_lowmoi,
#' grna_matrix = grna_matrix_lowmoi,
#' covariate_data_frame = covariate_data_frame_lowmoi,
#' grna_group_data_frame = grna_group_data_frame_lowmoi,
#' formula_object = formula_object,
#' response_grna_group_pairs = response_grna_group_pairs,
#' calibration_check = TRUE)
#'
#' # 4. plot the calibration result to ensure adequate calibration
#' plot_calibration_result(calibration_result)
#'
#' # 5. run the discovery analysis (NOTE: calibration_check set to FALSE)
#' discovery_result <- run_sceptre_lowmoi(response_matrix = response_matrix_lowmoi,
#' grna_matrix = grna_matrix_lowmoi,
#' covariate_data_frame = covariate_data_frame_lowmoi,
#' grna_group_data_frame = grna_group_data_frame_lowmoi,
#' formula_object = formula_object,
#' response_grna_group_pairs = response_grna_group_pairs,
#' calibration_check = FALSE)
#'
#' # 6. compare discovery p-values to the negative control p-values; make a volcano plot
#' compare_calibration_and_discovery_results(calibration_result, discovery_result)
#' make_volcano_plot(discovery_result)
#'
#' # 7. obtain the discovery set for downstream analysis
#' discovery_set <- obtain_discovery_set(discovery_result)
#' }
run_sceptre_lowmoi <- function(response_matrix,
           grna_matrix,
           covariate_data_frame,
           grna_group_data_frame,
           response_grna_group_pairs,
           formula_object,
           calibration_check,
           control_group = "nt_cells",
           resampling_mechanism = "permutations",
           n_nonzero_trt_thresh = 7L,
           n_nonzero_cntrl_thresh = 7L,
           side = "both",
           output_amount = 1L,
           fit_parametric_curve = TRUE,
           calibration_group_size = NULL,
           n_calibration_pairs = NULL,
           B1 = 499L,
           B2 = 4999L,
           B3 = 24999L,
           print_progress = TRUE) {
    run_sceptre(
      response_matrix = response_matrix,
      grna_matrix = grna_matrix,
      covariate_data_frame = covariate_data_frame,
      grna_group_data_frame = grna_group_data_frame,
      response_grna_group_pairs = response_grna_group_pairs,
      moi = "low",
      side = side,
      formula_object = formula_object,
      calibration_check = calibration_check,
      control_group = control_group,
      resampling_mechanism = resampling_mechanism,
      n_nonzero_trt_thresh = n_nonzero_trt_thresh,
      n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
      grna_assign_threshold = NA,
      output_amount = output_amount,
      fit_parametric_curve = fit_parametric_curve,
      calibration_group_size = calibration_group_size,
      n_calibration_pairs = n_calibration_pairs,
      B1 = B1,
      B2 = B2,
      B3 = B3
    )
  }


#' Run \code{sceptre} (high MOI; experimental)
#'
#' This is the new, experimental high MOI function. The user interface of this function is nearly identical to that of the low MOI function. Therefore, we recommend that users work through the [low MOI tutorial](https://katsevich-lab.github.io/sceptre/articles/lowmoi_tutorial.html) before using this function.
#'
#' @export
#' @inherit run_sceptre
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#'
#' # 0. load the data
#' data(response_matrix_highmoi_experimental)       # response-by-cell expression matrix
#' data(grna_matrix_highmoi_experimental)           # gRNA-by-cell expression matrix
#' data(covariate_data_frame_highmoi_experimental)  # cell-by-covariate data frame
#' data(grna_group_data_frame_highmoi_experimental) # gRNA group information
#' data(discovery_pairs_highmoi_experimental)       # discovery gRNA group / gene pairs
#'
#' # 1. set the formula object
#' formula_object <- formula(~log(grna_n_umis) + log(grna_n_nonzero) + log(gene_n_umis) + log(gene_n_nonzero) + batch + p_mito)
#'
#' # 2. run the calibration check
#' calibration_result <- run_sceptre_highmoi_experimental(
#'   response_matrix = response_matrix_highmoi_experimental,
#'   grna_matrix = grna_matrix_highmoi_experimental,
#'   covariate_data_frame = covariate_data_frame_highmoi_experimental,
#'   grna_group_data_frame = grna_group_data_frame_highmoi_experimental,
#'   response_grna_group_pairs = discovery_pairs_highmoi_experimental,
#'   formula_object = formula_object,
#'   side = "left",
#'   calibration_check = TRUE
#' )
#'
#' # 3. verify calibration
#' plot_calibration_result(calibration_result)
#'
#' # 4. run discovery analysis
#' discovery_result <- run_sceptre_highmoi_experimental(
#'   response_matrix = response_matrix_highmoi_experimental,
#'   grna_matrix = grna_matrix_highmoi_experimental,
#'   covariate_data_frame = covariate_data_frame_highmoi_experimental,
#'   grna_group_data_frame = grna_group_data_frame_highmoi_experimental,
#'   response_grna_group_pairs = discovery_pairs_highmoi_experimental,
#'   formula_object = formula_object,
#'   side = "left",
#'   calibration_check = FALSE
#' )
#'
#' # 5. compare discovery p-values to the negative control p-values; make a volcano plot
#' compare_calibration_and_discovery_results(calibration_result, discovery_result)
#' make_volcano_plot(discovery_result)
#'
#' # 6. obtain the discovery set
#' discovery_set <- obtain_discovery_set(discovery_result)
#'
#' # 7. optional: positive control analysis
#' pc_result <- run_sceptre_highmoi_experimental(
#'   response_matrix = response_matrix_highmoi_experimental,
#'   grna_matrix = grna_matrix_highmoi_experimental,
#'   covariate_data_frame = covariate_data_frame_highmoi_experimental,
#'   grna_group_data_frame = grna_group_data_frame_highmoi_experimental,
#'   response_grna_group_pairs = pc_pairs_highmoi_experimental,
#'   formula_object = formula_object,
#'   side = "left",
#'   calibration_check = FALSE
#' )
#' }
run_sceptre_highmoi_experimental <- function(response_matrix,
                                             grna_matrix,
                                             covariate_data_frame,
                                             grna_group_data_frame,
                                             formula_object,
                                             calibration_check,
                                             response_grna_group_pairs = NULL,
                                             resampling_mechanism = "default",
                                             n_nonzero_trt_thresh = 7L,
                                             n_nonzero_cntrl_thresh = 7L,
                                             side = "both",
                                             grna_assign_threshold = 5L,
                                             output_amount = 1L,
                                             fit_parametric_curve = TRUE,
                                             calibration_group_size = NULL,
                                             n_calibration_pairs = NULL,
                                             B1 = 499L,
                                             B2 = 4999L,
                                             B3 = 24999L,
                                             print_progress = TRUE

) {
    run_sceptre(
      response_matrix = response_matrix,
      grna_matrix = grna_matrix,
      covariate_data_frame = covariate_data_frame,
      grna_group_data_frame = grna_group_data_frame,
      side = side,
      response_grna_group_pairs = response_grna_group_pairs,
      moi = "high",
      formula_object = formula_object,
      calibration_check = calibration_check,
      control_group = "default",
      resampling_mechanism = resampling_mechanism,
      n_nonzero_trt_thresh = n_nonzero_trt_thresh,
      n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh,
      grna_assign_threshold = grna_assign_threshold,
      output_amount = output_amount,
      fit_parametric_curve = fit_parametric_curve,
      calibration_group_size = calibration_group_size,
      n_calibration_pairs = n_calibration_pairs,
      B1 = B1, B2 = B2, B3 = B3
    )
}
