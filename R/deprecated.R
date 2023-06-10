#' Run \code{sceptre} on low MOI data
#'
#' The function \code{run_sceptre_lowmoi()} runs \code{sceptre} on low multiplicity of infection (MOI) single-cell CRISPR screen data. \code{sceptre} tests for association (and estimates the log fold change in expression) for a user-specified set of response-gRNA group pairs. Applying \code{sceptre} consists of three main steps: (1) setting up the analysis, (2) carrying out a "calibration check" to verify that \code{sceptre} controls the rate of false discoveries on negative control data, and (3) conducting the "discovery analysis" to discover regulatory relationships.
#'
#' @param response_matrix (required) a matrix of raw expression counts. The responses (e.g., genes or proteins) should be in the rows, and the cells should be in the columns. The row names should be the unique IDs of the responses. The matrix can be a standard (dense) matrix or a sparse matrix of class \code{dgCMatrix}, \code{dgRMatrix}, or \code{dgTMatrix}.
#' @param grna_matrix (required) a matrix of gRNA expression counts. The gRNAs should be in the rows, and the cells should be in the columns. The row names should be the unique IDs of the gRNAs. The matrix can be a standard (dense) matrix or a sparse matrix of class \code{dgCMatrix}, \code{dgRMatrix}, or \code{dgTMatrix}. (See note about passing logical gRNA matrices.)
#' @param covariate_data_frame (required) a data frame containing the cell-specific covariates (e.g., sequencing batch, total UMI count, etc.), with cells in the rows and covariates in the columns. The column names should be the names of the covariates.
#' @param grna_group_data_frame (required) a data frame specifying the group to which each individual gRNA belongs. This data frame should contain columns \code{grna_id} and \code{grna_group}, with \code{grna_group} specifying the group membership of a given gRNA in the \code{grna_id} column. Typically, gRNAs that target the same site are grouped into the same gRNA group. IMPORTANT: Non-targeting gRNAs should be assigned a gRNA group label of "non-targeting".
#' @param formula_object (required) an R formula object specifying how \code{sceptre} is to adjust for the covariates (in the \code{covariate_data_frame}).
#' @param response_grna_group_pairs (required) a data frame specifying the set of response-gRNA group pairs to test for association. The data frame should contain columns \code{response_id} and \code{grna_group}.
#' @param calibration_check (required) a logical (i.e., \code{TRUE}/\code{FALSE}) value indicating whether to run a calibration check (\code{TRUE}) or a discovery analysis (\code{FALSE}). See the "Details" section for a more complete description of the calibration and discovery analyses.
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
#' @param discovery_test_stat (optional; lower importance; default "exact") the test statistic to use in the discovery analysis, either "exact" or "approximate".
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
#' If \code{return_resampling_dist} is set to \code{TRUE}, the null permutation distribution is returned, stored in columns named "z_{null_1}", ..., "z_{null_{B1}}". If \code{return_debugging_metrics} is set to \code{TRUE}, the additional columns \code{sn_fit_used}, \code{round}, and \code{regression_method} are included in the ouput data frame. These column are described in the "Notes" section.
#'
#' @export
#'
#' @details A typical \code{sceptre} analysis consists of three main steps. First, the user sets up the analysis, which entails specifying (i) the response expression matrix, (ii) the gRNA expression matrix, (iii) the cell covariate data frame, (iv) the gRNA group data frame, (v) the formula object, and (vi) the set of response-gRNA group pairs to test for association. Second, the user carries out the calibration check by passing the above objects to the \code{run_sceptre_lowmoi} function, setting \code{calibration_check} to \code{TRUE}. In carrying out the calibration check, \code{sceptre} automatically constructs a set of negative control response-gRNA group pairs to test for association. The negative control pairs are "matched" to the discovery pairs in three ways: (i) the negative control pairs are subjected to the same pair-wise QC as the discovery pairs, (ii) the number of negative control pairs equals the number of discovery pairs (after applying pair-wise QC), and (iii) the negative control pairs contain the same number of gRNAs per gRNA group as the discovery pairs. After running \code{run_sceptre_lowmoi}, the user verifies calibration of the p-values by plotting the calibration results via \code{plot_calibration_result}. Finally, the user runs a discovery analysis by again calling \code{run_sceptre_lowmoi}, this time setting \code{calibration_check} to \code{FALSE}. The user visualizes the discovery results via calls to \code{compare_calibration_and_discovery_results} and \code{make_volcano_plot} and obtains the final set of discoveries via a call to \code{obtain_discovery_set}.
#'
#' @note
#'
#' The notes are arranged roughly in order of most important to least important.
#'
#' \itemize{
#' \item{The function \code{compute_cell_covariates} can be used to help compute the \code{covariate_data_frame}, and the function \code{generate_all_pairs} can be used to help compute \code{response_grna_group_pairs}.}
#' \item{When constructing the \code{formula_object}, it is best practice to log-transform the integer count-based variables (e.g., total UMI count, number of expressed responses).}
#' \item{\code{grna_matrix} optionally can be a logical (i.e., \code{TRUE}/\code{FALSE}) matrix of cell-to-gRNA assignments. Each column should contain a single \code{TRUE} value, which indicates which gRNA has infected the given cell. The logical matrix should be a standard (dense) R matrix or a sparse matrix of type \code{lgCMatrix}, \code{lgRMatrix}, \code{lgTMatrix}.}
#' \item{\code{sceptre} uses three rounds of permutation testing to test a given pair. First, \code{sceptre} computes \code{B1} null statistics and calculates an initial empirical p-value \code{p1} using these null statistics. If \code{p1} is promising (i.e., \code{p1} < 0.02), then \code{sceptre} proceeds to round 2; otherwise, \code{sceptre} returns \code{p1}. In round 2, \code{sceptre} computes \code{B2} null statistics and fits a skew normal distribution to these null statistics. If the skew-normal fit is good, then \code{sceptre} returns \code{p2}. Otherwise, \code{sceptre} proceeds to round 3. In round 3, \code{sceptre} computes \code{B3} null statistics and calculates an empirical p-value \code{p3} using these null statistics. \code{sceptre} then returns \code{p3}.}
#' \item{When \code{return_debugging_metrics} is set to \code{TRUE}, the columns \code{sn_fit_used}, \code{round}, and \code{regression_method} are included in the output data frame. \code{sn_fit_used} indicates whether a p-value was calculated via the skew-normal distribution; \code{sn_fit_used} is \code{FALSE} if (i) \code{sceptre} computes the p-value in round 1 or round 3 or if (ii) \code{fit_skew_normal} is set to \code{FALSE} by the user. \code{round} indicates the round (among rounds 1, 2, and 3) in which the p-value was computed. Finally, \code{regression_method} is a string indicating the method that was used to regress the gene expression vector onto the cell covariate matrix and estimate the negative binomial theta parameter.}
#' \item{When \code{return_resampling_dist} is set to \code{TRUE}, \code{sceptre} uses only a single round of permutation testing, returning \code{B1} resampled statistics for each pair. \code{return_resampling_dist} should be set to \code{TRUE} for debugging purposes only.}
#' \item{The \code{discovery_test_stat} parameter controls the test statistic used. The "exact" statistic is more robust, while the "approximate" statistic is faster. In most cases the two test statistics yield similar results.}
#' }
#'
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
run_sceptre_lowmoi <- function(response_matrix, grna_matrix,
                               covariate_data_frame, grna_group_data_frame,
                               formula_object, response_grna_group_pairs,
                               calibration_check, n_nonzero_trt_thresh = 7L,
                               n_nonzero_cntrl_thresh = 7L, discovery_test_stat = "exact",
                               return_debugging_metrics = FALSE, return_resampling_dist = FALSE,
                               fit_skew_normal = TRUE, calibration_group_size = NULL,
                               n_calibration_pairs = NULL, B1 = 499L, B2 = 4999L, B3 = 24999L,
                               regression_method = "poisson_glm", print_progress = TRUE) {
  ###############
  # PART 1: SETUP
  ###############
  cat("Running setup. ")
  # 1. check function input arguments
  check_inputs(response_matrix, grna_matrix, covariate_data_frame, grna_group_data_frame,
               formula_object, calibration_check, response_grna_group_pairs, regression_method) |> invisible()

  # 1.5 order the response grna group pairs
  response_grna_group_pairs <- data.table::as.data.table(response_grna_group_pairs)
  data.table::setorderv(response_grna_group_pairs, cols = "response_id")

  # 2. harmonize arguments (called for side-effect)
  harmonize_arguments(return_resampling_dist, fit_skew_normal) |> invisible()

  # 3. cast and transpose response matrix; cast grna matrix
  response_matrix <- set_matrix_accessibility(response_matrix, TRUE)

  # 4. convert the cell covariate data frame into a design matrix
  covariate_matrix <- convert_covariate_df_to_design_matrix(covariate_data_frame, formula_object)
  rm(covariate_data_frame)

  # 5. assign the gRNAs to cells
  grna_assignments <- assign_grnas_to_cells_lowmoi_v2(grna_matrix, grna_group_data_frame, calibration_check, n_calibration_pairs)
  rm(grna_matrix)
  cat(crayon::green(' \u2713\n'))

  # 6. construct the undercover response_grna_group_pairs
  if (calibration_check) {
    cat("Constructing negative control pairs.")
    if (is.null(calibration_group_size)) calibration_group_size <- compute_calibration_group_size(grna_group_data_frame)
    response_grna_group_pairs <- construct_negative_control_pairs(n_calibration_pairs, calibration_group_size, grna_assignments, response_matrix, n_nonzero_trt_thresh, n_nonzero_cntrl_thresh, grna_group_data_frame, response_grna_group_pairs)
    cat(crayon::green(' \u2713\n'))
  }

  # 7. generate the set of synthetic indicator idxs
  cat("Generating permutation resamples.")
  synthetic_idxs <- get_synthetic_idxs_lowmoi(grna_assignments, B1 + B2 + B3, calibration_check, calibration_group_size)
  cat(crayon::green(' \u2713\n'))
  gc() |> invisible()

  ####################
  # PART 2: RUN METHOD
  ####################
  cat("Running differential expression analyses.\n")
  ret <- run_lowmoi_in_memory(response_matrix, grna_assignments,
                              covariate_matrix, response_grna_group_pairs,
                              synthetic_idxs, return_resampling_dist, fit_skew_normal,
                              B1, B2, B3, calibration_check, discovery_test_stat,
                              n_nonzero_trt_thresh, n_nonzero_cntrl_thresh,
                              return_debugging_metrics, regression_method, print_progress)
}


# data(covariate_matrix_highmoi)
# data(gene_matrix_highmoi)
# data(gRNA_matrix_highmoi)
# data(gRNA_groups_table_highmoi)
#'
# # modify the grna group table to get it into the "standard" format
# gRNA_groups_table_highmoi$gRNA_group[gRNA_groups_table_highmoi$gRNA_type == "non_target"] <- "non-targeting"
# gRNA_groups_table_highmoi <- gRNA_groups_table_highmoi |>
# dplyr::rename("grna_id" = "gRNA_id", "grna_group" = "gRNA_group") |>
# dplyr::mutate(gRNA_type = NULL)
#'
# response_matrix <- gene_matrix_highmoi
# grna_matrix <- gRNA_matrix_highmoi
# covariate_data_frame <- covariate_matrix_highmoi
# grna_group_data_frame <- gRNA_groups_table_highmoi
# formula_object <- formula(~ lg_gRNA_lib_size + lg_gene_lib_size + p_mito + batch)
# response_grna_group_pairs <- generate_all_pairs(response_matrix, grna_group_data_frame)
# calibration_check <- FALSE
.run_sceptre_highmoi <- function(response_matrix, grna_matrix,
                                covariate_data_frame, grna_group_data_frame,
                                formula_object, response_grna_group_pairs,
                                calibration_check, threshold = 5, n_nonzero_trt_thresh = 7L,
                                n_nonzero_cntrl_thresh = 7L, discovery_test_stat = "exact",
                                return_debugging_metrics = FALSE, return_resampling_dist = FALSE,
                                fit_skew_normal = TRUE, calibration_group_size = NULL,
                                n_calibration_pairs = NULL, B1 = 499L, B2 = 4999L, B3 = 24999L,
                                regression_method = "poisson_glm", print_progress = TRUE) {
  cat("Running setup. ")
  # 1. check function input arguments
  check_inputs(response_matrix, grna_matrix, covariate_data_frame, grna_group_data_frame,
               formula_object, calibration_check, response_grna_group_pairs, regression_method) |> invisible()

  # 2. order the response grna group pairs; cast grna group table to data.table
  response_grna_group_pairs <- data.table::as.data.table(response_grna_group_pairs)
  data.table::setorderv(response_grna_group_pairs, cols = "response_id")
  data.table::setDT(grna_group_data_frame)

  # 3. harmonize arguments (called for side-effects)
  harmonize_arguments(return_resampling_dist, fit_skew_normal) |> invisible()

  # 4. cast and transpose response matrix; cast grna matrix
  response_matrix <- set_matrix_accessibility(response_matrix, TRUE)

  # 5. convert the cell covariate data frame into a design matrix
  covariate_matrix <- convert_covariate_df_to_design_matrix(covariate_data_frame, formula_object)
  rm(covariate_data_frame)

  # 6. assign the gRNAs to cells
  grna_assignments <- assign_grnas_to_cells_highmoi(grna_matrix, threshold, grna_group_data_frame)
  # grna_matrix <- set_matrix_accessibility(grna_matrix, TRUE)
  cat(crayon::green(' \u2713\n'))

  # 7. construct the undercover response_grna_group_pairs (TO DO)
  # if (calibration_check) {
  #  cat("Constructing negative control pairs.")
  #  if (is.null(calibration_group_size)) calibration_group_size <- compute_calibration_group_size(grna_group_data_frame)
  #  response_grna_group_pairs <- construct_negative_control_pairs(n_calibration_pairs, calibration_group_size, grna_assignments, response_matrix, n_nonzero_trt_thresh, n_nonzero_cntrl_thresh, grna_group_data_frame, response_grna_group_pairs)
  #  cat(crayon::green(' \u2713\n'))
  # }

  # 8. generate the set of synthetic indicator idxs
  cat("Generating permutation resamples.")
  synthetic_idxs <- get_synthetic_idxs_highmoi(B = B1 + B2 + B3, grna_assignments = grna_assignments, n_cells = ncol(response_matrix))
  cat(crayon::green(' \u2713\n'))
  gc() |> invisible()

  ####################
  # PART 2: RUN METHOD
  ####################
  cat("Running differential expression analyses.\n")
  ret <- run_high_moi_in_memory(response_matrix, grna_assignments,
                                covariate_matrix, response_grna_group_pairs,
                                synthetic_idxs, return_resampling_dist, fit_skew_normal,
                                B1, B2, B3, calibration_check, discovery_test_stat,
                                n_nonzero_trt_thresh, n_nonzero_cntrl_thresh,
                                return_debugging_metrics, regression_method, print_progress)
}

.run_lowmoi_in_memory <- function(response_matrix, grna_assignments,
                                 covariate_matrix, response_grna_group_pairs,
                                 synthetic_idxs, return_resampling_dist, fit_skew_normal,
                                 B1, B2, B3, calibration_check, discovery_test_stat,
                                 n_nonzero_trt_thresh, n_nonzero_cntrl_thresh,
                                 return_debugging_metrics, regression_method, print_progress) {
  # 0. preliminary setup; initialize the args_to_pass, set the low_level_association_funct
  result_list_outer <- vector(mode = "list", length = 2 * length(unique(response_grna_group_pairs$response_id)))
  out_counter <- 1L
  args_to_pass <- list(synthetic_idxs = synthetic_idxs, B1 = B1, B2 = B2, B3 = B3,
                       fit_skew_normal = fit_skew_normal,
                       return_resampling_dist = return_resampling_dist,
                       grna_group_idxs = grna_assignments$grna_group_idxs,
                       covariate_matrix = covariate_matrix,
                       all_nt_idxs = grna_assignments$all_nt_idxs,
                       regression_method = regression_method,
                       indiv_nt_grna_idxs = if (calibration_check) grna_assignments$indiv_nt_grna_idxs else NA)
  if (calibration_check) {
    low_level_association_funct <- "lowmoi_undercover_stat"
  } else {
    low_level_association_funct <- if (discovery_test_stat == "exact") "lowmoi_exact_stat_discovery" else "lowmoi_approximate_stat_discovery"
  }
  run_outer_regression <- calibration_check || (discovery_test_stat == "approximate")

  # 1. obtain the subset of the covariate matrix corresponding to the NT cells and n_cells
  covariate_matrix_nt <- covariate_matrix[grna_assignments$all_nt_idxs,]
  n_cells <- ncol(response_matrix)

  # 2. loop over the response IDs
  response_ids <- unique(response_grna_group_pairs$response_id)
  for (response_idx in seq_along(response_ids)) {
    if ((response_idx == 1 || response_idx %% 5 == 0) && print_progress) {
      cat(paste0("Analyzing pairs containing response ", as.character(response_ids[response_idx]), " (", response_idx, " of ", length(response_ids), ")\n"))
    }
    if (response_idx %% 200 == 0) gc() |> invisible()
    response_id <- as.character(response_ids[response_idx])

    # 3. load the expressions of the current response; also get the nt expression vector
    expression_vector <- load_csr_row(j = response_matrix@j,
                                      p = response_matrix@p,
                                      x = response_matrix@x,
                                      row_idx = which(rownames(response_matrix) == response_id),
                                      n_cells = n_cells)
    expression_vector_nt <- expression_vector[grna_assignments$all_nt_idxs]

    # 4. obtain the gRNA groups to analyze
    l <- response_grna_group_pairs$response_id == response_id
    curr_df <- response_grna_group_pairs[l,]

    # 5. if running a discovery analysis, do QC
    if (!calibration_check) {
      n_nonzero_cntrl_curr <- sum(expression_vector_nt >= 1)
      grna_group_posits <- match(x = curr_df$grna_group, table = names(grna_assignments$grna_group_idxs))
      n_nonzero_trt_curr <- compute_n_nonzero_trt_vector(expression_vector = expression_vector,
                                                         grna_group_idxs = grna_assignments$grna_group_idxs,
                                                         grna_group_posits = grna_group_posits)
      curr_df$n_nonzero_trt <- n_nonzero_trt_curr
      curr_df$n_nonzero_cntrl <- n_nonzero_cntrl_curr
      pass_qc <- n_nonzero_trt_curr >= n_nonzero_trt_thresh

      # i. if n_nonzero_cntrl_curr is less than n_nonzero_cntrl_thresh, jump to next iteration
      if (n_nonzero_cntrl_curr < n_nonzero_cntrl_thresh || !any(pass_qc)) {
        result_list_outer[[out_counter]] <- curr_df
        out_counter <- out_counter + 1L
        next
      }

      # ii. remove any rows that have not passed qc; keep the rows that have passed qc
      if (!all(pass_qc)) {
        result_list_outer[[out_counter]] <- curr_df[!pass_qc,]
        out_counter <- out_counter + 1L
      }
      curr_df <- curr_df[pass_qc,]
    }


    # 6. perform the expression on technical factor regression
    if (run_outer_regression) {
      response_precomp <- perform_response_precomputation(expressions = expression_vector_nt,
                                                          covariate_matrix = covariate_matrix_nt,
                                                          regression_method = regression_method)

      # 7. obtain precomputation peices for NT cells
      pieces_precomp <- compute_precomputation_pieces(expression_vector_nt,
                                                      covariate_matrix_nt,
                                                      response_precomp$fitted_coefs,
                                                      response_precomp$theta,
                                                      full_test_stat = TRUE)
    } else {
      response_precomp <- pieces_precomp <- NA
    }

    # 8. update args_to_pass
    args_to_pass$grna_groups <- as.character(curr_df$grna_group)
    args_to_pass$pieces_precomp <- pieces_precomp
    args_to_pass$expression_vector_nt <- expression_vector_nt
    args_to_pass$expression_vector <- expression_vector
    args_to_pass$response_precomp <- response_precomp

    # 9. pass the arguments to the appropriate low-level association testing function
    curr_response_result <- do.call(what = low_level_association_funct, args = args_to_pass)

    # 10. combine the response-wise results into a data table; insert into list
    result_list_outer[[out_counter]] <- construct_data_frame_v2(curr_df, curr_response_result,
                                                                return_debugging_metrics, return_resampling_dist,
                                                                response_precomp$precomp_str, low_level_association_funct)
    out_counter <- out_counter + 1L
  }

  # combine and sort result
  ret <- data.table::rbindlist(result_list_outer, fill = TRUE)
  data.table::setorderv(ret, cols = c("p_value", "response_id"), na.last = TRUE)
  return(ret)
}


.run_high_moi_in_memory <- function(response_matrix, grna_assignments,
                                   covariate_matrix, response_grna_group_pairs,
                                   synthetic_idxs, return_resampling_dist, fit_skew_normal,
                                   B1, B2, B3, calibration_check, discovery_test_stat,
                                   n_nonzero_trt_thresh, n_nonzero_cntrl_thresh,
                                   return_debugging_metrics, regression_method, print_progress) {
  # 0. preliminary setup; initialize the args_to_pass, set the low_level_association_funct
  result_list_outer <- vector(mode = "list", length = 2 * length(unique(response_grna_group_pairs$response_id)))
  out_counter <- 1L
  args_to_pass <- list(synthetic_idxs = synthetic_idxs, B1 = B1, B2 = B2, B3 = B3,
                       fit_skew_normal = fit_skew_normal,
                       return_resampling_dist = return_resampling_dist,
                       grna_group_idxs = grna_assignments,
                       covariate_matrix = covariate_matrix)
  n_cells <- ncol(response_matrix)
  low_level_association_funct <- "highmoi_discovery_stat"

  # 2. loop over the response IDs
  response_ids <- unique(response_grna_group_pairs$response_id)
  for (response_idx in seq_along(response_ids)) {
    if ((response_idx == 1 || response_idx %% 5 == 0) && print_progress) {
      cat(paste0("Analyzing pairs containing response ", as.character(response_ids[response_idx]), " (", response_idx, " of ", length(response_ids), ")\n"))
    }
    if (response_idx %% 200 == 0) gc() |> invisible()
    response_id <- as.character(response_ids[response_idx])

    # 3. load the expressions of the current response; also get the nt expression vector
    expression_vector <- load_csr_row(j = response_matrix@j,
                                      p = response_matrix@p,
                                      x = response_matrix@x,
                                      row_idx = which(rownames(response_matrix) == response_id),
                                      n_cells = n_cells)

    # 4. obtain the gRNA groups to analyze
    l <- response_grna_group_pairs$response_id == response_id
    curr_df <- response_grna_group_pairs[l,]

    # 5. if running a discovery analysis, do QC
    if (!calibration_check) {
      n_nonzero_tot <- sum(expression_vector >= 1)
      grna_group_posits <- match(x = curr_df$grna_group, table = names(grna_assignments))
      n_nonzero_trt_curr <- compute_n_nonzero_trt_vector(expression_vector = expression_vector,
                                                         grna_group_idxs = grna_assignments,
                                                         grna_group_posits = grna_group_posits)
      curr_df$n_nonzero_trt <- n_nonzero_trt_curr
      curr_df$n_nonzero_cntrl <- n_nonzero_tot - n_nonzero_trt_curr
      pass_qc <- (curr_df$n_nonzero_trt >= n_nonzero_trt_thresh) & (curr_df$n_nonzero_cntrl >= n_nonzero_cntrl_thresh)

      # i. if no pair passes qc, jump to next iteration
      if (!any(pass_qc)) {
        result_list_outer[[out_counter]] <- curr_df
        out_counter <- out_counter + 1L
        next
      }

      # ii. remove any rows that have not passed qc; keep the rows that have passed qc
      if (!all(pass_qc)) {
        result_list_outer[[out_counter]] <- curr_df[!pass_qc,]
        out_counter <- out_counter + 1L
      }
      curr_df <- curr_df[pass_qc,]
    }

    # 6. perform the gene expression on technical factor regression
    response_precomp <- perform_response_precomputation(expressions = expression_vector,
                                                        covariate_matrix = covariate_matrix,
                                                        regression_method = regression_method)

    # 7. obtain precomputation peices for NT cells
    pieces_precomp <- compute_precomputation_pieces(expression_vector,
                                                    covariate_matrix,
                                                    response_precomp$fitted_coefs,
                                                    response_precomp$theta,
                                                    full_test_stat = TRUE)

    # 8. update args to pass
    args_to_pass$grna_groups <- as.character(curr_df$grna_group)
    args_to_pass$pieces_precomp <- pieces_precomp
    args_to_pass$expression_vector <- expression_vector

    # pass the arguments the the low-level association testing function
    curr_response_result <- do.call(what = low_level_association_funct, args = args_to_pass)

    # 10. combine the response-wise results into a data table; insert into list
    result_list_outer[[out_counter]] <- construct_data_frame_v2(curr_df, curr_response_result,
                                                                return_debugging_metrics, return_resampling_dist,
                                                                response_precomp$precomp_str, low_level_association_funct)
    out_counter <- out_counter + 1L
  }

  # combine and sort result
  ret <- data.table::rbindlist(result_list_outer, fill = TRUE)
  data.table::setorderv(ret, cols = c("p_value", "response_id"), na.last = TRUE)
  return(ret)
}

.highmoi_discovery_stat <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, grna_group_idxs, covariate_matrix, grna_groups, pieces_precomp, expression_vector) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  n_cells <- length(expression_vector)
  # compute the D matrix
  D <- compute_D_matrix(Zt_wZ = pieces_precomp$Zt_wZ, wZ = pieces_precomp$wZ)
  for (i in seq_along(grna_groups)) {
    curr_grna_group <- grna_groups[i]
    idxs <- get_discovery_idx_vector(curr_grna_group, grna_group_idxs, n_cells)
    result <- run_low_level_test_full_v2(y = expression_vector[idxs$idxs],
                                         mu = pieces_precomp$mu[idxs$idxs],
                                         a = pieces_precomp$a[idxs$idxs],
                                         w = pieces_precomp$w[idxs$idxs],
                                         D = D[,idxs$idxs],
                                         n_cntrl = idxs$n_cntrl,
                                         n_trt = idxs$n_trt,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}


.lowmoi_approximate_stat_discovery <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, grna_group_idxs, covariate_matrix, all_nt_idxs, regression_method, indiv_nt_grna_idxs, grna_groups, pieces_precomp, expression_vector_nt, expression_vector, response_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    grna_group <- grna_groups[i]
    expression_vector_trt <- expression_vector[grna_group_idxs[[grna_group]]]
    covariate_matrix_trt <- covariate_matrix[grna_group_idxs[[grna_group]],]

    # 1. compute the precomputation pieces for the treatment cells
    pieces_trt <- compute_precomputation_pieces(expression_vector_trt, covariate_matrix_trt, response_precomp$fitted_coefs, response_precomp$theta, TRUE)

    # 2. compute the shared weighted covariate matrix
    Zt_wZ <- pieces_trt$Zt_wZ + pieces_precomp$Zt_wZ

    # 3. compute the D matrix
    D <- compute_D_matrix(Zt_wZ = Zt_wZ, wZ = rbind(pieces_precomp$wZ, pieces_trt$wZ))

    # 4. create the combined mu, y, and a vectors
    y <- c(expression_vector_nt, expression_vector_trt)
    mu <- c(pieces_precomp$mu, pieces_trt$mu)
    a <- c(pieces_precomp$a, pieces_trt$a)
    w <- c(pieces_precomp$w, pieces_trt$w)

    # 5. call low-level estimation and testing function
    n_cntrl <- length(expression_vector_nt)
    n_trt <- length(expression_vector_trt)
    result <- run_low_level_test_full_v2(y = y,
                                         mu = mu,
                                         a = a,
                                         w = w,
                                         D = D,
                                         n_cntrl = n_cntrl,
                                         n_trt = n_trt,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}


.lowmoi_undercover_stat <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, grna_group_idxs, covariate_matrix, all_nt_idxs, regression_method, indiv_nt_grna_idxs, grna_groups, pieces_precomp, expression_vector_nt, expression_vector, response_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  # compute the D matrix
  D <- compute_D_matrix(Zt_wZ = pieces_precomp$Zt_wZ, wZ = pieces_precomp$wZ)
  for (i in seq_along(grna_groups)) {
    undercover_group <- grna_groups[i]
    idxs <- get_undercover_idx_vector(undercover_group, indiv_nt_grna_idxs)
    result <- run_low_level_test_full_v2(y = expression_vector_nt[idxs$idxs],
                                         mu = pieces_precomp$mu[idxs$idxs],
                                         a = pieces_precomp$a[idxs$idxs],
                                         w = pieces_precomp$w[idxs$idxs],
                                         D = D[,idxs$idxs],
                                         n_cntrl = idxs$n_cntrl,
                                         n_trt = idxs$n_trt,
                                         synthetic_idxs = synthetic_idxs,
                                         B1 = B1, B2 = B2, B3 = B3,
                                         fit_skew_normal = fit_skew_normal,
                                         return_resampling_dist = return_resampling_dist)
    result_list_inner[[i]] <- result
  }
  return(result_list_inner)
}


.lowmoi_exact_stat_discovery <- function(synthetic_idxs, B1, B2, B3, fit_skew_normal, return_resampling_dist, grna_group_idxs, covariate_matrix, all_nt_idxs, regression_method, indiv_nt_grna_idxs, grna_groups, pieces_precomp, expression_vector_nt, expression_vector, response_precomp) {
  result_list_inner <- vector(mode = "list", length = length(grna_groups))
  for (i in seq_along(grna_groups)) {
    grna_group <- grna_groups[i]
    trt_idxs <- grna_group_idxs[[grna_group]]
    n_trt <- length(trt_idxs)
    n_cntrl <- length(all_nt_idxs)
    combined_idxs <- c(all_nt_idxs, trt_idxs)

    # 1. get the expression vector and covariate matrix
    curr_expression_vector <- expression_vector[combined_idxs]
    curr_covariate_matrix <- covariate_matrix[combined_idxs,]

    # 2. perform the response precomputation
    response_precomp <- perform_response_precomputation(expressions = curr_expression_vector,
                                                        covariate_matrix = curr_covariate_matrix,
                                                        regression_method = regression_method)
    result_list_inner[[i]] <- tryCatch({
      # 2. get the precomp pieces
      precomp_pieces <- compute_precomputation_pieces(expression_vector = curr_expression_vector,
                                                      covariate_matrix = curr_covariate_matrix,
                                                      fitted_coefs = response_precomp$fitted_coefs,
                                                      theta = response_precomp$theta,
                                                      full_test_stat = TRUE)
      # 3. compute D matrix
      D <- compute_D_matrix(precomp_pieces$Zt_wZ, precomp_pieces$wZ)

      # 4. run the association test
      run_low_level_test_full_v2(y = curr_expression_vector,
                                 mu = precomp_pieces$mu,
                                 a = precomp_pieces$a,
                                 w = precomp_pieces$w,
                                 D = D,
                                 n_cntrl = n_cntrl,
                                 n_trt = n_trt,
                                 synthetic_idxs = synthetic_idxs,
                                 B1 = B1, B2 = B2, B3 = B3,
                                 fit_skew_normal = fit_skew_normal,
                                 return_resampling_dist = return_resampling_dist)
    }, error = function(e) backup_distilled(curr_expression_vector, curr_covariate_matrix,
                                            response_precomp, n_cntrl, n_trt, synthetic_idxs,
                                            B1, B2, B3, fit_skew_normal, return_resampling_dist),
    warning = function(w) backup_distilled(curr_expression_vector, curr_covariate_matrix,
                                           response_precomp, n_cntrl, n_trt, synthetic_idxs,
                                           B1, B2, B3, fit_skew_normal, return_resampling_dist))
  }
  return(result_list_inner)
}


.backup_distilled <- function(curr_expression_vector, curr_covariate_matrix,
                             response_precomp, n_cntrl, n_trt, synthetic_idxs,
                             B1, B2, B3, fit_skew_normal, return_resampling_dist) {
  # 2. get the precomp pieces
  precomp_pieces <- compute_precomputation_pieces(expression_vector = curr_expression_vector,
                                                  covariate_matrix = curr_covariate_matrix,
                                                  fitted_coefs = response_precomp$fitted_coefs,
                                                  theta = response_precomp$theta,
                                                  full_test_stat = FALSE)
  # 4. run the association test
  run_low_level_test_distilled(y = curr_expression_vector,
                               mu = precomp_pieces$mu,
                               a = precomp_pieces$a,
                               b = precomp_pieces$b,
                               n_cntrl = n_cntrl,
                               n_trt = n_trt,
                               synthetic_idxs = synthetic_idxs,
                               B1 = B1, B2 = B2, B3 = B3,
                               fit_skew_normal = fit_skew_normal,
                               return_resampling_dist = return_resampling_dist)
}


.compute_D_matrix <- function(Zt_wZ, wZ) {
  P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
  U <- P_decomp$vectors
  Lambda_minus_half <- 1/sqrt(P_decomp$values)
  D <- (Lambda_minus_half * t(U)) %*% t(wZ)
  return(D)
}


.get_undercover_idx_vector <- function(undercover_group, indiv_nt_grna_idxs) {
  undercover_nts <- strsplit(x = undercover_group, split = "&", fixed = TRUE)[[1]]
  control_cells <- indiv_nt_grna_idxs[setdiff(names(indiv_nt_grna_idxs), undercover_nts)] |>
    unlist() |> stats::setNames(NULL)
  undercover_cells <- indiv_nt_grna_idxs[undercover_nts] |> unlist() |> stats::setNames(NULL)
  idxs <- c(control_cells, undercover_cells)
  return(list(idxs = idxs, n_trt = length(undercover_cells), n_cntrl = length(control_cells)))
}


.get_discovery_idx_vector <- function(curr_grna_group, grna_group_idxs, n_cells) {
  trt_idxs <- grna_group_idxs[[curr_grna_group]]
  cntrl_idxs <- seq(1L, n_cells)[-trt_idxs]
  idxs <- c(cntrl_idxs, trt_idxs)
  return(list(idxs = idxs, n_trt = length(trt_idxs), n_cntrl = length(cntrl_idxs)))
}
