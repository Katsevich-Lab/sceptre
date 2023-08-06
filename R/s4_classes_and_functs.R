# response matrix class union
#' @import Matrix
setClassUnion("response_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix", "dgTMatrix"))
setClassUnion("grna_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix", "dgTMatrix", "lgCMatrix", "lgRMatrix", "lgTMatrix"))

# sceptre object class
setClass("sceptre_object",
         slots = list(# raw data
                      response_matrix = "response_matrix_class",
                      grna_matrix = "grna_matrix_class",
                      covariate_data_frame = "data.frame",
                      covariate_matrix = "matrix",
                      grna_group_data_frame = "data.frame",
                      low_moi = "logical",
                      user_specified_covariates = "character",

                      # analysis parameters
                      discovery_pairs = "data.frame",
                      positive_control_pairs = "data.frame",
                      formula_object = "formula",
                      side_code = "integer",
                      fit_parametric_curve = "logical",
                      control_group_complement = "logical",
                      run_permutations = "logical",
                      n_nonzero_trt_thresh = "integer",
                      n_nonzero_cntrl_thresh = "integer",
                      B1 = "integer", B2 = "integer", B3 = "integer",
                      grna_assign_threshold = "integer",
                      calibration_group_size = "integer",
                      n_calibration_pairs = "integer",

                      # computed objects
                      M_matrix = "matrix",
                      n_nonzero_tot_vector = "integer",
                      n_ok_discovery_pairs = "integer",
                      n_ok_positive_control_pairs = "integer",
                      discovery_pairs_with_info = "data.frame",
                      positive_control_pairs_with_info = "data.frame",

                      # cached objects
                      response_precomputations = "list",
                      grna_assignments = "list",
                      negative_control_pairs = "data.frame",

                      # flags
                      calibration_check_run = "logical",
                      power_check_run = "logical",
                      discovery_analysis_run = "logical",
                      analysis_prepared = "logical",
                      last_function_called = "character",

                      # results
                      calibration_result = "data.frame",
                      power_result = "data.frame",
                      discovery_result = "data.frame"))

# show method for sceptre class
setMethod("show", signature = signature("sceptre_object"), function(object) {
  n_cells <- ncol(object@response_matrix)
  n_responses <- nrow(object@response_matrix)
  n_nt_grnas <- object@grna_group_data_frame |>
    dplyr::filter(grna_group == "non-targeting") |>
    nrow()
  targeting_grnas_df <- object@grna_group_data_frame |>
    dplyr::filter(grna_group != "non-targeting")
  n_targeting_grna_groups <- length(unique(targeting_grnas_df$grna_group))
  n_targeting_grnas <- nrow(targeting_grnas_df)

  n_covariates <- ncol(object@covariate_data_frame)
  covariates <- paste0(sort(colnames(object@covariate_data_frame)), collapse = ", ")
  moi <- ifelse(object@low_moi, "Low", "High")
  cat(paste0("An object of class ", crayon::blue("sceptre_object"), ".\n\nAttributes of the data:\n\t\U2022 ",
             crayon::blue(n_cells), " cells\n\t\U2022 ",
             crayon::blue(n_responses), " responses\n\t\U2022 ",
             crayon::blue(moi), " multiplicity-of-infection \n\t\U2022 ",
             crayon::blue(n_nt_grnas), " non-targeting gRNAs \n\t\U2022 ",
             crayon::blue(n_targeting_grnas), " targeting gRNAs (distributed across ", crayon::blue(n_targeting_grna_groups), " gRNA groups) \n\t\U2022 ",
             crayon::blue(n_covariates), " covariates (", covariates, ")"))
})

# print method for sceptre class
setMethod("print", signature = signature("sceptre_object"), function(x, ...) {
  args <- list(...)
  print_full_output <- !is.null(args[["full_output"]]) && args[["full_output"]]
  show(x)
  get_mark <- function(bool) ifelse(bool, crayon::green("\u2713"), crayon::red("\u2717"))
  cat(paste0("\n\nAnalysis status:\n",
             "\t", get_mark(TRUE), " create_sceptre_object()\n",
             "\t", get_mark(x@analysis_prepared), " prepare_analysis()\n",
             "\t", get_mark(x@calibration_check_run), " run_calibration_check()\n",
             "\t", get_mark(x@power_check_run), " run_power_check()\n",
             "\t", get_mark(x@discovery_analysis_run), " run_discovery_analysis()"))
  n_discovery_pairs <- nrow(x@discovery_pairs)
  disc_pair_qc_performed <- !is.na(x@n_ok_discovery_pairs)
  n_pc_pairs <- nrow(x@positive_control_pairs)
  pc_pair_qc_performed <- !is.na(x@n_ok_positive_control_pairs)
  cat(paste0("\n\nUser-specified analysis parameters: \n",
             "\t\U2022 Discovery pairs:", if (n_discovery_pairs == 0) {" not specified"} else {paste0(" data frame with ", crayon::blue(n_discovery_pairs), " pairs",
                                                                                                      if (disc_pair_qc_performed) paste0(" (", crayon::blue(x@n_ok_discovery_pairs), " after pairwise QC)") else NULL)},
             "\n\t\U2022 Positive control pairs:", if (n_pc_pairs == 0) {" not specified"} else {paste0(" data frame with ", crayon::blue(n_pc_pairs), " pairs",
                                                                                                        if (pc_pair_qc_performed) paste0(" (", crayon::blue(x@n_ok_positive_control_pairs), " after pairwise QC)") else NULL)},
             "\n\t\U2022 Side: ", if (length(x@side_code) == 0L) "not specified" else crayon::blue(c("left", "both", "right")[x@side_code + 2L]),
             "\n\t\U2022 N nonzero treatment cells threshold: ", if (length(x@n_nonzero_trt_thresh) == 0L) "not specified" else crayon::blue(x@n_nonzero_trt_thresh),
             "\n\t\U2022 N nonzero control cells threshold: ", if (length(x@n_nonzero_cntrl_thresh) == 0L) "not specified" else crayon::blue(x@n_nonzero_cntrl_thresh),
             if (x@low_moi) NULL else {paste0("\n\t\U2022 gRNA assignment threshold: ", if (length(x@grna_assign_threshold) == 0L) "not specified" else crayon::blue(x@grna_assign_threshold))},
             if (!x@low_moi) NULL else {paste0("\n\t\U2022 Control group: ", if (length(x@control_group_complement) == 0L) "not specified" else crayon::blue(ifelse(x@control_group_complement, "complement set", "non-targeting cells")))},
             "\n\t\U2022 Formula object: ", if (length(x@formula_object) == 0L) "not specified" else crayon::blue(as.character(x@formula_object)[2]),
             if (!print_full_output) NULL else {
                paste0(
                  "\n\t\U2022 Resampling mechanism: ", if (length(x@run_permutations) == 0L) "not specified" else crayon::blue(ifelse(x@run_permutations, "permutations", "conditional resampling")),
                  "\n\t\U2022 Fit parametric curve: ", if (length(x@fit_parametric_curve) == 0L) "not specified" else crayon::blue(x@fit_parametric_curve),
                  "\n\t\U2022 B1: ", if (length(x@B1) == 0L) "not specified" else crayon::blue(x@B1), ", ",
                  "B2: ", if (length(x@B2) == 0L) "not specified" else crayon::blue(x@B2), ", ",
                  "B3: ", if (length(x@B3) == 0L) "not specified" else crayon::blue(x@B3),
                  "\n\t\U2022 Calibration check N pairs: ", if (length(x@n_calibration_pairs) == 0L) "not specified" else { if (is.na(x@n_calibration_pairs)) crayon::blue("match discovery pairs") else crayon::blue(x@n_calibration_pairs)},
                  "\n\t\U2022 Calibration check group size: ", if (length(x@calibration_group_size) == 0L) "not specified" else { if (is.na(x@calibration_group_size)) crayon::blue("match discovery pairs") else crayon::blue(x@calibration_group_size)}
                )
             }
  ))
})


# plot function for sceptre object
setMethod("plot", signature = signature("sceptre_object"), function(x) {
  last_function_called <- x@last_function_called
  if (last_function_called == "create_sceptre_object") {
    p <- plot_covariates(x)
  } else if (last_function_called == "prepare_analysis") {
    p <- plot_prepare_analysis(x)
  } else if (last_function_called == "run_calibration_check") {
    p <- plot_calibration_result(x)
  } else if (last_function_called == "run_discovery_analysis") {
    p <- plot_discovery_result(x)
  }
  return(p)
})

#' Create a `sceptre` object
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
#' sceptre_object <- prepare_analysis(
#' sceptre_object = sceptre_object,
#' discovery_pairs = "all")
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
#' sceptre_object <- prepare_analysis(
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
create_sceptre_object <- function(response_matrix, grna_matrix, grna_group_data_frame, moi, extra_covariates = NULL) {
  # 0. initialize output
  out <- new("sceptre_object")

  # 1. perform initial check
  check_create_sceptre_object_inputs(response_matrix, grna_matrix,
                                     grna_group_data_frame, moi, extra_covariates) |> invisible()

  # 2. compute the covariates
  covariate_data_frame <- auto_compute_cell_covariates(response_matrix, grna_matrix, extra_covariates)

  # 3. make the response matrix row accessible
  response_matrix <- set_matrix_accessibility(response_matrix, make_row_accessible = TRUE)

  # 4. update fields in output object and return
  # data fields
  out@response_matrix <- response_matrix
  out@grna_matrix <- grna_matrix
  out@covariate_data_frame <- covariate_data_frame
  out@grna_group_data_frame <- grna_group_data_frame
  out@low_moi <- (moi == "low")
  if (!is.null(extra_covariates)) {
    out@user_specified_covariates <- colnames(extra_covariates)
  }

  # cached fields
  out@calibration_check_run <- FALSE
  out@power_check_run <- FALSE
  out@discovery_analysis_run <- FALSE
  out@analysis_prepared <- FALSE
  out@last_function_called <- "create_sceptre_object"

  return(out)
}


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

  # update cached fields (four fields cached: response precomputations, grna assignments, negative control pairs, )
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
