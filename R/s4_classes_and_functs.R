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

                      # analysis parameters
                      response_grna_group_pairs = "data.frame",
                      side_code = "integer",
                      fit_parametric_curve = "logical",
                      control_group_complement = "logical",
                      run_permutations = "logical",
                      n_nonzero_trt_thresh = "integer",
                      n_nonzero_cntrl_thresh = "integer",
                      B1 = "integer", B2 = "integer", B3 = "integer",
                      grna_assign_threshold = "integer",

                      # cached objects
                      response_precomputations = "list",
                      negative_control_pairs = "data.frame"))

#' Create a `sceptre` object
#'
#' @inheritParams run_sceptre
#' @return an object of class `sceptre_object`
#' @export
#'
#' @examples
#' # Low MOI example
#' data(response_matrix_lowmoi) # response-by-cell expression matrix
#' data(grna_matrix_lowmoi) # gRNA-by-cell expression matrix
#' data(covariate_data_frame_lowmoi) # cell-by-covariate data frame
#' data(grna_group_data_frame_lowmoi) # gRNA group information
#'
#' # 1. create the sceptre object
#' sceptre_object <- create_sceptre_object(
#' response_matrix = response_matrix_lowmoi,
#' grna_matrix = grna_matrix_lowmoi,
#' covariate_data_frame = covariate_data_frame_lowmoi,
#' grna_group_data_frame = grna_group_data_frame_lowmoi,
#' moi = "low")
#'
#' # 2. plan analysis; this includes setting the formula and obtaining the pairs to analyze
#' formula_object <- formula(~log(response_n_umis) + log(response_n_nonzero) + bio_rep + p_mito)
#' response_grna_group_pairs <- generate_all_pairs(response_matrix_lowmoi,
#' grna_group_data_frame_lowmoi)
#' sceptre_object <- plan_analysis(
#' sceptre_object = sceptre_object,
#' formula_object = formula_object,
#' response_grna_group_pairs = response_grna_group_pairs)
#'
#' # 3. run calibration check
#'
create_sceptre_object <- function(response_matrix, grna_matrix,
                                  covariate_data_frame, grna_group_data_frame, moi) {
  # 0. initialize output
  out <- new("sceptre_object")

  # 1. perform initial check
  check_create_sceptre_object_inputs(response_matrix, grna_matrix, covariate_data_frame,
                                     grna_group_data_frame, moi) |> invisible()

  # 2. make the response matrix row accessible
  response_matrix <- set_matrix_accessibility(response_matrix, make_row_accessible = TRUE)

  # 3. update fields in output object and return
  out@response_matrix <- response_matrix
  out@grna_matrix <- grna_matrix
  out@covariate_data_frame <- covariate_data_frame
  out@grna_group_data_frame <- grna_group_data_frame
  out@low_moi <- (moi == "low")
  return(out)
}


plan_analysis <- function(sceptre_object,
                          formula_object,
                          response_grna_group_pairs = NULL,
                          side = "both",
                          fit_parametric_curve = TRUE,
                          control_group = "default",
                          resampling_mechanism = "default",
                          n_nonzero_trt_thresh = 7L,
                          n_nonzero_cntrl_thresh = 7L,
                          grna_assign_threshold = 5L,
                          B1 = 499L, B2 = 4999L, B3 = 24999L) {
  # 1. check inputs
  check_plan_analysis_inputs(response_matrix = sceptre_object@response_matrix,
                             grna_matrix = sceptre_object@grna_matrix,
                             covariate_data_frame = sceptre_object@covariate_data_frame,
                             grna_group_data_frame = sceptre_object@grna_group_data_frame,
                             formula_object = formula_object,
                             response_grna_group_pairs = response_grna_group_pairs,
                             control_group = control_group,
                             resampling_mechanism = resampling_mechanism,
                             side = side, low_moi = sceptre_object@low_moi) |> invisible()

  # 2. update fields of sceptre object
  sceptre_object <- set_fields_plan_analysis(sceptre_object, fit_parametric_curve, sceptre_object@low_moi,
                                             control_group, resampling_mechanism, side, B1, B2, B3,
                                             n_nonzero_trt_thresh, n_nonzero_cntrl_thresh, response_grna_group_pairs,
                                             grna_assign_threshold)

  # 3. create the covariate matrix, check it for correctness, and insert it into object
  sceptre_object@covariate_matrix <- convert_covariate_df_to_design_matrix(
    sceptre_object@covariate_data_frame,
    formula_object)

  return(sceptre_object)
}


run_calibration_check <- function(sceptre_object, calibration_group_size = NULL,
                                  n_calibration_pairs = NULL, output_amount = 1, print_progress = TRUE) {
  # 1. perform argument check
  check_calibration_check_inputs(grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                 control_group_complement = sceptre_object@control_group_complement) |> invisible()

  # 2. assign gRNAs to cells
  grna_assignments <- assign_grnas_to_cells(grna_matrix = sceptre_object@grna_matrix,
                                            grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                            grna_assign_threshold = sceptre_object@grna_assign_threshold,
                                            low_moi = sceptre_object@low_moi,
                                            control_group_complement = sceptre_object@control_group_complement,
                                            calibration_check = TRUE)
  # 3. construct the negative control pairs
  cat("Constructing negative control pairs.")
  response_grna_group_pairs <- construct_negative_control_pairs(n_calibration_pairs = n_calibration_pairs,
                                                                calibration_group_size = calibration_group_size,
                                                                grna_assignments = grna_assignments,
                                                                response_matrix = sceptre_object@response_matrix,
                                                                n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                                                n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                                                grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                                                response_grna_group_pairs = sceptre_object@response_grna_group_pairs,
                                                                control_group_complement = sceptre_object@control_group_complement,
                                                                low_moi = sceptre_object@low_moi)
  cat(crayon::green(' \u2713\n'))

  # 7. generate the set of synthetic indicator idxs (if running permutations)
  if (run_permutations) {
    cat("Generating permutation resamples.")
    n_cells <- nrow(covariate_matrix)
    synthetic_idxs <- get_synthetic_permutation_idxs(grna_assignments, B1 + B2 + B3, calibration_check,
                                                     control_group_complement, calibration_group_size, n_cells)
    cat(crayon::green(' \u2713\n'))
  }
  gc() |> invisible()

}
