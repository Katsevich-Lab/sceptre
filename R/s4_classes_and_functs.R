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

                      # cached objects
                      response_precomputations = "list",
                      grna_assignments = "list",
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
  covariates <- paste0(colnames(object@covariate_data_frame), collapse = ", ")
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
setMethod("print", signature = signature("sceptre_object"), function(x) {
  show(x)
  get_mark <- function(bool) ifelse(bool, crayon::green("\u2713"), crayon::red("\u2717"))
  cat(paste0("\n\nAnalysis status:\n",
             "\t", get_mark(x@analysis_prepared), " prepare_analysis()\n",
             "\t", get_mark(x@calibration_check_run), " run_calibration_check()\n",
             "\t", get_mark(x@power_check_run), " run_power_check()\n",
             "\t", get_mark(x@discovery_analysis_run), " run_discovery_analysis()"))

  n_discovery_pairs <- nrow(x@discovery_pairs)
  n_pc_pairs <- nrow(x@positive_control_pairs)
  cat(paste0("\n\nUser-specified analysis parameters: \n",
             "\t\U2022 Discovery pairs:", if (n_discovery_pairs == 0) {" not provided"} else {paste0(" data frame with ", crayon::blue(n_discovery_pairs), " pairs")},
             "\n\t\U2022 Positive control pairs:", if (n_pc_pairs == 0) {" not provided"} else {paste0(" data frame with ", crayon::blue(n_pc_pairs), " pairs")},
             "\n\t\U2022 Formula object: ", if (length(x@formula_object) == 0L) "not provided" else crayon::blue(as.character(x@formula_object)[2]),
             "\n\t\U2022 Side: ", if (length(x@side_code) == 0L) "not provided" else crayon::blue(c("left", "both", "right")[x@side_code + 2L]),
             "\n\t\U2022 Control group: ", if (length(x@control_group_complement) == 0L) "not provided" else crayon::blue(ifelse(x@control_group_complement, "complement set", "non-targeting cells")),
             "\n\t\U2022 Resampling mechanism: ", if (length(x@run_permutations) == 0L) "not provided" else crayon::blue(ifelse(x@run_permutations, "permutations", "conditional resampling")),
             "\n\t\U2022 N nonzero treatment cells QC threshold: ", if (length(x@n_nonzero_trt_thresh) == 0L) "not provided" else crayon::blue(x@n_nonzero_trt_thresh),
             "\n\t\U2022 N nonzero control cells QC threshold: ", if (length(x@n_nonzero_cntrl_thresh) == 0L) "not provided" else crayon::blue(x@n_nonzero_cntrl_thresh),
             if (x@low_moi) NULL else {paste0("\n\t\U2022 gRNA assignment threshold: ", if (length(x@grna_assign_threshold) == 0L) "not provided" else crayon::blue(x@grna_assign_threshold))},
             "\n\t\U2022 Fit parametric curve: ", if (length(x@fit_parametric_curve) == 0L) "not provided" else crayon::blue(x@fit_parametric_curve),
             "\n\t\U2022 B1: ", if (length(x@B1) == 0L) "not provided" else crayon::blue(x@B1), ", ",
             "B2: ", if (length(x@B2) == 0L) "not provided" else crayon::blue(x@B2), ", ",
             "B3: ", if (length(x@B3) == 0L) "not provided" else crayon::blue(x@B3)
  ))
})

# plot function for sceptre object
setMethod("plot", signature = signature("sceptre_object"), function(x) {
  last_function_called <- x@last_function_called
  if (last_function_called == "run_calibration_check") {
    p <- plot_calibration_result(x)
  }
  return(p)
})

#' Create a `sceptre` object
#'
#' @inheritParams run_sceptre
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
#' # 2. obtain the formula object and pairs to analyze
#' formula_object <- formula(~log(response_n_umis) + log(response_n_nonzero) + bio_rep + p_mito)
#' discovery_pairs <- generate_all_pairs(response_matrix_lowmoi,
#' grna_group_data_frame_lowmoi)
#'
#' # 3. prepare the analysis
#' sceptre_object <- prepare_analysis(
#' sceptre_object = sceptre_object,
#' formula_object = formula_object,
#' discovery_pairs = discovery_pairs)
#'
#' # 4. run the calibration check
#' sceptre_object <- run_calibration_check(sceptre_object)
#' plot(sceptre_object)
#'
#' # 5. (optional) run the power check
#' sceptre_object <- run_power_check(sceptre_object)
#' plot(sceptre_object)
#'
#' # 6. run discovery analysis
#' sceptre_object <- run_discovery_analysis(sceptre_object)
#' plot(sceptre_object)
#'
#' ##################
#' # High MOI example
#' ##################
#' # 0. obtain the data required for a single-cell screen analysis
#' data(response_matrix_highmoi_experimental)
#' data(grna_matrix_highmoi_experimental)
#' data(covariate_data_frame_highmoi_experimental)
#' data(grna_group_data_frame_highmoi_experimental)
#'
#' # 1. create the sceptre object
#' sceptre_object <- create_sceptre_object(
#' response_matrix = response_matrix_highmoi_experimental,
#' grna_matrix = grna_matrix_highmoi_experimental,
#' covariate_data_frame = covariate_data_frame_highmoi_experimental,
#' grna_group_data_frame = grna_group_data_frame_highmoi_experimental,
#' moi = "high")
#'
#' # 2. set the formula, positive control pairs, and discovery pairs
#' formula_object <- formula(~log(grna_n_umis) + log(grna_n_nonzero) + log(gene_n_umis) + log(gene_n_nonzero) + batch + p_mito)
#' data(discovery_pairs_highmoi_experimental)
#' data(pc_pairs_highmoi_experimental)
#'
#' # 3. prepare the analysis
#' sceptre_object <- prepare_analysis(
#' sceptre_object = sceptre_object,
#' formula_object = formula_object,
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
  # data fields
  out@response_matrix <- response_matrix
  out@grna_matrix <- grna_matrix
  out@covariate_data_frame <- covariate_data_frame
  out@grna_group_data_frame <- grna_group_data_frame
  out@low_moi <- (moi == "low")

  # cached fields
  out@response_precomputations <- list()
  out@calibration_check_run <- FALSE
  out@power_check_run <- FALSE
  out@discovery_analysis_run <- FALSE
  out@analysis_prepared <- FALSE
  out@calibration_result <- data.frame()
  out@discovery_result <- data.frame()
  out@last_function_called <- "create_sceptre_object"

  return(out)
}


prepare_analysis <- function(sceptre_object,
                             formula_object,
                             discovery_pairs = data.frame(),
                             positive_control_pairs = data.frame(),
                             side = "both",
                             fit_parametric_curve = TRUE,
                             control_group = "default",
                             resampling_mechanism = "default",
                             n_nonzero_trt_thresh = 7L,
                             n_nonzero_cntrl_thresh = 7L,
                             grna_assign_threshold = 5L,
                             B1 = 499L, B2 = 4999L, B3 = 24999L) {
  # 1. check inputs
  check_prepare_analysis_inputs(response_matrix = sceptre_object@response_matrix,
                                grna_matrix = sceptre_object@grna_matrix,
                                covariate_data_frame = sceptre_object@covariate_data_frame,
                                grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                formula_object = formula_object,
                                response_grna_group_pairs_list = list(discovery_pairs, positive_control_pairs),
                                control_group = control_group,
                                resampling_mechanism = resampling_mechanism,
                                side = side, low_moi = sceptre_object@low_moi) |> invisible()

  # 2. check whether to update cached objects
  discard_response_precomputations <- !identical(sceptre_object@formula_object, formula_object)

  # 3. update analysis parameter fields of sceptre object
  if (!sceptre_object@low_moi) {
    control_group <- "complement"
    if (resampling_mechanism == "default") resampling_mechanism <- "crt"
  }
  if (sceptre_object@low_moi) {
    if (control_group == "default") control_group <- "nt_cells"
    if (resampling_mechanism == "default") resampling_mechanism <- "permutations"
  }
  sceptre_object@control_group_complement <- control_group == "complement"
  sceptre_object@run_permutations <- resampling_mechanism == "permutations"
  sceptre_object@side_code <- which(side == c("left", "both", "right")) - 2L
  sceptre_object@fit_parametric_curve <- fit_parametric_curve
  sceptre_object@B1 <- B1
  sceptre_object@B2 <- B2
  sceptre_object@B3 <- B3
  sceptre_object@n_nonzero_trt_thresh <- n_nonzero_trt_thresh
  sceptre_object@n_nonzero_cntrl_thresh <- n_nonzero_cntrl_thresh
  sceptre_object@discovery_pairs <- discovery_pairs
  sceptre_object@positive_control_pairs <- positive_control_pairs
  sceptre_object@grna_assign_threshold <- grna_assign_threshold
  sceptre_object@formula_object <- formula_object

  # 4. update cached fields
  if (discard_response_precomputations) sceptre_object@response_precomputations <- list()
  sceptre_object@calibration_check_run <- FALSE
  sceptre_object@power_check_run <- FALSE
  sceptre_object@discovery_analysis_run <- FALSE
  sceptre_object@analysis_prepared <- TRUE
  sceptre_object@calibration_result <- data.frame()
  sceptre_object@power_result <- data.frame()
  sceptre_object@discovery_result <- data.frame()
  sceptre_object@last_function_called <- "prepare_analysis"

  # 5. create the covariate matrix, check it for correctness, and insert it into object
  sceptre_object@covariate_matrix <- convert_covariate_df_to_design_matrix(
    sceptre_object@covariate_data_frame, formula_object)

  # 6. assign grnas
  grna_assignments <- assign_grnas_to_cells(grna_matrix = sceptre_object@grna_matrix,
                                            grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                            grna_assign_threshold = sceptre_object@grna_assign_threshold,
                                            low_moi = sceptre_object@low_moi,
                                            control_group_complement = sceptre_object@control_group_complement)
  sceptre_object@grna_assignments <- grna_assignments

  return(sceptre_object)
}


run_calibration_check <- function(sceptre_object, calibration_group_size = NULL,
                                  n_calibration_pairs = NULL, output_amount = 1, print_progress = TRUE) {
  # 1. do argument check
  check_calibration_check_inputs(analysis_prepared = sceptre_object@analysis_prepared,
                                 grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                 control_group_complement = sceptre_object@control_group_complement) |> invisible()

  # 2. get grna assignments
  grna_assignments <- sceptre_object@grna_assignments

  # 3. construct the negative control pairs
  cat("Constructing negative control pairs.")
  if (is.null(calibration_group_size)) calibration_group_size <- compute_calibration_group_size(sceptre_object@grna_group_data_frame)
  response_grna_group_pairs <- construct_negative_control_pairs(n_calibration_pairs = n_calibration_pairs,
                                                                calibration_group_size = calibration_group_size,
                                                                grna_assignments = grna_assignments,
                                                                response_matrix = sceptre_object@response_matrix,
                                                                n_nonzero_trt_thresh = sceptre_object@n_nonzero_trt_thresh,
                                                                n_nonzero_cntrl_thresh = sceptre_object@n_nonzero_cntrl_thresh,
                                                                grna_group_data_frame = sceptre_object@grna_group_data_frame,
                                                                response_grna_group_pairs = sceptre_object@discovery_pairs,
                                                                control_group_complement = sceptre_object@control_group_complement,
                                                                low_moi = sceptre_object@low_moi)
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
                                   print_progress = print_progress)
  } else {
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
  sceptre_object@calibration_result <- out$ret
  sceptre_object@response_precomputations <- out$response_precomputations
  sceptre_object@calibration_check_run <- TRUE
  sceptre_object@last_function_called <- "run_calibration_check"
  return(sceptre_object)
}


run_power_check <- function(sceptre_object, output_amount = 1, print_progress = TRUE) {
  # 1. do argument check
  check_discovery_analysis_inputs(response_grna_group_pairs = sceptre_object@positive_control_pairs,
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
                                   response_grna_group_pairs = sceptre_object@positive_control_pairs,
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
                                   print_progress = print_progress)
  } else {
    ret <- run_crt_in_memory_v2(response_matrix = sceptre_object@response_matrix,
                                grna_assignments = grna_assignments,
                                covariate_matrix = sceptre_object@covariate_matrix,
                                response_grna_group_pairs = sceptre_object@positive_control_pairs,
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

  # update fields of sceptre object
  sceptre_object@power_result <- out$ret
  sceptre_object@response_precomputations <- out$response_precomputations
  sceptre_object@power_check_run <- TRUE
  sceptre_object@last_function_called <- "run_power_check"
  return(sceptre_object)
}


run_discovery_analysis <- function(sceptre_object, output_amount = 1, print_progress = TRUE) {
  # 1. do argument check
  check_discovery_analysis_inputs(response_grna_group_pairs = sceptre_object@discovery_pairs,
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

  # 5. run the method
  if (sceptre_object@run_permutations) {
    out <- run_perm_test_in_memory(response_matrix = sceptre_object@response_matrix,
                                   grna_assignments = grna_assignments,
                                   covariate_matrix = sceptre_object@covariate_matrix,
                                   response_grna_group_pairs = sceptre_object@discovery_pairs,
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
                                   print_progress = print_progress)
  } else {
    ret <- run_crt_in_memory_v2(response_matrix = sceptre_object@response_matrix,
                                grna_assignments = grna_assignments,
                                covariate_matrix = sceptre_object@covariate_matrix,
                                response_grna_group_pairs = sceptre_object@discovery_pairs,
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

  # update fields of sceptre object
  sceptre_object@discovery_result <- out$ret
  sceptre_object@response_precomputations <- out$response_precomputations
  sceptre_object@discovery_analysis_run <- TRUE
  sceptre_object@last_function_called <- "run_discovery_analysis"
  return(sceptre_object)
}
