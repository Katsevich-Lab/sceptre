# response matrix class union
#' @import Matrix
setClassUnion("response_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix", "dgTMatrix"))
setClassUnion("grna_matrix_class", c("matrix", "dgCMatrix", "dgRMatrix", "dgTMatrix", "lgCMatrix", "lgRMatrix", "lgTMatrix"))


# sceptre object class
setClass("sceptre_object",
         slots = list(
           # raw data
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
           grna_assignment_method = "character",
           grna_assignment_hyperparameters = "list",
           multiple_testing_alpha = "numeric",
           multiple_testing_method = "character",

           # computed objects
           M_matrix = "matrix",
           n_nonzero_tot_vector = "integer",
           n_ok_discovery_pairs = "integer",
           n_ok_positive_control_pairs = "integer",
           discovery_pairs_with_info = "data.frame",
           positive_control_pairs_with_info = "data.frame",
           negative_control_pairs = "data.frame",
           initial_grna_assignment_list = "list",
           grna_assignments_raw = "list",
           grna_assignments = "list",
           cells_w_multiple_grnas = "integer",
           cells_per_grna = "integer",
           grnas_per_cell = "integer",
           cells_per_targeting_grna_group = "integer",
           cells_in_use = "integer",
           calibration_group_size = "integer",
           n_calibration_pairs = "integer",

           # cached objects
           response_precomputations = "list",

           # flags
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
             crayon::blue(n_cells), " cells", if (length(object@cells_in_use) >= 1) {
               paste0(" (", crayon::blue(length(object@cells_in_use)), " after cellwise QC)")
             } else NULL, "\n\t\U2022 ",
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

  # 0. determine the functions that have been run
  function_rank_vector <- get_function_rank_vector()
  current_function <- x@last_function_called
  curr_rank <- function_rank_vector[[current_function]]
  if (curr_rank <= 4) {
    completed_functs <- names(function_rank_vector)[function_rank_vector <= curr_rank]
  } else {
    completed_functs <- names(function_rank_vector)[function_rank_vector <= 4]
    if (nrow(x@calibration_result) >= 1L) completed_functs <- c(completed_functs, "run_calibration_check")
    if (nrow(x@power_result) >= 1L) completed_functs <- c(completed_functs, "run_power_check")
    if (nrow(x@discovery_result) >= 1L) completed_functs <- c(completed_functs, "run_discovery_analysis")
  }
  funct_run_vect <- sapply(names(function_rank_vector), function(funct) funct %in% completed_functs)

  # 1. print analysis status
  get_mark <- function(bool) ifelse(bool, crayon::green("\u2713"), crayon::red("\u2717"))
  cat("\n\nAnalysis status:\n")
  for (i in seq_along(funct_run_vect)) {
    cat(paste0("\t", get_mark(funct_run_vect[i]), " ", names(funct_run_vect)[i], "()\n"))
  }

  # 2. print analysis parameters
  n_discovery_pairs <- nrow(x@discovery_pairs)
  disc_pair_qc_performed <- length(x@n_ok_discovery_pairs) >= 1
  n_pc_pairs <- nrow(x@positive_control_pairs)
  pc_pair_qc_performed <- length(x@n_ok_positive_control_pairs) >= 1
  cat(paste0("\nUser-specified analysis parameters: \n",
             "\t\U2022 Discovery pairs:", if (n_discovery_pairs == 0) {" not specified"} else {paste0(" data frame with ", crayon::blue(n_discovery_pairs), " pairs",
                                                                                                      if (disc_pair_qc_performed) paste0(" (", crayon::blue(x@n_ok_discovery_pairs), " after pairwise QC)") else NULL)},
             "\n\t\U2022 Positive control pairs:", if (n_pc_pairs == 0) {" not specified"} else {paste0(" data frame with ", crayon::blue(n_pc_pairs), " pairs",
                                                                                                        if (pc_pair_qc_performed) paste0(" (", crayon::blue(x@n_ok_positive_control_pairs), " after pairwise QC)") else NULL)},
             "\n\t\U2022 Sidedness of test: ", if (length(x@side_code) == 0L) "not specified" else crayon::blue(c("left", "both", "right")[x@side_code + 2L]),
             "\n\t\U2022 N nonzero treatment cells threshold: ", if (length(x@n_nonzero_trt_thresh) == 0L) "not specified" else crayon::blue(x@n_nonzero_trt_thresh),
             "\n\t\U2022 N nonzero control cells threshold: ", if (length(x@n_nonzero_cntrl_thresh) == 0L) "not specified" else crayon::blue(x@n_nonzero_cntrl_thresh),
             if (!x@low_moi) NULL else {paste0("\n\t\U2022 Control group: ", if (length(x@control_group_complement) == 0L) "not specified" else crayon::blue(ifelse(x@control_group_complement, "complement set", "non-targeting cells")))},
             "\n\t\U2022 Resampling mechanism: ", if (length(x@run_permutations) == 0L) "not specified" else crayon::blue(ifelse(x@run_permutations, "permutations", "conditional resampling")),
             if (!print_full_output) NULL else {
               paste0(
                 "\n\t\U2022 Fit parametric curve: ", if (length(x@fit_parametric_curve) == 0L) "not specified" else crayon::blue(x@fit_parametric_curve),
                 "\n\t\U2022 B1: ", if (length(x@B1) == 0L) "not specified" else crayon::blue(x@B1), ", ",
                 "B2: ", if (length(x@B2) == 0L) "not specified" else crayon::blue(x@B2), ", ",
                 "B3: ", if (length(x@B3) == 0L) "not specified" else crayon::blue(x@B3))
             },
             "\n\t\U2022 Formula object: ", if (length(x@formula_object) == 0L) "not specified" else crayon::blue(as.character(x@formula_object)[2])
  ))

  # 3. print the gRNA-to-cell assignment information
  grna_assignment_run <- funct_run_vect[["assign_grnas"]]
  if (grna_assignment_run) {
    cat(paste0("\n\ngRNA-to-cell assignment information:",
               "\n\t\U2022 Assignment method: ", crayon::blue(x@grna_assignment_method),
               "\n\t\U2022 Mean N cells per gRNA: ", crayon::blue(x@cells_per_grna |> mean() |> round(2)),
               "\n\t\U2022 Mean N cells per targeting gRNA group: ", crayon::blue(x@cells_per_targeting_grna_group |> mean() |> round(2)),
               "\n\t\U2022 Mean N gRNAs per cell: ", crayon::blue(x@grnas_per_cell |> mean() |> round(2))))
  }

  # 4. print the results summary
  calib_check_run <- funct_run_vect[["run_calibration_check"]]
  discovery_analysis_run <- funct_run_vect[["run_discovery_analysis"]]
  power_check_run <- funct_run_vect[["run_power_check"]]
  if (calib_check_run || discovery_analysis_run) cat("\n\nSummary of results:")
  if (calib_check_run) {
    n_false_rejections <- sum(x@calibration_result$reject)
    mean_lfc <- signif(mean(x@calibration_result$log_2_fold_change), 2)
    n_calib_pairs <- nrow(x@calibration_result)
    cat(paste0("\n\t\U2022 N", crayon::red(" negative control "), "pairs rejected: ", crayon::blue(paste0(n_false_rejections, "/", n_calib_pairs)),
               "\n\t\U2022 Mean log-2 FC for", crayon::red(" negative control "), "pairs: ", crayon::blue(mean_lfc)))
  }
  if (power_check_run) {
    median_pc_p_value <- median(x@power_result$p_value, na.rm = TRUE)
    cat(paste0("\n\t\U2022 Median", crayon::green(" positive control "), "p-value: ", crayon::blue(signif(median_pc_p_value, 2))))
  }
  if (discovery_analysis_run) {
    n_discoveries <- sum(x@discovery_result$reject, na.rm = TRUE)
    n_discovery_pairs <- x@n_ok_discovery_pairs
    cat(paste0("\n\t\U2022 N", crayon::yellow(" discovery pairs "), "rejected: ", crayon::blue(paste0(n_discoveries, "/", n_discovery_pairs))))
  }
})


# plot function for sceptre object
setMethod("plot", signature = signature("sceptre_object"), function(x, ...) {
  args <- list(...)
  args[["sceptre_object"]] <- x

  last_function_called <- x@last_function_called
  if (last_function_called %in% c("import_data", "set_analysis_parameters", "run_power_check")) {
    stop("There is no generic plot function configured for the ", last_function_called, "() step.")
  }

  funct_to_call <- paste0("plot_", last_function_called)
  p <- do.call(what = funct_to_call, args = args)
  return(p)
})


reset_results <- function(sceptre_object) {
  sceptre_object@calibration_result <- data.frame()
  sceptre_object@discovery_result <- data.frame()
  sceptre_object@power_result <- data.frame()
  sceptre_object@n_ok_discovery_pairs <- integer()
  sceptre_object@n_ok_positive_control_pairs <- integer()
  return(sceptre_object)
}
