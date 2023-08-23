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

           # computed objects
           M_matrix = "matrix",
           n_nonzero_tot_vector = "integer",
           n_ok_discovery_pairs = "integer",
           n_ok_positive_control_pairs = "integer",
           discovery_pairs_with_info = "data.frame",
           positive_control_pairs_with_info = "data.frame",
           negative_control_pairs = "data.frame",
           grna_assignments = "list",
           grna_assignment_extra_info = "list",

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

  # first, print analysis status; determine the functions that have been run
  get_mark <- function(bool) ifelse(bool, crayon::green("\u2713"), crayon::red("\u2717"))
  ordered_functs <- names(get_function_rank_vector())
  curr_funct <- x@last_function_called
  functs_called <- ordered_functs[seq(1, which(curr_funct == ordered_functs))]
  cat("\n\nAnalysis status:\n")
  # special case: has run_power_check been called?
  if (curr_funct == "run_discovery_analysis" && nrow(sceptre_object@power_result) == 0L) {
    functs_called <- functs_called[functs_called != "run_power_check"]
  }
  functs_not_called <- ordered_functs[!(ordered_functs %in% functs_called)]
  for (funct_called in functs_called) {
    cat(paste0("\t", get_mark(TRUE), " ", funct_called, "()\n"))
  }
  for (funct_not_called in functs_not_called) {
    cat(paste0("\t", get_mark(FALSE), " ", funct_not_called, "()\n"))
  }

  # second, print analysis parameters
  n_discovery_pairs <- nrow(x@discovery_pairs)
  disc_pair_qc_performed <- length(x@n_ok_discovery_pairs) >= 1
  n_pc_pairs <- nrow(x@positive_control_pairs)
  pc_pair_qc_performed <- length(x@n_ok_positive_control_pairs) >= 1
  cat(paste0("\nUser-specified analysis parameters: \n",
             "\t\U2022 Discovery pairs:", if (n_discovery_pairs == 0) {" not specified"} else {paste0(" data frame with ", crayon::blue(n_discovery_pairs), " pairs",
                                                                                                      if (disc_pair_qc_performed) paste0(" (", crayon::blue(x@n_ok_discovery_pairs), " after pairwise QC)") else NULL)},
             "\n\t\U2022 Positive control pairs:", if (n_pc_pairs == 0) {" not specified"} else {paste0(" data frame with ", crayon::blue(n_pc_pairs), " pairs",
                                                                                                        if (pc_pair_qc_performed) paste0(" (", crayon::blue(x@n_ok_positive_control_pairs), " after pairwise QC)") else NULL)},
             "\n\t\U2022 Side: ", if (length(x@side_code) == 0L) "not specified" else crayon::blue(c("left", "both", "right")[x@side_code + 2L]),
             "\n\t\U2022 N nonzero treatment cells threshold: ", if (length(x@n_nonzero_trt_thresh) == 0L) "not specified" else crayon::blue(x@n_nonzero_trt_thresh),
             "\n\t\U2022 N nonzero control cells threshold: ", if (length(x@n_nonzero_cntrl_thresh) == 0L) "not specified" else crayon::blue(x@n_nonzero_cntrl_thresh),
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


reset_results <- function(sceptre_object) {
  sceptre_object@calibration_result <- data.frame()
  sceptre_object@discovery_result <- data.frame()
  sceptre_object@power_result <- data.frame()
  return(sceptre_object)
}
