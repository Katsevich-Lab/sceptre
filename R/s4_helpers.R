# show method
setMethod("show", signature = signature("sceptre_object"), function(object) {
  # 0. determine the functions that have been run
  funct_run_vect <- object@functs_called

  # 1. obtain the basic information
  n_cells <- ncol(object@response_matrix)
  n_responses <- nrow(object@response_matrix)
  n_nt_grnas <- object@grna_target_data_frame |>
    dplyr::filter(grna_target == "non-targeting") |> nrow()
  targeting_grnas_df <- object@grna_target_data_frame |>
    dplyr::filter(grna_target != "non-targeting")
  n_targeting_grna_targets <- length(unique(targeting_grnas_df$grna_target))
  n_targeting_grnas <- nrow(targeting_grnas_df)
  n_covariates <- ncol(object@covariate_data_frame)
  covariates <- paste0(sort(colnames(object@covariate_data_frame)), collapse = ", ")
  moi <- ifelse(object@low_moi, "Low", "High")
  cat(paste0("An object of class ", crayon::blue("sceptre_object"), ".\n\nAttributes of the data:\n\t\U2022 ",
             crayon::blue(n_cells), " cells", if (funct_run_vect["run_qc"]) {
               paste0(" (", crayon::blue(length(object@cells_in_use)), " after cellwise QC)")
             } else NULL, "\n\t\U2022 ",
             crayon::blue(n_responses), " responses\n\t\U2022 ",
             crayon::blue(moi), " multiplicity-of-infection \n\t\U2022 ",
             crayon::blue(n_targeting_grnas), " targeting gRNAs (distributed across ", crayon::blue(n_targeting_grna_targets), " targets) \n\t\U2022 ",
             crayon::blue(n_nt_grnas), " non-targeting gRNAs \n\t\U2022 ",
             crayon::blue(n_covariates), " covariates (", covariates, ")"))
})


#' Print
#'
#' `print()` prints information about the dataset and the status of the analysis to the console.
#'
#' @param x a `sceptre_object`
#' @return NULL
#' @export
setMethod("print", signature = signature("sceptre_object"), function(x) {
  show(x)

  # 1. print analysis status
  funct_run_vect <- x@functs_called
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
  cat(paste0("\nAnalysis parameters: \n",
             "\t\U2022 Discovery pairs:", if (n_discovery_pairs == 0) {" not specified"} else {paste0(" data frame with ", crayon::blue(n_discovery_pairs), " pairs",
                                                                                                      if (funct_run_vect["run_qc"]) paste0(" (", crayon::blue(x@n_ok_discovery_pairs), " after pairwise QC)") else NULL)},
             "\n\t\U2022 Positive control pairs:", if (n_pc_pairs == 0) {" not specified"} else {paste0(" data frame with ", crayon::blue(n_pc_pairs), " pairs",
                                                                                                        if (funct_run_vect["run_qc"]) paste0(" (", crayon::blue(x@n_ok_positive_control_pairs), " after pairwise QC)") else NULL)},
             "\n\t\U2022 Sidedness of test: ", if (length(x@side_code) == 0L) "not specified" else crayon::blue(c("left", "both", "right")[x@side_code + 2L]),
             "\n\t\U2022 N nonzero treatment cells threshold: ", if (length(x@n_nonzero_trt_thresh) == 0L) "not specified" else crayon::blue(x@n_nonzero_trt_thresh),
             "\n\t\U2022 N nonzero control cells threshold: ", if (length(x@n_nonzero_cntrl_thresh) == 0L) "not specified" else crayon::blue(x@n_nonzero_cntrl_thresh),
             if (!x@low_moi) NULL else {paste0("\n\t\U2022 Control group: ", if (length(x@control_group_complement) == 0L) "not specified" else crayon::blue(ifelse(x@control_group_complement, "complement set", "non-targeting cells")))},
             "\n\t\U2022 Resampling mechanism: ", if (length(x@run_permutations) == 0L) "not specified" else crayon::blue(ifelse(x@run_permutations, "permutations", "conditional resampling")),
             "\n\t\U2022 gRNA integration strategy: ", if (length(x@grna_integration_strategy) == 0L) "not specified" else crayon::blue(x@grna_integration_strategy),
             "\n\t\U2022 Formula object: ", if (length(x@formula_object) == 0L) "not specified" else crayon::blue(as.character(x@formula_object)[2])
  ))

  # 3. print the gRNA-to-cell assignment information
  grna_assignment_run <- funct_run_vect[["assign_grnas"]]
  if (grna_assignment_run) {
    mean_cells_per_grna <- sapply(x@initial_grna_assignment_list, length) |> mean()
    mean_cells_per_grna_group <- sapply(x@grna_assignments_raw$grna_group_idxs, length) |> mean()
    cat(paste0("\n\ngRNA-to-cell assignment information:",
               "\n\t\U2022 Assignment method: ", crayon::blue(x@grna_assignment_method),
               "\n\t\U2022 Mean N cells per gRNA: ", crayon::blue(mean_cells_per_grna |> round(2)),
               "\n\t\U2022 Mean N gRNAs per cell (MOI): ", if (x@grna_assignment_method == "maximum") "not computed when using \"maximum\" assignment method" else crayon::blue(x@grnas_per_cell |> mean() |> round(2))))
  }

  # 4. print the results summary
  calib_check_run <- funct_run_vect[["run_calibration_check"]]
  discovery_analysis_run <- funct_run_vect[["run_discovery_analysis"]]
  power_check_run <- funct_run_vect[["run_power_check"]]
  if (calib_check_run || discovery_analysis_run) cat("\n\nSummary of results:")
  if (calib_check_run) {
    n_false_rejections <- sum(x@calibration_result$significant)
    mean_lfc <- signif(mean(x@calibration_result$log_2_fold_change), 2)
    n_calib_pairs <- nrow(x@calibration_result)
    cat(paste0("\n\t\U2022 N", crayon::red(" negative control "), "pairs called as significant: ", crayon::blue(paste0(n_false_rejections, "/", n_calib_pairs)),
               "\n\t\U2022 Mean log-2 FC for", crayon::red(" negative control "), "pairs: ", crayon::blue(mean_lfc)))
  }
  if (power_check_run) {
    median_pc_p_value <- stats::median(x@power_result$p_value, na.rm = TRUE)
    cat(paste0("\n\t\U2022 Median", crayon::green(" positive control "), "p-value: ", crayon::blue(signif(median_pc_p_value, 2))))
  }
  if (discovery_analysis_run) {
    n_discoveries <- sum(x@discovery_result$significant, na.rm = TRUE)
    n_discovery_pairs <- x@n_ok_discovery_pairs
    cat(paste0("\n\t\U2022 N", crayon::yellow(" discovery pairs "), "called as significant: ", crayon::blue(paste0(n_discoveries, "/", n_discovery_pairs))))
  }
})


#' Plot
#'
#' `plot()` creates a plot depicting the current state of a `sceptre_object`.
#'
#' `plot()` is "generic" in the sense that it dispatches a specific plotting function based on the pipeline function that was most recently called on the `sceptre_object`. For example, if `run_assign_grnas()` is the most recently called pipeline function, then `plot()` dispatches `plot_run_assign_grnas()`. Similarly, if `run_power_check()` is the most recently called pipeline function, then `plot()` dispatches `plot_run_power_check()`, and so on. Users can pass arguments to the function dispatched by `plot()` as named arguments to `plot()`.
#'
#' @param x a `sceptre_object`
#' @param y ignored argument
#' @param ... arguments passed to the plotting function dispatched by `plot()`
#' @return a single \code{cowplot} object containing the combined panels (if \code{return_indiv_plots} is set to \code{TRUE}) or a list of the individual panels (if \code{return_indiv_plots} is set to \code{FALSE})
#'
#' @export
#' @examples
#' # A full example can be found at ?sceptre
setMethod("plot", signature = signature("sceptre_object"), function(x, y, ...) {
  args <- list(...)
  args[["sceptre_object"]] <- x

  last_function_called <- x@last_function_called
  if (last_function_called %in% c("import_data", "set_analysis_parameters")) {
    stop("There is no generic plot function configured for the ", last_function_called, "() step.")
  }

  funct_to_call <- paste0("plot_", last_function_called)
  p <- do.call(what = funct_to_call, args = args)
  return(p)
})


perform_status_check_and_update <- function(sceptre_object, curr_funct) {
  rank_vector <- c(import_data = 1L, set_analysis_parameters = 2L, assign_grnas = 3L,
                   run_qc = 4L, run_power_check = 5L, run_discovery_analysis = 5L, run_calibration_check = 5L)
  functs_called <- sceptre_object@functs_called
  curr_rank <- rank_vector[[curr_funct]]
  direct_upstream_functs <- names(rank_vector)[rank_vector == curr_rank - 1L]
  downstream_functs <- names(rank_vector)[rank_vector > curr_rank]
  # verify that the direct upstream funct has been called
  if (!all(functs_called[direct_upstream_functs])) {
    stop(paste0("The function", if (length(direct_upstream_functs) >= 2L) "s " else " ",
                paste0("`", direct_upstream_functs, "()`", collapse = ", "),
                " must be called before `", curr_funct, "()` is called."))
  }
  # set current funct to true
  functs_called[curr_funct] <- TRUE
  # set all downstream functions to false
  functs_called[downstream_functs] <- FALSE
  # update sceptre_object and return
  sceptre_object@functs_called <- functs_called
  sceptre_object@last_function_called <- curr_funct
  return(sceptre_object)
}


load_row <- function(mat, id) {
  if (methods::is(mat, "odm")) {
    mat[id,]
  } else if (methods::is(mat, "dgRMatrix")) {
    load_csr_row(j = mat@j, p = mat@p, x = mat@x, row_idx = which(id == rownames(mat)), n_cells = ncol(mat))
  }
}
