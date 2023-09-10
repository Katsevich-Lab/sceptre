# show method
setMethod("show", signature = signature("sceptre_object"), function(object) {
  # 0. determine the functions that have been run
  funct_run_vect <- object@functs_called

  # 1. obtain the basic information
  n_cells <- ncol(object@response_matrix)
  n_responses <- nrow(object@response_matrix)
  n_nt_grnas <- object@grna_group_data_frame |>
    dplyr::filter(grna_group == "non-targeting") |> nrow()
  targeting_grnas_df <- object@grna_group_data_frame |>
    dplyr::filter(grna_group != "non-targeting")
  n_targeting_grna_groups <- length(unique(targeting_grnas_df$grna_group))
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
             crayon::blue(n_nt_grnas), " non-targeting gRNAs \n\t\U2022 ",
             crayon::blue(n_targeting_grnas), " targeting gRNAs (distributed across ", crayon::blue(n_targeting_grna_groups), " gRNA groups) \n\t\U2022 ",
             crayon::blue(n_covariates), " covariates (", covariates, ")"))
})


# print method
#' @export
setMethod("print", signature = signature("sceptre_object"), function(x, ...) {
  args <- list(...)
  print_full_output <- !is.null(args[["full_output"]]) && args[["full_output"]]
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
               "\n\t\U2022 Mean N gRNAs per cell (MOI): ", crayon::blue(x@grnas_per_cell |> mean() |> round(2))))
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
    median_pc_p_value <- median(x@power_result$p_value, na.rm = TRUE)
    cat(paste0("\n\t\U2022 Median", crayon::green(" positive control "), "p-value: ", crayon::blue(signif(median_pc_p_value, 2))))
  }
  if (discovery_analysis_run) {
    n_discoveries <- sum(x@discovery_result$significant, na.rm = TRUE)
    n_discovery_pairs <- x@n_ok_discovery_pairs
    cat(paste0("\n\t\U2022 N", crayon::yellow(" discovery pairs "), "called as significant: ", crayon::blue(paste0(n_discoveries, "/", n_discovery_pairs))))
  }
})


# plot method
#' @export
setMethod("plot", signature = signature("sceptre_object"), function(x, ...) {
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


# write outputs to directory
#' @export
write_outputs_to_directory <- function(sceptre_object, directory) {
  # 0. create directory
  if (!dir.exists(directory)) {
    dir.create(path = directory, recursive = TRUE)
  }

  # 1. create analysis_summary.txt file
  summary_file_fp <- paste0(directory, "/analysis_summary.txt")
  sink(file = summary_file_fp, append = FALSE)
  print(sceptre_object)
  sink(NULL)

  # 2. determine the functions that have been called
  plotting_params <- list(device = "png", scale = 1.1, width = 5, height = 4, dpi = 330)
  functs_run <- get_funct_run_vect(sceptre_object)
  functs_run_plots <- functs_run[names(functs_run) %in% c("assign_grnas", "run_qc", "run_calibration_check", "run_power_check", "run_discovery_analysis")]
  functs_run_plots <- c(functs_run_plots, "grna_count_distributions" = TRUE)
  for (funct in names(functs_run_plots)) {
    if (functs_run_plots[[funct]]) {
      plotting_params$plot <- do.call(what = paste0("plot_", funct), args = list(sceptre_object))
      plotting_params$filename <- paste0(directory, "/plot_", funct, ".png")
      do.call(what = ggplot2::ggsave, args = plotting_params)
    }
  }

  # 3. save results
  for (analysis in c("run_calibration_check", "run_power_check", "run_discovery_analysis")) {
    if (functs_run_plots[[analysis]]) {
      res <- get_result(sceptre_object = sceptre_object, analysis = analysis)
      saveRDS(object = res, file = paste0(directory, "/results_", analysis, ".rds"))
    }
  }

  return(NULL)
}


perform_status_check_and_update <- function(sceptre_object, curr_funct) {
  rank_vector <- c(import_data = 1L, set_analysis_parameters = 2L, assign_grnas = 2L,
                   run_qc = 3L, run_power_check = 4L, run_discovery_analysis = 4L, run_calibration_check = 4L)
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
