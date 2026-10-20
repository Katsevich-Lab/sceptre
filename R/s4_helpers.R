# Show the original matrix dimensions and which workflow steps have finished.
setMethod("show", signature = signature("sceptre_object"), function(object) {
    object <- migrate_sceptre_object(object)
    # 0. determine the functions that have been run
    funct_run_vect <- workflow_status(object)

    # 1. obtain the basic information
    n_cells <- ncol(get_response_matrix(object))
    n_responses <- nrow(get_response_matrix(object))
    grna_target_data_frame <- guide_group_metadata(object)
    guide_counts <- screen_guide_counts(grna_target_data_frame)
    n_covariates <- length(names(object@data@covariate_data_frame))
    covariates <- paste0(
        sort(names(object@data@covariate_data_frame)),
        collapse = ", "
    )
    moi <- ifelse(object@data@low_moi, "Low", "High")
    cat(paste0(
        "An object of class ",
        crayon::blue("sceptre_object"),
        ".\n\nAttributes of the data:\n\t\U2022 ",
        crayon::blue(n_cells),
        " cells",
        if (funct_run_vect["run_qc"]) {
            paste0(
                " (",
                crayon::blue(length(object@preparation@analysis_cells)),
                " after cellwise QC)"
            )
        } else {
            NULL
        },
        "\n\t\U2022 ",
        crayon::blue(n_responses),
        " responses\n\t\U2022 ",
        crayon::blue(moi),
        " multiplicity-of-infection \n\t\U2022 ",
        crayon::blue(guide_counts$n_targeting_grnas),
        " targeting gRNAs (distributed across ",
        crayon::blue(guide_counts$n_targeting_grna_targets),
        " targets) \n\t\U2022 ",
        crayon::blue(guide_counts$n_nt_grnas),
        " non-targeting gRNAs \n\t\U2022 ",
        crayon::blue(n_covariates),
        " covariates (",
        covariates,
        ")"
    ))
})


#' Print
#'
#' `print()` prints information about the dataset and the status of the analysis
#' to the console. The output contains several fields: Attributes of the data
#' (summarizing key features of the data), Analysis status (indicating the
#' analysis functions that have been called), Analysis parameters (summarizing
#' the analysis parameters set in `set_analysis_parameters()`), gRNA-to-cell
#' assignment information (summarizing the outcome of the gRNA-to-cell
#' assignment step), and Summary of results (summarizing the key analysis
#' results). A subset of these fields may be printed, depending on the status of
#' the analysis.
#'
#' @param x a `sceptre_object`
#' @return the value NULL
#' @export
#' @examples
#'
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' print(sceptre_object)
setMethod("print", signature = signature("sceptre_object"), function(x) {
    x <- migrate_sceptre_object(x)
    show(x)
    print_alternative_inference(x)

    # 1. print analysis status
    funct_run_vect <- workflow_status(x)
    get_mark <- function(bool) {
        ifelse(bool, crayon::green("\u2713"), crayon::red("\u2717"))
    }
    cat("\n\nAnalysis status:\n")
    for (i in seq_along(funct_run_vect)) {
        cat(paste0(
            "\t",
            get_mark(funct_run_vect[i]),
            " ",
            names(funct_run_vect)[i],
            "()\n"
        ))
    }

    # 2. print analysis parameters
    n_discovery_pairs <- if (funct_run_vect[["set_analysis_parameters"]]) {
        pair_count(x, "discovery")
    } else {
        integer()
    }
    n_pc_pairs <- if (funct_run_vect[["set_analysis_parameters"]]) {
        pair_count(x, "positive_control")
    } else {
        integer()
    }
    cat(paste0(
        "\nAnalysis parameters: \n",
        "\t\U2022 Discovery pairs:",
        if (length(n_discovery_pairs) == 0) {
            " not specified"
        } else {
            paste0(
                " data frame with ",
                crayon::blue(n_discovery_pairs),
                " pairs",
                if (funct_run_vect["run_qc"] && n_discovery_pairs >= 1L) {
                    paste0(
                        " (",
                        crayon::blue(pair_count(
                            x,
                            "discovery",
                            passing = TRUE
                        )),
                        " after pairwise QC)"
                    )
                } else {
                    NULL
                }
            )
        },
        "\n\t\U2022 Positive control pairs:",
        if (length(n_pc_pairs) == 0) {
            " not specified"
        } else {
            paste0(
                " data frame with ",
                crayon::blue(n_pc_pairs),
                " pairs",
                if (funct_run_vect["run_qc"] && n_pc_pairs >= 1L) {
                    paste0(
                        " (",
                        crayon::blue(pair_count(
                            x,
                            "positive_control",
                            passing = TRUE
                        )),
                        " after pairwise QC)"
                    )
                } else {
                    NULL
                }
            )
        },
        "\n\t\U2022 Sidedness of test: ",
        if (length(x@spec@side_code) == 0L) {
            "not specified"
        } else {
            crayon::blue(c("left", "both", "right")[x@spec@side_code + 2L])
        },
        if (!x@data@low_moi) {
            NULL
        } else {
            paste0(
                "\n\t\U2022 Control group: ",
                if (length(x@spec@control_group_complement) == 0L) {
                    "not specified"
                } else {
                    crayon::blue(ifelse(
                        x@spec@control_group_complement,
                        "complement set",
                        "non-targeting cells"
                    ))
                }
            )
        },
        "\n\t\U2022 Resampling mechanism: ",
        if (length(x@spec@run_permutations) == 0L) {
            "not specified"
        } else {
            crayon::blue(ifelse(
                x@spec@run_permutations,
                "permutations",
                "conditional resampling"
            ))
        },
        "\n\t\U2022 gRNA integration strategy: ",
        if (length(x@spec@grna_integration_strategy) == 0L) {
            "not specified"
        } else {
            crayon::blue(x@spec@grna_integration_strategy)
        },
        "\n\t\U2022 Resampling approximation: ",
        if (length(x@spec@resampling_approximation) == 0L) {
            "not specified"
        } else {
            crayon::blue(gsub(
                pattern = "_",
                replacement = " ",
                fixed = TRUE,
                x = x@spec@resampling_approximation
            ))
        },
        "\n\t\U2022 Multiple testing adjustment: ",
        if (length(x@spec@multiple_testing_method) == 0L) {
            "none"
        } else {
            paste0(
                crayon::blue(x@spec@multiple_testing_method),
                " at level ",
                crayon::blue(x@spec@multiple_testing_alpha)
            )
        },
        "\n\t\U2022 N nonzero treatment cells threshold: ",
        if (length(x@spec@n_nonzero_trt_thresh) == 0L) {
            "not specified"
        } else {
            crayon::blue(x@spec@n_nonzero_trt_thresh)
        },
        "\n\t\U2022 N nonzero control cells threshold: ",
        if (length(x@spec@n_nonzero_cntrl_thresh) == 0L) {
            "not specified"
        } else {
            crayon::blue(x@spec@n_nonzero_cntrl_thresh)
        },
        "\n\t\U2022 Formula object: ",
        if (length(x@spec@formula_object) == 0L) {
            "not specified"
        } else {
            crayon::blue(as.character(x@spec@formula_object)[2])
        }
    ))

    # 3. print the gRNA-to-cell assignment information
    grna_assignment_run <- funct_run_vect[["assign_grnas"]]
    if (grna_assignment_run) {
        mean_cells_per_grna <- mean_assignment_size(x)
        cat(
            paste0(
                "\n\ngRNA-to-cell assignment information:",
                "\n\t\U2022 Assignment method: ",
                crayon::blue(x@spec@grna_assignment_method),
                "\n\t\U2022 Mean N cells per gRNA: ",
                crayon::blue(mean_cells_per_grna |> round(2)),
                "\n\t\U2022 Mean N gRNAs per cell (MOI): ",
                if (x@spec@grna_assignment_method == "maximum") {
                    "not computed when using \"maximum\" assignment method"
                } else {
                    crayon::blue(
                        x@preparation@guide_counts_per_cell |>
                            mean() |>
                            round(2)
                    )
                }
            ),
            if (x@spec@grna_assignment_method == "mixture") {
                paste0(
                    "\n\t\U2022 gRNA assignment formula object: ",
                    crayon::blue(as.character(
                        x@spec@grna_assignment_hyperparameters$formula_object
                    )[2])
                )
            } else {
                NULL
            }
        )
    }

    # 4. print the results summary
    calib_check_run <- nrow(current_result(
        x,
        "run_calibration_check",
        required = FALSE
    )) >
        0L
    discovery_analysis_run <- nrow(current_result(
        x,
        "run_discovery_analysis",
        required = FALSE
    )) >
        0L
    power_check_run <- nrow(current_result(
        x,
        "run_power_check",
        required = FALSE
    )) >
        0L
    if (calib_check_run || discovery_analysis_run || power_check_run) {
        cat("\n\nSummary of results:")
    }
    if (calib_check_run) {
        n_false_rejections <- sum(
            current_result(x, "run_calibration_check")$significant
        )
        mean_lfc <- signif(
            mean(current_result(x, "run_calibration_check")$log_2_fold_change),
            2
        )
        n_calib_pairs <- nrow(current_result(x, "run_calibration_check"))
        cat(paste0(
            "\n\t\U2022 N",
            crayon::red(" negative control "),
            "pairs called as significant: ",
            crayon::blue(paste0(n_false_rejections, "/", n_calib_pairs)),
            "\n\t\U2022 Mean log-2 FC for",
            crayon::red(" negative control "),
            "pairs: ",
            crayon::blue(mean_lfc)
        ))
    }
    if (power_check_run) {
        median_pc_p_value <- stats::median(
            current_result(x, "run_power_check")$p_value,
            na.rm = TRUE
        )
        cat(paste0(
            "\n\t\U2022 Median",
            crayon::green(" positive control "),
            "p-value: ",
            crayon::blue(signif(median_pc_p_value, 2))
        ))
    }
    if (discovery_analysis_run) {
        n_discoveries <- sum(
            current_result(x, "run_discovery_analysis")$significant,
            na.rm = TRUE
        )
        n_discovery_pairs <- pair_count(x, "discovery", passing = TRUE)
        cat(paste0(
            "\n\t\U2022 N",
            crayon::yellow(" discovery pairs "),
            "called as significant: ",
            crayon::blue(paste0(n_discoveries, "/", n_discovery_pairs))
        ))
    }
})


#' Plot
#'
#' `plot()` creates a plot depicting the current state of a `sceptre_object`.
#'
#' `plot()` is "generic" in the sense that it dispatches a specific plotting
#' function based on completed current components. It prefers discovery,
#' power, then calibration results, followed by QC and guide assignments.
#' Historical results do not select a plot. Users
#' can pass arguments to the function dispatched by `plot()` as named arguments
#' to `plot()`.
#'
#' @param x a `sceptre_object`
#' @param y ignored argument
#' @param ... arguments passed to the plotting function dispatched by `plot()`
#' @return a single \code{cowplot} object containing the combined panels (if
#' \code{return_indiv_plots} is set to \code{FALSE}) or a list of the individual
#' panels (if \code{return_indiv_plots} is set to \code{TRUE})
#'
#' @export
#' @examples
#'
#' data(highmoi_example_data)
#' data(grna_target_data_frame_highmoi)
#' # import data
#' sceptre_object <- import_data(
#'   response_matrix = highmoi_example_data$response_matrix,
#'   grna_matrix = highmoi_example_data$grna_matrix,
#'   grna_target_data_frame = grna_target_data_frame_highmoi,
#'   moi = "high",
#'   extra_covariates = highmoi_example_data$extra_covariates,
#'   response_names = highmoi_example_data$gene_names
#' )
#' # set analysis parameters, assign grnas
#' sceptre_object <- sceptre_object |>
#'   set_analysis_parameters() |>
#'   assign_grnas(method = "thresholding") |>
#'   plot()
setMethod("plot", signature = signature("sceptre_object"), function(x, y, ...) {
    x <- migrate_sceptre_object(x)
    args <- list(...)
    args[["sceptre_object"]] <- x

    status <- workflow_status(x)
    candidates <- c(
        "run_discovery_analysis",
        "run_power_check",
        "run_calibration_check",
        "run_qc",
        "assign_grnas"
    )
    available <- candidates[status[candidates]]
    if (!length(available)) {
        stop(
            "There is no generic plot configured for the current prepared ",
            "state."
        )
    }
    funct_to_call <- paste0("plot_", available[[1L]])
    p <- do.call(what = funct_to_call, args = args)
    return(p)
})


#' Read all cells for one named gene or guide
#' @param mat Row-accessible sparse matrix (dgRMatrix) or open ondisc matrix
#'   (odm)
#'   with unique row names.
#' @param id Gene/guide name present in rownames(mat).
#' @return Dense numeric vector in original cell-column order. Other matrix
#'   classes return NULL; callers must prepare row access first. Used by mixture
#'   fitting and plotting; the supplied matrix is left unchanged.
#' @noRd
load_row <- function(mat, id) {
    if (methods::is(mat, "odm") || methods::is(mat, "dgRMatrix")) {
        storage_read_response(mat, storage_response_index(mat, id))
    }
}

screen_guide_counts <- function(grna_target_data_frame) {
    n_nt_grnas <- grna_target_data_frame |>
        dplyr::filter(grna_target == "non-targeting") |>
        nrow()
    targeting_grnas_df <- grna_target_data_frame |>
        dplyr::filter(grna_target != "non-targeting")
    n_targeting_grna_targets <- length(unique(targeting_grnas_df$grna_target))
    n_targeting_grnas <- nrow(targeting_grnas_df)
    list(
        n_nt_grnas = n_nt_grnas,
        n_targeting_grna_targets = n_targeting_grna_targets,
        n_targeting_grnas = n_targeting_grnas
    )
}
