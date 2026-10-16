args <- commandArgs(trailingOnly = TRUE)
candidate <- args[[1L]]
evidence <- args[[2L]]
source(file.path(candidate, "R/s4_class.R"))
core <- data.frame(
    response_id = c("response_a", "response_b"),
    grna_target = c("target_a", "target_b"),
    n_nonzero_trt = c(15L, 18L), n_nonzero_cntrl = c(20L, 21L),
    pass_qc = TRUE, p_value = c(0.0001, 0.5),
    fold_change = c(0.5, 1.1), se_fold_change = c(0.1, 0.2),
    log_2_fold_change = log2(c(0.5, 1.1)),
    significant = c(TRUE, FALSE)
)
for (revision in c("baseline21", "cc01683")) {
    env <- new.env(parent = globalenv())
    file <- if (revision == "baseline21") {
        file.path(evidence, "baseline-plotting_functions.R")
    } else file.path(candidate, "R/plotting_functions.R")
    qq_file <- if (revision == "baseline21") {
        file.path(evidence, "baseline-qq_plot_helpers.R")
    } else file.path(candidate, "R/qq_plot_helpers.R")
    sys.source(qq_file, envir = env)
    sys.source(file, envir = env)
    for (amount in 1:3) {
        object <- methods::new("sceptre_object")
        object@functs_called <- c(run_discovery_analysis = TRUE)
        object@multiple_testing_alpha <- 0.1
        object@discovery_result <- core
        if (amount >= 2L) {
            # This is the native stage-1/empirical-stage-3 diagnostic schema.
            object@discovery_result$stage <- c(1L, 3L)
            object@discovery_result$z_orig <- c(-4, 0.25)
            object@discovery_result$xi <- NA_real_
            object@discovery_result$omega <- NA_real_
            object@discovery_result$alpha <- NA_real_
        }
        if (amount == 3L) {
            object@discovery_result$z_null_1 <- c(-1.1, 0.8)
            object@discovery_result$z_null_2 <- c(NA_real_, -0.2)
        }
        plots <- env$plot_run_discovery_analysis(object, return_indiv_plots = TRUE)
        volcano <- plots[[3L]]$data
        label <- ggplot2::ggplot_build(plots[[4L]])$data[[1L]]$label
        cat(revision, "output_amount=", amount,
            "result_rows=", nrow(object@discovery_result),
            "result_significant=", sum(object@discovery_result$significant),
            "volcano_rows=", nrow(volcano),
            "volcano_significant=", sum(volcano$significant), "\n", sep = " ")
        cat("summary=", label, "\n", sep = "")
        stopifnot(nrow(volcano) == if (amount == 1L) 2L else 0L)
        stopifnot(grepl(if (amount == 1L) "1 of 2" else "0 of 0", label, fixed = TRUE))
    }
}
cat("Source-evaluated reproduction completed; no package install, build, fixture regeneration or scientific analysis was performed.\n")
print(sessionInfo())
