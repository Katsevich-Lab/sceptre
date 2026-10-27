plot_diagnostic_object <- function(output_amount = 1L, as_data_table = FALSE) {
    object <- methods::new("sceptre_object")
    object@functs_called <- c(run_discovery_analysis = TRUE)
    object@multiple_testing_alpha <- 0.1
    object@discovery_result <- data.frame(
        response_id = c("response_a", "response_b"),
        grna_target = c("target_a", "target_b"),
        n_nonzero_trt = c(15L, 18L),
        n_nonzero_cntrl = c(20L, 21L),
        pass_qc = TRUE,
        p_value = c(0.0001, 0.5),
        fold_change = c(0.5, 1.1),
        se_fold_change = c(0.1, 0.2),
        log_2_fold_change = log2(c(0.5, 1.1)),
        significant = c(TRUE, FALSE)
    )
    object@calibration_result <- data.frame(p_value = c(0.1, 0.8))
    if (output_amount >= 2L) {
        # Stage 1 and empirical stage 3 legitimately have no fitted SN curve.
        object@discovery_result$stage <- c(1L, 3L)
        object@discovery_result$z_orig <- c(-4, 0.25)
        object@discovery_result$xi <- NA_real_
        object@discovery_result$omega <- NA_real_
        object@discovery_result$alpha <- NA_real_
    }
    if (output_amount >= 3L) {
        object@discovery_result$z_null_1 <- c(-1.1, 0.8)
        object@discovery_result$z_null_2 <- c(NA_real_, -0.2)
    }
    if (as_data_table) {
        object@discovery_result <- data.table::as.data.table(
            object@discovery_result
        )
    }
    object
}

plot_diagnostic_data <- function(panels, core) {
    expected_qq <- data.frame(
        p_value = c(0.1, 0.8, core$p_value),
        lab = factor(
            c(rep("Negative control", 2L), rep("Discovery", nrow(core))),
            levels = c("Discovery", "Negative control")
        )
    )
    expect_identical(panels[[1L]]$data, expected_qq)
    expect_identical(panels[[2L]]$data, expected_qq)
    volcano <- as.data.frame(panels[[3L]]$data)
    expect_setequal(volcano$response_id, core$response_id)
    expect_identical(nrow(volcano), nrow(core))
    columns <- c("response_id", "p_value", "log_2_fold_change", "significant")
    actual <- volcano[match(core$response_id, volcano$response_id), columns]
    rownames(actual) <- NULL
    expected <- core[, columns]
    rownames(expected) <- NULL
    expect_identical(actual, expected)
    expect_identical(volcano$reject, volcano$significant)
    expect_identical(panels[[3L]]$layers[[2L]]$data$yintercept, 0.0001)
    built <- lapply(panels, function(plot) ggplot2::ggplot_build(plot)$data)
    expect_identical(
        built[[4L]][[1L]]$label,
        "Number of discovery pairs\ncalled as significant\n(at alpha 0.1): 1 of 2"
    )
    built
}

test_that("diagnostic output amounts retain the same discovery plot family", {
    withr::local_seed(581L)
    seed <- .Random.seed
    plot_discovery <- plot_run_discovery_analysis
    core <- plot_diagnostic_object()@discovery_result
    built_by_class <- list()
    for (as_data_table in c(FALSE, TRUE)) {
        built <- vector("list", 3L)
        for (amount in 1:3) {
            object <- plot_diagnostic_object(amount, as_data_table)
            retained_result <- object@discovery_result
            before <- serialize(object, NULL)
            result_before <- serialize(retained_result, NULL)
            panels <- plot_discovery(object, return_indiv_plots = TRUE)
            expect_length(panels, 4L)
            built[[amount]] <- plot_diagnostic_data(panels, core)
            expect_identical(serialize(object, NULL), before)
            expect_identical(serialize(retained_result, NULL), result_before)
            expect_identical(.Random.seed, seed)
        }
        # Compare actual QQ points/bands, volcano points/threshold and text data,
        # not just the input result count or the unbuilt ggplot definitions.
        expect_identical(built[[2L]], built[[1L]])
        expect_identical(built[[3L]], built[[1L]])
        built_by_class[[length(built_by_class) + 1L]] <- built
    }
    expect_identical(built_by_class[[2L]], built_by_class[[1L]])
})

test_that("missing primary plot fields still exclude unplottable discovery rows", {
    withr::local_seed(852L)
    seed <- .Random.seed
    core <- plot_diagnostic_object(3L)@discovery_result
    extra <- core[rep(2L, 4L), ]
    extra$response_id <- c(
        "missing_p",
        "missing_effect",
        "missing_flag",
        "failed_qc"
    )
    extra$p_value[c(1L, 4L)] <- NA_real_
    extra$log_2_fold_change[c(2L, 4L)] <- NA_real_
    extra$significant[c(3L, 4L)] <- NA
    extra$pass_qc[4L] <- FALSE
    for (as_data_table in c(FALSE, TRUE)) {
        object <- plot_diagnostic_object(3L)
        object@discovery_result <- rbind(core, extra)
        if (as_data_table) {
            object@discovery_result <- data.table::as.data.table(
                object@discovery_result
            )
        }
        retained_result <- object@discovery_result
        before <- serialize(object, NULL)
        result_before <- serialize(retained_result, NULL)
        panels <- plot_run_discovery_analysis(
            object,
            return_indiv_plots = TRUE
        )
        plot_diagnostic_data(panels, core)
        expect_identical(serialize(object, NULL), before)
        expect_identical(serialize(retained_result, NULL), result_before)
        expect_identical(.Random.seed, seed)
    }
})
