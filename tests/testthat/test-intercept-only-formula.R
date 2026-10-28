# Score statistic of an intercept-only fit, computed directly: with design
# matrix 1, the fitted mean is the sample mean and the score for treatment
# cells trt among n cells is sum(a[trt]) / sqrt(w (k - k^2 / n)).
intercept_only_score <- function(y, trt) {
    mu <- mean(y)
    theta <- estimate_theta(
        y,
        rep(mu, length(y)),
        length(y) - 1L,
        limit = 50,
        eps = .Machine$double.eps^(1 / 4)
    )[[1L]]
    theta <- max(min(theta, 1000), 0.01)
    denom <- 1 + mu / theta
    k <- length(trt)
    sum((y[trt] - mu) / denom) / sqrt(mu / denom * (k - k^2 / length(y)))
}

test_that("intercept-only formula works with NT controls (calibration + discovery)", {
    inputs <- refactor_make_inputs()
    for (mechanism in c("permutations", "crt")) {
        case <- refactor_case_definitions()$low_nt_perm
        case$mechanism <- mechanism
        object <- refactor_import(inputs, case)
        object <- sceptre::set_analysis_parameters(
            object,
            discovery_pairs = inputs$discovery_pairs,
            positive_control_pairs = inputs$positive_control_pairs,
            formula_object = ~1,
            control_group = "nt_cells",
            resampling_mechanism = mechanism
        )
        expect_identical(ncol(object@covariate_matrix), 1L)
        object <- sceptre::assign_grnas(
            object,
            method = "thresholding",
            threshold = 40,
            print_progress = FALSE
        )
        object <- sceptre::run_qc(
            object,
            response_n_umis_range = c(0, 1),
            response_n_nonzero_range = c(0, 1),
            n_nonzero_trt_thresh = 7L,
            n_nonzero_cntrl_thresh = 7L,
            additional_cells_to_remove = inputs$remove_cells
        )
        object <- suppressMessages(sceptre::run_calibration_check(
            object,
            n_calibration_pairs = 4L,
            calibration_group_size = 2L,
            print_progress = FALSE
        ))
        expect_gt(nrow(object@calibration_result), 0L)
        object <- suppressMessages(sceptre::run_discovery_analysis(
            object,
            print_progress = FALSE,
            output_amount = 2L
        ))
        result <- object@discovery_result
        result <- result[result$pass_qc, ]
        expect_gt(nrow(result), 0L)
        expect_true(all(is.finite(result$p_value) & result$p_value <= 1))
        # each score equals the intercept-only score over the gRNA's treatment
        # cells followed by the NT cells
        assignments <- object@grna_assignments
        for (row in seq_len(nrow(result))) {
            trt <- assignments$grna_group_idxs[[result$grna_target[row]]]
            cells <- c(trt, assignments$all_nt_idxs)
            y <- inputs$response_matrix[
                result$response_id[row],
                object@cells_in_use
            ][cells]
            expect_equal(
                result$z_orig[row],
                intercept_only_score(y, seq_along(trt)),
                tolerance = 1e-8
            )
        }
    }
})
