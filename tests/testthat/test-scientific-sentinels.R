# A small actual-outcome power/direction sentinel. Offline statistical precision
# belongs in tools/validation; this fixed dataset is not a calibration estimate.
test_that("independent count effects retain power and direction with default inference", {
    withr::with_seed(20261020, {
        n <- 400L
        ids <- paste0("cell", seq_len(n))
        groups <- sample(rep(1:4, each = n / 4L))
        guides <- matrix(0, 4L, n, dimnames = list(c("a", "nt1", "nt2", "b"), ids))
        guides[cbind(groups, seq_len(n))] <- 50
        treated <- groups == 1L
        response <- rbind(up = stats::rpois(n, ifelse(treated, 12, 3)),
                          down = stats::rpois(n, ifelse(treated, 1, 8)),
                          background = stats::rpois(n, 10))
        colnames(response) <- ids
        object <- suppressMessages(import_data(response, guides,
            data.frame(grna_id = rownames(guides),
                grna_target = c("A", "non-targeting", "non-targeting", "B")),
            moi = "low", use_ondisc = FALSE))
        object <- suppressMessages(set_analysis_parameters(object,
            discovery_pairs = data.frame(grna_target = "A", response_id = c("up", "down")),
            formula_object = ~ 1))
        object <- suppressMessages(assign_grnas(object, method = "thresholding", print_progress = FALSE))
        object <- suppressMessages(run_qc(object, response_n_umis_range = c(0, 1),
            response_n_nonzero_range = c(0, 1)))
        object <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE))
        result <- get_result(object, "run_discovery_analysis")
        expect_true(all(result$pass_qc))
        expect_true(all(is.finite(result$p_value) & result$p_value >= 0 & result$p_value < .01))
        expect_gt(result$log_2_fold_change[result$response_id == "up"], 0)
        expect_lt(result$log_2_fold_change[result$response_id == "down"], 0)
    })
})
