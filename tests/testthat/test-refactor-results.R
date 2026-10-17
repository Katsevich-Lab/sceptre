test_that("result finalization integrates guides before complete-family correction", {
    raw <- data.frame(response_id = c("r1", "r1", "r2"),
        grna_group = c("a1", "a2", "b1"),
        n_nonzero_trt = c(11L, 12L, 13L), n_nonzero_cntrl = c(21L, 22L, 23L),
        log_2_fold_change = c(-1, -2, 1), pass_qc = TRUE,
        p_value = c(0.01, 0.04, 0.2), alpha = NA_real_)
    targets <- data.frame(grna_id = c("a1", "a2", "b1"),
                          grna_target = c("A", "A", "B"))
    before <- serialize(list(raw, targets), NULL)
    result <- finalize_association_result(raw, "bonferroni", targets, "BH", 0.05)
    expect_identical(result$grna_target, c("A", "B"))
    expect_equal(result$p_value, c(0.02, 0.2))
    expect_identical(result$significant, c(TRUE, FALSE))
    expect_identical(result$n_nonzero_trt, c(11L, 13L))
    expect_equal(result$log_2_fold_change, c(-1, 1))
    expect_identical(serialize(list(raw, targets), NULL), before)
    raw$pass_qc[2L] <- FALSE
    raw$p_value[2L] <- NA_real_
    partial <- finalize_association_result(raw, "bonferroni", targets, "BH", 0.05)
    expect_equal(partial$p_value, c(0.01, 0.2))
    expect_identical(serialize(targets, NULL), serialize(unserialize(before)[[2L]], NULL))
    # Adjusting each worker's output alone would call p = 0.04 significant.
    # Across both pairs, BH adjusts it to 0.08; recompute significance after
    # collection instead of retaining that worker's flag.
    union <- raw[c(1L, 3L), ]
    union$p_value <- c(0.04, 0.2)
    union$significant <- c(TRUE, FALSE)
    final <- finalize_association_result(union, "union", targets, "BH", 0.05)
    expect_identical(final$significant, c(FALSE, FALSE))
})

test_that("raw assembly retains failed QC identity and native output fields", {
    requested <- data.frame(response_id = c("r2", "r1", "r3"),
        grna_group = c("B", "A", "C"), pass_qc = c(TRUE, TRUE, FALSE))
    record <- list(p = 0.04, fc = 2, se = 0.2, stage = 3L,
                   z_orig = 2, sn_params = rep(NA_real_, 3L),
                   resampling_dist = c(-1, 0, 1))
    raw <- pair_results_from_native(requested[1:2, ], list(record, record), 3L)
    result <- assemble_pair_results(raw, requested)
    expect_identical(result$response_id, c("r1", "r2", "r3"))
    expect_identical(result$pass_qc, c(TRUE, TRUE, FALSE))
    expect_true(is.na(result$p_value[3L]))
    expect_true(is.na(result$fold_change[3L]))
    expect_equal(result$se_fold_change[1:2], c(0.2, 0.2))
    expect_equal(result$log_2_fold_change[1:2], c(1, 1))
    expect_identical(result$stage[1:2], c(3L, 3L))
    expect_named(result, c(names(requested), "p_value", "fold_change",
        "se_fold_change", "log_2_fold_change", "stage", "z_orig", "xi",
        "omega", "alpha", "z_null_1", "z_null_2", "z_null_3"))
})
