# Hand-calculated small families exercise integration branches without a DGP.
# These are regression/algebra checks, not statistical calibration evidence.
coverage_math_input <- function() {
    data.frame(response_id = rep(c("a", "b", "c"), each = 3L),
        grna_group = paste0("g", 1:9),
        p_value = c(.4, .7, NA, NA, NA, NA, .01, .04, NA),
        pass_qc = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE),
        n_nonzero_trt = c(10L, 20L, 99L, 4L, 7L, 5L, 8L, 9L, 100L),
        n_nonzero_cntrl = c(30L, 40L, 99L, 9L, 6L, 2L, 11L, 12L, 100L),
        log_2_fold_change = c(-1, -2, NA, NA, NA, NA, -3, -4, NA))
}

coverage_math_map <- function() data.frame(grna_id = paste0("g", 1:9),
    grna_target = rep(c("A", "B", "C"), each = 3L))

test_that("guide correction caps p-values and retains all-failed QC counts", {
    raw <- coverage_math_input()
    before <- serialize(raw, NULL)
    actual <- integrate_guides(raw, "bonferroni", coverage_math_map())
    actual <- actual[match(c("a", "b", "c"), actual$response_id), ]
    # Only passing guides contribute to m: 2*.4=.8; 2*.01=.02.
    expect_equal(actual$p_value, c(.8, NA, .02))
    expect_identical(actual$pass_qc, c(TRUE, FALSE, TRUE))
    expect_equal(actual$n_nonzero_trt, c(10, 7, 8))
    expect_equal(actual$n_nonzero_cntrl, c(30, 9, 11))
    expect_equal(actual$log_2_fold_change, c(-1, NA, -3))
    expect_identical(serialize(raw, NULL), before)
    raw$p_value[1:2] <- c(.8, .9)
    capped <- integrate_guides(raw, "bonferroni", coverage_math_map())
    expect_equal(capped$p_value[match("a", capped$response_id)], 1)
    calibration <- integrate_guides(raw, "bonferroni", coverage_math_map(), TRUE)
    expect_equal(nrow(calibration), 9L)
    expect_equal(calibration$p_value[match(raw$grna_group, calibration$grna_id)], raw$p_value)
})

test_that("complete-family BH uses the scalar cutoff and strict boundary", {
    # Sorted finite p=(1/32,1/8,1/2), finite family n=3 (one NA).
    # BH adjusted values are (3/32,3/16,1/2,NA). At alpha=3/16,
    # the second row must fail the strict '<', even with tail-shape alpha=1.
    raw <- data.frame(response_id = letters[1:4], grna_group = LETTERS[1:4],
        p_value = c(1/32, 1/8, 1/2, NA), alpha = c(1, 1, 1, 1), significant = TRUE)
    actual <- finalize_association_result(raw, "union", coverage_math_map(), "BH", 3/16)
    check <- function(x) expect_identical(x$significant, c(TRUE, FALSE, FALSE, NA),
        info = "whole family, scalar cutoff and strict boundary")
    check(actual)
    # A per-worker/raw-p decision and a <= boundary decision both look plausible;
    # these assertions must demonstrably reject them.
    bad <- actual
    bad$significant <- c(TRUE, TRUE, FALSE, NA)
    expect_failure(check(bad), "whole family")
    bad$significant <- c(TRUE, TRUE, TRUE, NA)
    expect_failure(check(bad), "whole family")
    bonf <- finalize_association_result(raw, "union", coverage_math_map(), "bonferroni", 3/16)
    expect_identical(bonf$significant, c(TRUE, FALSE, FALSE, NA))
})
