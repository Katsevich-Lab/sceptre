test_that("legacy eligibility preserves boundaries missing counts and input order", {
    counts <- data.frame(n_nonzero_trt = c(6L, 7L, 9L, NA, NA, 9L, 9L),
                         n_nonzero_cntrl = c(9L, 7L, 6L, 6L, 9L, NA, 9L))
    before <- serialize(counts, NULL)
    options <- list(n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L)
    expected <- c(FALSE, TRUE, FALSE, FALSE, NA, NA, TRUE)
    expect_identical(evaluate_pair_eligibility(counts, options), expected)
    expect_identical(evaluate_pair_eligibility(counts[7:1, ], options), rev(expected))
    expect_identical(evaluate_pair_eligibility(counts[FALSE, ], options), logical())
    expect_identical(serialize(counts, NULL), before)
})

test_that("calibration candidate acceptance uses the same count policy", {
    # Use one candidate guide group so its output label is unambiguous. Keep
    # the sampler's random candidate order and accept counts equal to thresholds.
    treatment <- matrix(c(6L, 7L, 10L, 7L), nrow = 1L)
    total <- c(16L, 14L, 16L, 16L)
    groups <- matrix(0L, 1L, 1L)
    before <- serialize(list(treatment, total, groups), NULL)
    sampled <- sample_calibration_pairs(treatment, total, groups,
        n_genes = 4L, n_calibration_pairs = 4L,
        n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L)
    metrics <- list(n_nonzero_trt = as.integer(treatment[1L, ]),
                    n_nonzero_cntrl = total - treatment[1L, ])
    eligible <- evaluate_pair_eligibility(metrics,
        list(n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L))
    expect_identical(sort(sampled$response_idxs), which(eligible))
    expect_identical(sampled$n_nonzero_trt_v, metrics$n_nonzero_trt[sampled$response_idxs])
    expect_identical(sampled$n_nonzero_cntrl_v, metrics$n_nonzero_cntrl[sampled$response_idxs])
    expect_identical(serialize(list(treatment, total, groups), NULL), before)
})
