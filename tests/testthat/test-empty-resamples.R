# A conditional randomization resample may contain no treated cells (the
# number of treated cells per resample is random). The C++ test statistic code
# took the address of the first element of such an empty vector, which is
# undefined behavior and aborts R where C++ library assertions are enabled
# (for example, packages built with -D_GLIBCXX_ASSERTIONS, as on Fedora). This
# test exercises that path; it fails only under such a build or a sanitizer.

test_that("CRT resamples with no treated cells are handled", {
    set.seed(1)
    n <- 200L
    y <- stats::rpois(n, 3)
    mu <- rep(3, n)
    a <- stats::rnorm(n)
    w <- stats::runif(n) + 1
    D <- matrix(stats::rnorm(2L * n), 2L, n)
    # P(a resample has no treated cell) = 0.995^200, about 0.37
    synthetic_idxs <- crt_index_sampler_fast(rep(0.005, n), B = 600L, seed = 1L)
    out <- run_low_level_test_full_v4(
        y, mu, a, w, D,
        trt_idxs = 1:3, n_trt = 3L, use_all_cells = TRUE,
        synthetic_idxs = synthetic_idxs,
        B1 = 100L, B2 = 0L, B3 = 500L,
        fit_parametric_curve = FALSE, return_resampling_dist = TRUE,
        side_code = 0L
    )
    expect_true(is.finite(out$p))
    expect_true(out$p >= 0 && out$p <= 1)
})
