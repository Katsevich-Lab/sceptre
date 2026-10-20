#!/usr/bin/env Rscript
# Test internal C++ functions directly from the current source checkout. The
# installed-package suite also checks stages, effect formulas and sampling laws.
# These tests use no old implementation or saved historical outputs.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
stopifnot(file.exists(file.path(root, "src", "inference.h")))
environment <- new.env(parent = globalenv())
bridge <- tempfile("sceptre-native-components-", fileext = ".cpp")
modules <- c(SCORES = "inference_scores.cpp", TAIL = "inference_tail.cpp",
             STAGES = "inference_stages.cpp", LIFETIME = "result_lifetime.h")
prefix <- sprintf("#define SCEPTRE_COMPONENT_%s %s", names(modules),
                  vapply(file.path(root, "src", modules), encodeString,
                         character(1), quote = '"'))
writeLines(c(prefix, readLines(file.path(root, "tools", "checks", "native-components.cpp"))),
           bridge)
Rcpp::sourceCpp(bridge, env = environment, rebuild = TRUE, showOutput = FALSE)
unlink(bridge)

testthat::test_that("effect and scores have independent explicit inputs", {
    y <- c(4, 12, 6, 2)
    mu <- c(3, 7, 5, 1)
    treatment <- c(4L, 2L)
    selected_y <- y[treatment]
    selected_mu <- mu[treatment]
    fc <- mean(selected_y) / mean(selected_mu)
    se <- sqrt(sum((selected_y - fc * selected_mu)^2) / sum(selected_mu)^2)
    testthat::expect_equal(environment$component_effect(y, mu, treatment - 1L, 2L), c(fc, se))
    a <- c(0.4, 1.1, -0.5, 0.2)
    w <- c(2, 3, 2, 4)
    D <- rbind(c(0.1, 0.2, 0.1, -0.2), c(0.2, 0.1, -0.1, 0.1))
    score <- function(index) sum(a[index]) /
        sqrt(sum(w[index]) - sum(rowSums(D[, index, drop = FALSE])^2))
    testthat::expect_equal(environment$component_observed(a, w, D, treatment - 1L), score(treatment))
    draws <- list(c(0L, 2L, 3L), c(3L, 1L), c(1L, 2L, 0L), c(2L, 0L))
    before <- serialize(list(a, w, D, treatment, draws), NULL)
    # Draw offsets and cell positions both start at 0. A hybrid draw can serve
    # several group sizes; each test reads only its requested first k positions.
    fixed <- environment$component_segment(a, w, D, draws, 1L, 2L, 2L, FALSE)
    variable <- environment$component_segment(a, w, D, draws, 1L, 2L, 2L, TRUE)
    testthat::expect_equal(fixed, vapply(draws[2:3], function(x) score(x[1:2] + 1L), numeric(1)))
    testthat::expect_equal(variable, vapply(draws[2:3], function(x) score(x + 1L), numeric(1)))
    testthat::expect_length(environment$component_segment(a, w, D, draws, 2L, 0L, 2L, TRUE), 0L)
    testthat::expect_identical(serialize(list(a, w, D, treatment, draws), NULL), before)
})

testthat::test_that("empirical tails include ties and bound doubled tails at one", {
    scores <- c(-1, 0, 0, 2)
    testthat::expect_equal(environment$component_empirical(scores, 0, -1L), 4/5)
    testthat::expect_equal(environment$component_empirical(scores, 0, 1L), 4/5)
    testthat::expect_equal(environment$component_empirical(scores, 0, 0L), 1)
})

# Independent tail-count oracle, deliberately including zero/one draw, all ties,
# mixed ties, both extremes, and both tails. The bound changes no rejection for
# alpha < 1 and does not use an RNG or mutate the scores.
testthat::test_that("empirical bounds preserve counts, inputs, RNG and rejections", {
    set.seed(237L)
    seed <- .Random.seed
    for (scores in list(numeric(), 0, rep(0, 101), c(-1, 0, 0, 2),
                        seq(-3, 3, length.out = 100))) {
        before <- serialize(scores, NULL)
        for (observed in c(-100, -1, 0, 1, 100)) {
            left <- (1 + sum(scores <= observed)) / (1 + length(scores))
            right <- (1 + sum(scores >= observed)) / (1 + length(scores))
            expected <- c(left, min(1, 2 * min(left, right)), right)
            actual <- vapply(-1:1, function(side)
                environment$component_empirical(scores, observed, side), numeric(1))
            testthat::expect_identical(actual, expected)
            testthat::expect_true(all(actual >= 0 & actual <= 1))
            for (alpha in c(0, 0.01, 0.05, 0.1, 0.99)) {
                testthat::expect_identical(actual[2] <= alpha,
                                          2 * min(left, right) <= alpha)
            }
        }
        testthat::expect_identical(serialize(scores, NULL), before)
    }
    testthat::expect_identical(.Random.seed, seed)
})

testthat::test_that("all-tied staged fallback preserves selected draws", {
    a <- rep(0, 4)
    w <- rep(1, 4)
    D <- matrix(0, 1, 4)
    draws <- rep(list(0L, 1L, 2L, 3L), 30)
    before <- serialize(draws, NULL)
    set.seed(238L)
    seed <- .Random.seed
    for (side in -1:1) {
        out <- environment$component_tail_record(0, a, w, D, draws, 1L, FALSE,
            10L, 100L, 10L, FALSE, side, 1, 0, TRUE)
        testthat::expect_equal(out$p, 1)
        testthat::expect_equal(out$resampling_dist, rep(0, 10))
    }
    testthat::expect_identical(serialize(draws, NULL), before)
    testthat::expect_identical(.Random.seed, seed)
})
