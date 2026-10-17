#!/usr/bin/env Rscript
# Direct tests for private native components, run from a source checkout. The
# installed-package suite tests stages, effect mathematics and draw laws.
# No historical implementation or expected output is loaded here.
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
    # Segment offsets are zero-based; cell indices in these banks are already
    # zero-based. Fixed-cardinality hybrid draws use only the requested prefix.
    fixed <- environment$component_segment(a, w, D, draws, 1L, 2L, 2L, FALSE)
    variable <- environment$component_segment(a, w, D, draws, 1L, 2L, 2L, TRUE)
    testthat::expect_equal(fixed, vapply(draws[2:3], function(x) score(x[1:2] + 1L), numeric(1)))
    testthat::expect_equal(variable, vapply(draws[2:3], function(x) score(x + 1L), numeric(1)))
    testthat::expect_length(environment$component_segment(a, w, D, draws, 2L, 0L, 2L, TRUE), 0L)
    testthat::expect_identical(serialize(list(a, w, D, treatment, draws), NULL), before)
})

testthat::test_that("empirical tails retain ties and the existing two-sided convention", {
    scores <- c(-1, 0, 0, 2)
    testthat::expect_equal(environment$component_empirical(scores, 0, -1L), 4/5)
    testthat::expect_equal(environment$component_empirical(scores, 0, 1L), 4/5)
    testthat::expect_equal(environment$component_empirical(scores, 0, 0L), 8/5)
})
