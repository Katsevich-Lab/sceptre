# The comparison rejects small, deliberate changes to a recorded observation.
# These tests use the strict tolerance whatever the platform's tier.
strict <- consistency_tolerances[["strict"]]

perturb_table <- function(observation, table, f) {
    plain <- observation$value$tables[[table]]
    plain$columns <- f(plain$columns)
    observation$value$tables[[table]] <- plain
    observation
}

test_that("the comparison detects perturbed numbers, decisions and rows", {
    expected <- consistency_expected("workflow_low_complement_crt")
    expect_length(consistency_differences(expected, expected, strict), 0L)
    detects <- function(f) {
        changed <- perturb_table(expected, "discovery_result", f)
        length(consistency_differences(changed, expected, strict)) > 0L
    }
    expect_true(detects(function(x) {
        x$p_value[1L] <- x$p_value[1L] * (1 + 1e-6)
        x
    }))
    expect_true(detects(function(x) {
        x$significant[1L] <- !x$significant[1L]
        x
    }))
    expect_true(detects(function(x) lapply(x, function(v) v[-1L])))
    expect_true(detects(function(x) lapply(x, rev)))
    expect_true(detects(function(x) {
        x$response_id[1L] <- "r8"
        x
    }))
    changed <- expected
    changed$warnings <- "Unexpected warning"
    expect_gt(length(consistency_differences(changed, expected, strict)), 0L)
    changed <- expected
    changed$rng_after[3L] <- changed$rng_after[3L] + 1L
    expect_gt(length(consistency_differences(changed, expected, strict)), 0L)
    # Differences within the strict tolerance are accepted.
    expect_false(detects(function(x) {
        x$p_value <- x$p_value * (1 + 1e-10)
        x
    }))
})

test_that("tiny p-values are compared relatively, not against a floor", {
    expected <- consistency_expected("example_high_defaults")
    p <- expected$value$tables$power_result$columns$p_value
    tiny <- which(p < 1e-12)
    expect_gt(length(tiny), 0L)
    for (factor in c(10, 1 + 1e-6)) {
        changed <- perturb_table(expected, "power_result", function(x) {
            x$p_value[tiny[[1L]]] <- x$p_value[tiny[[1L]]] * factor
            x
        })
        expect_gt(
            length(consistency_differences(changed, expected, strict)),
            0L
        )
    }
})

test_that("factor levels, classes and row names are compared", {
    expected <- consistency_expected("workflow_low_complement_crt")
    covariates <- expected$value$tables$cell_covariates
    expect_true(is.factor(covariates$columns$batch))
    changed <- expected
    batch <- covariates$columns$batch
    changed$value$tables$cell_covariates$columns$batch <- factor(
        as.character(batch),
        levels = rev(levels(batch))
    )
    expect_gt(length(consistency_differences(changed, expected, strict)), 0L)
    changed <- expected
    changed$value$tables$cell_covariates$columns$batch <- as.character(batch)
    expect_gt(length(consistency_differences(changed, expected, strict)), 0L)
    changed <- expected
    changed$value$tables$discovery_result$class <- c(
        "modified",
        expected$value$tables$discovery_result$class
    )
    expect_gt(length(consistency_differences(changed, expected, strict)), 0L)
})
