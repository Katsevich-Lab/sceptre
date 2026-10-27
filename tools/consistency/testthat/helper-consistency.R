# Frozen behavior suite: test helpers. Run through tools/consistency/consistency.R.
suite_dir <- normalizePath(
    testthat::test_path(".."),
    winslash = "/",
    mustWork = TRUE
)
suite_files <- c("harness.R", "fixture-io.R", "export.R")
for (file in suite_files) {
    source(file.path(suite_dir, file), local = TRUE)
}

consistency_fixture_dir <- function() {
    directory <- Sys.getenv(
        "SCEPTRE_CONSISTENCY_FIXTURES",
        file.path(suite_dir, "fixtures")
    )
    normalizePath(directory, winslash = "/", mustWork = TRUE)
}

consistency_manifest <- local({
    cached <- NULL
    function() {
        if (is.null(cached)) {
            cached <<- jsonlite::read_json(
                file.path(consistency_fixture_dir(), "manifest.json"),
                simplifyVector = FALSE
            )
        }
        cached
    }
})

# strict on the recording platform (same R version and platform string),
# portable elsewhere, unless SCEPTRE_CONSISTENCY_TIER says otherwise.
consistency_tier <- function() {
    tier <- Sys.getenv("SCEPTRE_CONSISTENCY_TIER")
    if (nzchar(tier)) {
        return(match.arg(tier, names(consistency_tolerances)))
    }
    recorded <- consistency_manifest()$environment
    same <- identical(recorded$R, R.version.string) &&
        identical(recorded$platform, R.version$platform)
    if (same) "strict" else "portable"
}

consistency_internal_waived <- function() {
    identical(Sys.getenv("SCEPTRE_CONSISTENCY_WAIVE_INTERNAL"), "true")
}

consistency_expected <- function(name) {
    readRDS(file.path(consistency_fixture_dir(), "cases", paste0(name, ".rds")))
}

# Differences between two observations, as readable paths. Doubles agree when
# |actual - expected| <= tolerance * |expected| elementwise (an expected exact
# zero allows |actual| <= 1e-12), with the same NA/NaN/Inf positions. The
# comparison is relative at every magnitude, so a p-value of 1e-25 that moves
# to 1e-14 is a difference. Every other value must be identical.
consistency_differences <- function(
    actual,
    expected,
    tolerance = consistency_tolerances[[consistency_tier()]],
    path = "observation"
) {
    if (!identical(class(actual), class(expected))) {
        return(paste0(path, ": class differs"))
    }
    if (is.list(expected)) {
        if (!identical(names(actual), names(expected))) {
            return(paste0(path, ": names differ"))
        }
        if (length(actual) != length(expected)) {
            return(paste0(path, ": length differs"))
        }
        keys <- if (is.null(names(expected))) {
            paste0("[[", seq_along(expected), "]]")
        } else {
            paste0("$", names(expected))
        }
        return(unlist(lapply(seq_along(expected), function(k) {
            consistency_differences(
                actual[[k]],
                expected[[k]],
                tolerance,
                paste0(path, keys[[k]])
            )
        })))
    }
    if (!is.double(expected)) {
        if (identical(actual, expected)) {
            return(character())
        }
        return(paste0(path, ": value differs"))
    }
    if (!identical(attributes(actual), attributes(expected))) {
        return(paste0(path, ": attributes differ"))
    }
    special <- function(x) {
        ifelse(
            is.nan(x),
            "NaN",
            ifelse(
                is.na(x),
                "NA",
                ifelse(is.infinite(x), ifelse(x > 0, "Inf", "-Inf"), "finite")
            )
        )
    }
    if (!identical(special(actual), special(expected))) {
        return(paste0(path, ": NA/NaN/Inf pattern differs"))
    }
    finite <- is.finite(expected)
    a <- actual[finite]
    e <- expected[finite]
    gap <- abs(a - e)
    allowed <- ifelse(e == 0, 1e-12, tolerance * abs(e))
    bad <- gap > allowed
    if (any(bad)) {
        worst <- which.max(gap / allowed)
        return(sprintf(
            "%s: %d value(s) differ (worst: actual %.17g, expected %.17g)",
            path,
            sum(bad),
            a[worst],
            e[worst]
        ))
    }
    character()
}

consistency_expect_same <- function(actual, expected, name, part) {
    differences <- consistency_differences(actual, expected)
    bug <- consistency_known_bugs[[name]]
    note <- if (is.null(bug)) {
        paste0(part, " output differs from the frozen reference.")
    } else {
        paste0(
            part,
            " output differs. This case records known-bug behavior (",
            bug$fix_pr,
            ": ",
            bug$bugs,
            "); see README.md before regenerating."
        )
    }
    testthat::expect(
        length(differences) == 0L,
        paste0(
            name,
            " [",
            consistency_tier(),
            " tier]: ",
            note,
            "\n",
            paste(utils::head(differences, 20L), collapse = "\n"),
            if (length(differences) > 20L) {
                paste0("\n... and ", length(differences) - 20L, " more")
            }
        )
    )
    invisible(differences)
}

# Run every case once per test session; public and internal tests share it.
consistency_actual <- local({
    cache <- list()
    function(name) {
        if (is.null(cache[[name]])) {
            cache[[name]] <<- consistency_run_case(
                name,
                consistency_fixture_dir()
            )
        }
        cache[[name]]
    }
})
