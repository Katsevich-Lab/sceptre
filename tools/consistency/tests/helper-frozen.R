# Frozen behavior suite: snapshots of what every case outputs, so that a change
# meant to keep results the same can show that it does (see ../README.md).
# Each test-<group>.R file snapshots one group of cases into _snaps/<group>/.
source(test_path("..", "harness.R"), local = TRUE)
input_dir <- test_path("inputs")
waive_internal <- identical(Sys.getenv("SCEPTRE_FROZEN_WAIVE_INTERNAL"), "true")
tolerance <- as.numeric(Sys.getenv("SCEPTRE_FROZEN_TOLERANCE", "1e-8"))

# Snapshot `x` as R code with 17 significant digits, readable in a diff (R's
# parser can read a double back one unit in the last place off, far below any
# tolerance). consistency_differences() compares it with the recorded
# snapshot: testthat compares numbers relative to their average size, which
# would not notice a p-value moving from 1e-25 to 1e-14.
expect_frozen <- function(x, name) {
    path <- tempfile(fileext = ".R")
    writeLines(
        deparse(
            x,
            control = c(
                "keepNA",
                "keepInteger",
                "niceNames",
                "showAttributes",
                "digits17"
            )
        ),
        path
    )
    expect_snapshot_file(path, paste0(name, ".R"), compare = compare_frozen)
}

compare_frozen <- function(old, new) {
    read <- function(file) eval(parse(file, keep.source = FALSE), baseenv())
    differences <- consistency_differences(read(new), read(old), tolerance)
    if (length(differences)) {
        message(paste(utils::head(differences, 20L), collapse = "\n"))
    }
    length(differences) == 0L
}

# Public snapshots that differ by platform for a documented reason
# (../README.md, "Platforms"), skipped on that platform: case = c(os, reason).
platform_differences <- list(
    edge_intercept_only_nt = c(
        "linux",
        "permutation statistics tie the observed one; tie counts differ by platform"
    )
)

snapshot_cases <- function(group) {
    cases <- consistency_cases()
    for (name in names(cases)[vapply(cases, `[[`, "", "group") == group]) {
        local({
            case_name <- name
            observation <- consistency_run_case(case_name, input_dir)
            test_that(paste("public outputs:", case_name), {
                # announced first, so a skipped run does not delete it as unused
                announce_snapshot_file(paste0(case_name, ".R"))
                difference <- platform_differences[[case_name]]
                skip_if(
                    !is.null(difference) &&
                        tolower(Sys.info()[["sysname"]]) == difference[[1L]],
                    difference[[2L]]
                )
                expect_frozen(consistency_public(observation), case_name)
            })
            # A restructuring that moves slots may waive these snapshots, and
            # must then say so in its pull request.
            test_that(paste("internal state:", case_name), {
                # announced first, so a waived run does not delete it as unused
                announce_snapshot_file(paste0(case_name, "-internal.R"))
                skip_if(waive_internal, "internal snapshots waived")
                expect_frozen(
                    consistency_internal(observation),
                    paste0(case_name, "-internal")
                )
            })
        })
    }
}
