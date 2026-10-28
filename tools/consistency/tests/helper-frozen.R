# Frozen behavior suite: snapshots of what every case outputs, so that a change
# meant to keep results the same can show that it does (see ../README.md).
# Each test-<group>.R file snapshots one group of cases into _snaps/<group>/.
source(test_path("..", "harness.R"), local = TRUE)
input_dir <- test_path("inputs")
waive_internal <- identical(Sys.getenv("SCEPTRE_FROZEN_WAIVE_INTERNAL"), "true")
tolerance <- as.numeric(Sys.getenv("SCEPTRE_FROZEN_TOLERANCE", "1e-8"))

# Snapshot `x` as R code, exact to the last bit and readable in a diff, and
# compare it with the recorded snapshot by consistency_differences(): testthat
# compares numbers relative to their average size, which would not notice a
# p-value moving from 1e-25 to 1e-14.
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

snapshot_cases <- function(group) {
    cases <- consistency_cases()
    for (name in names(cases)[vapply(cases, `[[`, "", "group") == group]) {
        local({
            case_name <- name
            observation <- consistency_run_case(case_name, input_dir)
            test_that(paste("public outputs:", case_name), {
                expect_frozen(consistency_public(observation), case_name)
            })
            # A restructuring that moves slots may waive these snapshots, and
            # must then say so in its pull request.
            test_that(paste("internal state:", case_name), {
                skip_if(waive_internal, "internal snapshots waived")
                expect_frozen(
                    consistency_internal(observation),
                    paste0(case_name, "-internal")
                )
            })
        })
    }
}
