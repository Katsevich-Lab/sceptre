# Frozen behavior suite: snapshots of what every case outputs, so that a change
# meant to keep results the same can show that it does (see ../README.md).
# Each test-<group>.R file snapshots one group of cases into _snaps/<group>.md.
source(test_path("..", "harness.R"), local = TRUE)
input_dir <- test_path("inputs")
waive_internal <- identical(Sys.getenv("SCEPTRE_FROZEN_WAIVE_INTERNAL"), "true")

snapshot_cases <- function(group) {
    cases <- consistency_cases()
    for (name in names(cases)[vapply(cases, `[[`, "", "group") == group]) {
        local({
            case_name <- name
            observation <- consistency_run_case(case_name, input_dir)
            test_that(paste("public outputs:", case_name), {
                expect_snapshot_value(
                    consistency_round(consistency_public(observation)),
                    style = "deparse",
                    tolerance = 0
                )
            })
            # A restructuring that moves slots may waive these snapshots, and
            # must then say so in its pull request.
            test_that(paste("internal state:", case_name), {
                skip_if(waive_internal, "internal snapshots waived")
                expect_snapshot_value(
                    consistency_round(consistency_internal(observation)),
                    style = "deparse",
                    tolerance = 0
                )
            })
        })
    }
}
