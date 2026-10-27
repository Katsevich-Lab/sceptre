test_that("serial and two-worker runs give the same workflow results", {
    skip_on_os("windows")
    name <- "workflow_high_complement_crt"
    setting <- consistency_workflow_settings()$high_complement_crt
    actual <- consistency_seeded(function() {
        consistency_run_workflow(setting, parallel = TRUE)
    })
    actual <- consistency_sort_precomputations(actual)
    expected <- consistency_sort_precomputations(consistency_expected(name))
    consistency_expect_same(
        consistency_public(actual),
        consistency_public(expected),
        name,
        "Public (two workers)"
    )
    if (!consistency_internal_waived()) {
        consistency_expect_same(
            consistency_internal(actual),
            consistency_internal(expected),
            name,
            "Internal (two workers)"
        )
    }
})
