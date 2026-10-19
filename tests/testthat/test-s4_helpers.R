# `perform_status_check_and_update` is used in each of the main sceptre_object
# functions aside from `import_data`
test_that("perform_status_check_and_update", {
    # `import_data` initializes @functs_called with the "import_data" value
    # at TRUE so we'll never see an all-FALSE @functs_called
    initial_functs_called <- c(
        import_data = TRUE,
        set_analysis_parameters = FALSE,
        assign_grnas = FALSE,
        run_qc = FALSE,
        run_calibration_check = FALSE,
        run_power_check = FALSE,
        run_discovery_analysis = FALSE
    )

    blank_object <- new("sceptre_object")
    blank_object@functs_called <- initial_functs_called

    ## call `set_analysis_parameters`
    results <- perform_status_check_and_update(
        blank_object,
        "set_analysis_parameters"
    )
    expect_identical(results@last_function_called, blank_object@last_function_called)
    expect_equal(
        results@functs_called |> as.logical(),
        rep(c(TRUE, FALSE), c(1, 6))
    )

    ## try to use a function before calling the downstream ones
    expect_error(
        perform_status_check_and_update(blank_object, "run_qc"),
        regex = "The function `assign_grnas\\(\\)` must be called before `run_qc\\(\\)` is called"
    )

    ## Beginning a stage is insufficient to unlock a dependent stage.
    expect_error(perform_status_check_and_update(results, "assign_grnas"),
        "set_analysis_parameters")
    # Model completed computation explicitly; public-stage success paths are
    # exercised by the workflow tests and component lifecycle regressions.
    results <- blank_object
    for (stage in c("set_analysis_parameters", "assign_grnas", "run_qc")) {
        results <- perform_status_check_and_update(results, stage)
        expect_false(results@functs_called[[stage]])
        results@functs_called[stage] <- TRUE
        results@last_function_called <- stage
    }
    expect_equal(results@last_function_called, "run_qc")
    expect_equal(results@functs_called |> as.logical(), rep(c(TRUE, FALSE), c(4, 3)))

    ## confirming if we rerun a lower level one after having run a higher level one
    results <- results |> perform_status_check_and_update("assign_grnas")
    expect_equal(results@last_function_called, "run_qc")
    expect_false(results@functs_called["assign_grnas"])
    expect_false(results@functs_called["run_qc"]) # this was TRUE before
    expect_true(results@functs_called["set_analysis_parameters"]) # this should still be TRUE
})
