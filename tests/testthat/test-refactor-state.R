test_that("workflow transitions preserve sibling analyses and caller state", {
    state <- workflow_initial_state()
    original <- state
    for (step in c("set_analysis_parameters", "assign_grnas", "run_qc")) {
        state <- workflow_transition(state, step)
    }
    expect_identical(original, c(import_data = TRUE,
        set_analysis_parameters = FALSE, assign_grnas = FALSE, run_qc = FALSE,
        run_calibration_check = FALSE, run_power_check = FALSE,
        run_discovery_analysis = FALSE))
    ready <- state
    for (analysis in c("run_discovery_analysis", "run_calibration_check", "run_power_check")) {
        one <- workflow_transition(ready, analysis)
        expect_true(one[[analysis]])
        siblings <- setdiff(c("run_discovery_analysis", "run_calibration_check", "run_power_check"), analysis)
        expect_false(any(one[siblings]))
        state <- workflow_transition(state, analysis)
    }
    expect_true(all(state))
    expect_identical(workflow_status(state), state)
    expect_identical(workflow_transition(state, "run_qc"), ready)
    expect_identical(workflow_initial_state(), original)
    before <- serialize(original, NULL)
    expect_error(workflow_transition(original, "run_qc"),
                 "`assign_grnas\\(\\)` must be called")
    expect_identical(serialize(original, NULL), before)
})
