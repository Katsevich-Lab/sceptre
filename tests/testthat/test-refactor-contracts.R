# Permanent tests check formulas and compare the package's current execution
# modes. Original-versus-refactored comparisons run separately under tools/.
test_that("supported workflows produce valid complete results reproducibly", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    before <- serialize(inputs, NULL)
    for (name in names(cases)) {
        actual <- refactor_run_case(inputs, cases[[name]])
        repeat_run <- refactor_run_case(inputs, cases[[name]])
        expect_identical(actual, repeat_run, info = name)
        expect_workflow_results(actual, inputs, cases[[name]])
        final <- actual$value$final
        expect_equal(length(final$analysis_cells), ncol(inputs$response_matrix) - length(inputs$remove_cells))
        expect_identical(serialize(inputs, NULL), before)
    }
})

test_that("repeated analyses and cache invalidation preserve current fit views", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    actual <- refactor_cache_observations(inputs, cases$low_complement_perm)
    expect_equal(actual$value$warm$discovery_result, actual$value$cold$discovery_result, tolerance = 1e-8)
    expect_equal(actual$value$cleared$discovery_result, actual$value$cold$discovery_result, tolerance = 1e-8)
    expect_identical(actual$value$retained_side_cache, actual$value$cold$response_precomputations)
    expect_length(actual$value$formula_cache, 0L)
    expect_length(actual$value$changed_qc_cache, 0L)
    expect_false(identical(actual$value$reassigned$analysis_groups, actual$value$cold$analysis_groups))
    expect_equal(actual$value$reassigned$discovery_result, actual$value$reassigned_cold$discovery_result, tolerance = 1e-8)
})

test_that("serial and two-worker execution preserve pair identity and results", {
    skip_on_os("windows")
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    serial <- refactor_run_case(inputs, cases$high_complement_crt)
    parallel <- refactor_run_case(inputs, cases$high_complement_crt, parallel = TRUE)
    # Workers can append fits in different orders; compare them by gene name.
    for (stage in c("calibrated", "final")) {
        for (out in c("serial", "parallel")) {
            value <- get(out)
            fits <- value$value[[stage]]$response_precomputations
            value$value[[stage]]$response_precomputations <- fits[sort(names(fits))]
            assign(out, value)
        }
    }
    expect_equal(parallel, serial, tolerance = 1e-8)
})

test_that("backed permutation workflows and independent reattachment match sparse", {
    if (identical(Sys.getenv("SCEPTRE_REQUIRE_ONDISC"), "true")) {
        expect_true(requireNamespace("ondisc", quietly = TRUE))
    }
    skip_if_not_installed("ondisc")
    inputs <- readRDS(file.path(refactor_fixture_dir(), "backed-inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    directory <- tempfile("refactor-ondisc-")
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    actual <- refactor_run_case(inputs, cases$high_perm, backed = TRUE, directory = directory)
    expected <- refactor_run_case(inputs, cases$high_perm)
    actual$value$rng_after <- expected$value$rng_after <- NULL
    expect_identical(actual$warnings, character())
    expect_identical(expected$warnings, character())
    expect_equal(actual, expected, tolerance = 1e-8)

    object <- refactor_capture(refactor_prepare(inputs, cases$high_perm, TRUE, directory))$value
    sceptre::write_ondisc_backed_sceptre_object(object, directory)
    rm(object)
    gc()
    reopened <- sceptre::read_ondisc_backed_sceptre_object(
        file.path(directory, "sceptre_object.rds"), file.path(directory, "response.odm"), file.path(directory, "grna.odm"))
    reopened_counts <- do.call(rbind, lapply(rownames(inputs$response_matrix), function(id) {
        refactor_private("load_row")(sceptre::get_response_matrix(reopened), id)
    }))
    expect_equal(unname(reopened_counts), unname(inputs$response_matrix))
    continued <- refactor_capture(sceptre::run_discovery_analysis(reopened, print_progress = FALSE, output_amount = 3L))
    sparse <- refactor_capture(refactor_prepare(inputs, cases$high_perm))$value
    discovery_only <- refactor_capture(sceptre::run_discovery_analysis(
        sparse, print_progress = FALSE, output_amount = 3L))
    discovery_only$value <- refactor_observe(discovery_only$value)
    expect_equal(refactor_observe(continued$value)$discovery_result, discovery_only$value$discovery_result, tolerance = 1e-8)

    other <- tempfile("refactor-ondisc-other-")
    on.exit(unlink(other, recursive = TRUE), add = TRUE)
    invisible(refactor_capture(refactor_import(inputs, cases$high_perm, TRUE, other)))
    expect_error(sceptre::read_ondisc_backed_sceptre_object(
        file.path(directory, "sceptre_object.rds"), file.path(other, "response.odm"), file.path(directory, "grna.odm")), "distinct IDs")
    expect_error(sceptre::read_ondisc_backed_sceptre_object(
        file.path(directory, "sceptre_object.rds"), file.path(directory, "response.odm"), file.path(other, "grna.odm")), "distinct IDs")
})

test_that("current saved objects resume without changing discovery results", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_complement_perm
    prepared <- refactor_capture(refactor_prepare(inputs, case))$value
    file <- tempfile(fileext = ".rds")
    on.exit(unlink(file), add = TRUE)
    saveRDS(prepared, file)
    restored <- readRDS(file)
    before <- serialize(restored, NULL)
    first <- refactor_capture(sceptre::run_discovery_analysis(prepared, print_progress = FALSE))
    second <- refactor_capture(sceptre::run_discovery_analysis(restored, print_progress = FALSE))
    expect_identical(refactor_observe(first$value), refactor_observe(second$value))
    expect_identical(serialize(restored, NULL), before)
})

test_that("a versioned saved object remains a supported input", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_complement_perm
    restored <- refactor_private("migrate_sceptre_object")(readRDS(
        test_path("fixtures", "saved-objects", "sceptre-0.99.0-21f9ea0.rds")))
    before <- serialize(restored, NULL)
    current <- refactor_capture(refactor_prepare(inputs, case))
    old_run <- refactor_capture(sceptre::run_discovery_analysis(restored, print_progress = FALSE))
    new_run <- refactor_capture(sceptre::run_discovery_analysis(current$value, print_progress = FALSE))
    expect_identical(current$warnings, character())
    expect_identical(old_run$warnings, character())
    expect_identical(new_run$warnings, character())
    expect_identical(refactor_observe(old_run$value), refactor_observe(new_run$value))
    expect_identical(serialize(restored, NULL), before)
})
