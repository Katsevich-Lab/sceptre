test_that("reference provenance and every fixture digest are intact", {
    manifest <- readRDS(file.path(refactor_fixture_dir(), "manifest.rds"))
    expect_identical(manifest$source_commit, "21f9ea098b69c444a884a49bc254de711968e3b5")
    expect_true(manifest$reference_source_verified)
    expect_true(manifest$clean_native_rebuild)
    files <- file.path(refactor_fixture_dir(), names(manifest$artifacts_md5))
    expect_true(all(file.exists(files)))
    expect_identical(unname(tools::md5sum(files)), unname(manifest$artifacts_md5))
})

test_that("refactor preserves the exact public API and existing S4 facade", {
    expected <- readRDS(file.path(refactor_fixture_dir(), "api.rds"))
    expect_identical(refactor_api(), expected)
})

test_that("supported scientific workflows reproduce the pinned baseline", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    expected <- readRDS(file.path(refactor_fixture_dir(), "workflows.rds"))
    for (name in names(cases)) {
        actual <- refactor_run_case(inputs, cases[[name]])
        # Fitted/floating fields allow the declared cross-toolchain tolerance.
        expect_equal(actual, expected[[name]], tolerance = 1e-8, info = name)
        for (stage in c("prepared", "calibrated", "final")) {
            a <- actual$value[[stage]]
            b <- expected[[name]]$value[[stage]]
            for (field in c("cells_in_use", "grna_assignments", "B1", "B2", "B3", "n_ok_discovery_pairs")) {
                expect_identical(a[[field]], b[[field]], info = paste(name, stage, field))
            }
            for (table in c("discovery_pairs_with_info", "positive_control_pairs_with_info", "negative_control_pairs", "calibration_result", "power_result", "discovery_result")) {
                for (column in intersect(c("grna_target", "grna_group", "response_id", "pass_qc", "significant", "stage"), names(a[[table]]$columns))) {
                    expect_identical(a[[table]]$columns[[column]], b[[table]]$columns[[column]], info = paste(name, stage, table, column))
                }
            }
        }
        if (identical(name, "low_bonferroni")) {
            info <- actual$value$prepared$discovery_pairs_with_info$columns
            selected <- info$response_id == "r3" & info$grna_group %in% c("a1", "a2")
            expect_identical(sort(info$pass_qc[selected]), c(FALSE, TRUE))
        }
    }
})

test_that("ordered assignment remapping preserves source and fit views", {
    actual <- refactor_index_observations()
    expect_identical(actual, readRDS(file.path(refactor_fixture_dir(), "indexing.rds")))
    expect_identical(actual$treatment_source, c(4L, 5L))
    expect_identical(actual$nt_source, c(8L, 2L, 11L))
    expect_identical(actual$combined_source, c(4L, 5L, 8L, 2L, 11L))
    expect_true(actual$caller_unchanged)
})

test_that("repeated analyses and existing cache invalidation match baseline", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    actual <- refactor_cache_observations(inputs, cases$low_complement_perm)
    expected <- readRDS(file.path(refactor_fixture_dir(), "cache.rds"))
    expect_equal(actual, expected, tolerance = 1e-8)
    expect_equal(actual$value$warm$discovery_result, actual$value$cold$discovery_result, tolerance = 1e-8)
    expect_equal(actual$value$cleared$discovery_result, actual$value$cold$discovery_result, tolerance = 1e-8)
    expect_identical(actual$value$retained_side_cache, actual$value$cold$response_precomputations)
    expect_length(actual$value$formula_cache, 0L)
    expect_length(actual$value$changed_qc_cache, 0L)
    expect_false(identical(actual$value$reassigned$grna_assignments, actual$value$cold$grna_assignments))
    expect_equal(actual$value$reassigned$discovery_result, actual$value$reassigned_cold$discovery_result, tolerance = 1e-8)
})

test_that("baseline saved objects remain analyzable", {
    object <- readRDS(file.path(refactor_fixture_dir(), "saved-object.rds"))
    result <- refactor_capture(sceptre::run_discovery_analysis(object, print_progress = FALSE, output_amount = 2L))
    actual <- list(value = refactor_observe(result$value), warnings = result$warnings)
    expect_equal(actual, readRDS(file.path(refactor_fixture_dir(), "saved-object-result.rds")), tolerance = 1e-8)
})

test_that("rejected requests and degenerate baseline behavior remain explicit", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_complement_perm
    actual <- refactor_edge_observations(inputs, case)
    expect_equal(actual, readRDS(file.path(refactor_fixture_dir(), "edges.rds")), tolerance = 1e-8)
    expect_match(actual$duplicate_requests$error, "duplicate rows", fixed = TRUE)
})

test_that("serial and two-worker execution preserve pair identity and results", {
    skip_on_os("windows")
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    cases <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))
    serial <- refactor_run_case(inputs, cases$high_complement_crt)
    parallel <- refactor_run_case(inputs, cases$high_complement_crt, parallel = TRUE)
    # Existing caches may have worker-specific insertion order; compare keyed fits.
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
    expected <- readRDS(file.path(refactor_fixture_dir(), "backed-workflow.rds"))
    # Backed import generates backing-file IDs using R's RNG. Scientific parity
    # does not require that import to consume the same draws as sparse import.
    backend_metadata <- readRDS(file.path(refactor_fixture_dir(), "backed-metadata.rds"))
    if (identical(as.character(packageVersion("ondisc")), backend_metadata$ondisc_version)) {
        expect_identical(actual$warnings, backend_metadata$warnings)
        expect_identical(actual$value$rng_after, backend_metadata$rng_after)
    }
    actual$value$rng_after <- expected$value$rng_after <- NULL
    actual$warnings <- expected$warnings <- NULL
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
    discovery_only <- readRDS(file.path(refactor_fixture_dir(), "backed-discovery.rds"))
    expect_equal(refactor_observe(continued$value)$discovery_result, discovery_only$value$discovery_result, tolerance = 1e-8)

    other <- tempfile("refactor-ondisc-other-")
    on.exit(unlink(other, recursive = TRUE), add = TRUE)
    invisible(refactor_capture(refactor_import(inputs, cases$high_perm, TRUE, other)))
    expect_error(sceptre::read_ondisc_backed_sceptre_object(
        file.path(directory, "sceptre_object.rds"), file.path(other, "response.odm"), file.path(directory, "grna.odm")), "distinct IDs")
    expect_error(sceptre::read_ondisc_backed_sceptre_object(
        file.path(directory, "sceptre_object.rds"), file.path(directory, "response.odm"), file.path(other, "grna.odm")), "distinct IDs")
})
