test_that("backed save preserves imported descriptors and reopen verifies content", {
    skip_if_not_installed("ondisc")
    skip_if_not_installed("rhdf5")
    inputs <- readRDS(test_path("fixtures", "refactor", "backed-inputs.rds"))
    case <- readRDS(test_path("fixtures", "refactor", "cases.rds"))$high_perm
    directory <- tempfile("immutable-backed-")
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    object <- suppressMessages(refactor_prepare(inputs, case, backed = TRUE, directory = directory))
    expect_true(length(object@data@backend) > 0L) # anchored at import, before save
    object <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE))
    expect_gt(length(object@fits@response_precomputations), 0L)
    expect_true(result_record_is_current(object, "run_discovery_analysis"))
    write_ondisc_backed_sceptre_object(object, directory)
    original_file <- readBin(file.path(directory, "sceptre_object.rds"), "raw",
        n = file.info(file.path(directory, "sceptre_object.rds"))$size)
    reopened <- read_ondisc_backed_sceptre_object(file.path(directory, "sceptre_object.rds"),
        file.path(directory, "response.odm"), file.path(directory, "grna.odm"))
    pairs <- reopened@preparation@discovery_pairs_with_info
    pairs <- pairs[pairs$pass_qc, ]
    expect_gt(length(prepare_analysis_context(reopened, "discovery_analysis", pairs)$response_fits), 0L)
    for (modality in c("response", "grna")) {
        path <- file.path(directory, paste0(modality, ".odm"))
        original <- readBin(path, "raw", n = file.info(path)$size)
        counts <- rhdf5::h5read(path, "x")
        counts[1L] <- counts[1L] + 1000L
        rhdf5::h5write(counts, path, "x")
        # Active stores are immutable by contract. Routine getters/reporting
        # do not monitor unsupported external edits. Reopen and save still
        # guard the original persisted descriptor rather than rebaselining it.
        for (candidate in list(object, reopened)) {
            expect_error(write_ondisc_backed_sceptre_object(candidate, directory), "content does not match")
        }
        expect_error(read_ondisc_backed_sceptre_object(file.path(directory, "sceptre_object.rds"),
            file.path(directory, "response.odm"), file.path(directory, "grna.odm")), "content does not match")
        expect_identical(readBin(file.path(directory, "sceptre_object.rds"), "raw",
            n = length(original_file)), original_file)
        writeBin(original, path)
        expect_true(result_record_is_current(reopened, "run_discovery_analysis"))
    }
    # File location can change after exact content verification.
    moved <- tempfile("relocated-backed-")
    dir.create(moved)
    on.exit(unlink(moved, recursive = TRUE), add = TRUE)
    expect_true(all(file.copy(file.path(directory, c("response.odm", "grna.odm")), moved)))
    relocated <- read_ondisc_backed_sceptre_object(file.path(directory, "sceptre_object.rds"),
        file.path(moved, "response.odm"), file.path(moved, "grna.odm"))
    expect_identical(relocated@data@identity, reopened@data@identity)
    expect_true(result_record_is_current(relocated, "run_discovery_analysis"))
    expect_identical(relocated@fits@response_precomputations, reopened@fits@response_precomputations)
})
