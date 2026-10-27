# The fixtures are complete, unmodified and recorded by the current harness.

test_that("every fixture file has its recorded SHA256", {
    manifest <- consistency_manifest()
    directory <- consistency_fixture_dir()
    files <- names(manifest$files)
    present <- list.files(directory, recursive = TRUE)
    present <- setdiff(present, c("manifest.json", ".gitattributes"))
    expect_setequal(present, files)
    expect_identical(
        consistency_sha256(file.path(directory, files)),
        unname(unlist(manifest$files))
    )
})

test_that("fixtures were recorded by this harness", {
    manifest <- consistency_manifest()
    expect_identical(manifest$schema_version, consistency_schema_version)
    expect_identical(names(manifest$harness), suite_files)
    for (file in suite_files) {
        expect_identical(
            manifest$harness[[file]],
            consistency_sha256(file.path(suite_dir, file)),
            info = paste(
                file,
                "changed after the fixtures were recorded. Regenerate them",
                "from the reference commit (README.md) and review the diff."
            )
        )
    }
    expect_match(manifest$reference$commit, "^[0-9a-f]{40}$")
})

test_that("the manifest lists exactly the registered cases and known bugs", {
    manifest <- consistency_manifest()
    expect_identical(names(manifest$cases), names(consistency_cases()))
    expect_identical(
        unlist(manifest$known_bug_cases),
        consistency_known_bug_cases()
    )
    expect_true(all(consistency_known_bug_cases() %in% names(manifest$cases)))
    flagged <- names(Filter(function(x) !is.null(x$known_bug), manifest$cases))
    expect_setequal(flagged, consistency_known_bug_cases())
})

test_that("the saved-object input is the released 0.99.0 file", {
    path <- file.path(
        consistency_fixture_dir(),
        "saved-objects",
        consistency_saved_object_file
    )
    expect_identical(consistency_sha256(path), consistency_saved_object_sha256)
})

test_that("inputs read from the package are unchanged", {
    expected <- consistency_manifest()$inputs
    actual <- consistency_input_digests()
    expect_identical(actual$synthetic_md5, expected$synthetic_md5)
    expect_identical(actual$example_data_md5, expected$example_data_md5)
    expect_identical(actual$extdata_sha256, expected$extdata_sha256)
})
