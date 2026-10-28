# The inputs and the list of cases are part of what the snapshots record: a
# changed input would otherwise look like changed behavior, and a dropped case
# would silently stop being checked.
test_that("inputs and cases are unchanged", {
    extdata <- c(
        "highmoi_example/gem_group_1/barcodes.tsv.gz",
        "highmoi_example/gem_group_1/features.tsv.gz",
        "highmoi_example/gem_group_1/matrix.mtx.gz",
        "highmoi_example/gem_group_2/barcodes.tsv.gz",
        "highmoi_example/gem_group_2/features.tsv.gz",
        "highmoi_example/gem_group_2/matrix.mtx.gz",
        "parse_example/all_genes.csv",
        "parse_example/all_grnas.csv",
        "parse_example/gene_mat.mtx",
        "parse_example/grna_mat.mtx"
    )
    saved_object <- file.path(
        input_dir,
        "saved-objects",
        consistency_saved_object_file
    )
    expect_frozen(
        list(
            cases = vapply(consistency_cases(), `[[`, "", "group"),
            synthetic = consistency_synthetic_inputs(),
            example_data = lapply(consistency_example_data(), consistency_md5),
            extdata = stats::setNames(
                unname(tools::md5sum(vapply(extdata, consistency_extdata, ""))),
                extdata
            ),
            saved_object = unname(tools::md5sum(saved_object))
        ),
        "inputs"
    )
})

test_that("two workers give the serial results", {
    skip_on_os("windows")
    setting <- consistency_workflow_settings()$high_complement_crt
    run <- function(parallel) {
        observation <- consistency_seeded(function() {
            consistency_run_workflow(setting, parallel = parallel)
        })
        observation <- consistency_sort_precomputations(observation)
        if (waive_internal) consistency_public(observation) else observation
    }
    expect_equal(
        consistency_differences(run(parallel = TRUE), run(parallel = FALSE), 0),
        character()
    )
})

# Each group of cases has its test file, so deleting one cannot silently drop
# its cases.
test_that("every group of cases has a test file", {
    groups <- unique(vapply(consistency_cases(), `[[`, "", "group"))
    expect_setequal(
        paste0("test-", gsub("_", "-", groups), ".R"),
        setdiff(list.files(test_path(), "^test-.*[.]R$"), "test-inputs.R")
    )
})
