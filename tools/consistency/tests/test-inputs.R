# The inputs are part of what the snapshots record: a changed input would
# otherwise look like changed behavior.
test_that("inputs are unchanged", {
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
    expect_snapshot_value(
        list(
            synthetic = consistency_round(consistency_synthetic_inputs()),
            example_data = lapply(consistency_example_data(), consistency_md5),
            extdata = stats::setNames(
                unname(tools::md5sum(vapply(extdata, consistency_extdata, ""))),
                extdata
            )
        ),
        style = "deparse",
        tolerance = 0
    )
})

test_that("two workers give the serial results", {
    skip_on_os("windows")
    setting <- consistency_workflow_settings()$high_complement_crt
    run <- function(parallel) {
        observation <- consistency_seeded(function() {
            consistency_run_workflow(setting, parallel = parallel)
        })
        observation <- consistency_round(
            consistency_sort_precomputations(observation)
        )
        if (waive_internal) consistency_public(observation) else observation
    }
    expect_identical(run(parallel = TRUE), run(parallel = FALSE))
})
