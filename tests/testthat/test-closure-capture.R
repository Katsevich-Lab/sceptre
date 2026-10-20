test_that("treatment lookups capture their chosen context at construction", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$high_complement_crt
    object <- refactor_capture(refactor_prepare(inputs, case))$value
    pairs <- object@preparation@discovery_pairs_with_info
    context <- refactor_private("prepare_analysis_context")(object,
        "discovery_analysis", pairs[pairs$pass_qc, , drop = FALSE])
    index <- context$index
    group <- pairs$grna_group[which(pairs$pass_qc)[1L]]
    expected <- refactor_private("group_treatment_analysis")(index, group)
    lookup <- refactor_private("make_treatment_index_lookup")(index)
    index <- NULL
    expect_identical(lookup(group), list(trt_idxs = expected, n_trt = length(expected)))
})
