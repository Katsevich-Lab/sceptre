test_that("installed native entry points expose the coherent final private names only", {
    namespace <- asNamespace("sceptre")
    retired <- c("run_low_level_test_full_v4", "crt_index_sampler_fast", "crt_index_sampler",
                 "sample_combinations_v2", "sample_undercover_pairs_v2", "increment_matrix",
                 "increment_vector", "construct_negative_control_pairs_v2")
    retained <- c("run_association_test", "sample_crt_indices", "sample_calibration_groups",
                  "sample_calibration_pairs", "increment_group_indices_in_place")
    for (name in retired) expect_false(exists(name, namespace, inherits = FALSE), info = name)
    for (name in c(retained, "construct_negative_control_pairs"))
        expect_true(is.function(getFromNamespace(name, "sceptre")), info = name)
    registrations <- names(getDLLRegisteredRoutines(getLoadedDLLs()[["sceptre"]])[[".Call"]])
    expect_true(all(paste0("_sceptre_", retained) %in% registrations))
    expect_false(any(paste0("_sceptre_", retired) %in% registrations))
})

test_that("calibration coordinate increment mutates only explicitly owned scratch", {
    source <- matrix(c(0L, 2L, 1L, 0L), 2L)
    before <- serialize(source, NULL)
    owned <- unserialize(before)
    refactor_private("increment_group_indices_in_place")(owned)
    expect_identical(owned, source + 1L)
    expect_identical(serialize(source, NULL), before)
})
