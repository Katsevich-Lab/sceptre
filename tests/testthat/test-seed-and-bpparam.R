seed_test_object <- function(case = "low_nt_perm") {
    refactor_import(refactor_make_inputs(), refactor_case_definitions()[[case]])
}

seed_test_analysis <- function(object, ...) {
    object |>
        sceptre::assign_grnas(
            method = "thresholding",
            threshold = 40,
            print_progress = FALSE,
            ...
        ) |>
        sceptre::run_qc(
            response_n_umis_range = c(0, 1),
            response_n_nonzero_range = c(0, 1),
            n_nonzero_trt_thresh = 7L,
            n_nonzero_cntrl_thresh = 7L
        ) |>
        sceptre::run_calibration_check(
            n_calibration_pairs = 4L,
            calibration_group_size = 2L,
            print_progress = FALSE,
            ...
        ) |>
        sceptre::run_discovery_analysis(print_progress = FALSE, ...)
}

test_that("the seed is checked, recorded and printed", {
    object <- seed_test_object()
    for (bad in list(-1, 1.5, NA, c(1, 2), "4", .Machine$integer.max + 1)) {
        expect_error(
            set_analysis_parameters(object, seed = bad),
            "`seed` must be NULL or a single whole number",
            fixed = TRUE
        )
    }
    object <- set_analysis_parameters(object, seed = 4)
    expect_identical(object@seed, 4L)
    expect_output(
        print(object),
        "Seed: 4",
        fixed = TRUE
    )
})

test_that("a NULL seed is drawn from R's random number generator", {
    object <- seed_test_object()
    set.seed(11)
    first <- set_analysis_parameters(object)
    set.seed(11)
    second <- set_analysis_parameters(object)
    expect_identical(first@seed, second@seed)
    expect_identical(
        first@seed,
        withr::with_seed(11, sample.int(.Machine$integer.max, 1L))
    )
    third <- withr::with_seed(12, set_analysis_parameters(object))
    expect_false(identical(first@seed, third@seed))
})

test_that("set.seed() reproduces an analysis run with the default seed", {
    pairs <- refactor_make_inputs()$discovery_pairs
    run <- function(k) {
        set.seed(k)
        suppressMessages(
            seed_test_object() |>
                set_analysis_parameters(discovery_pairs = pairs) |>
                seed_test_analysis()
        )
    }
    first <- run(5)
    second <- run(5)
    expect_identical(first@seed, second@seed)
    expect_identical(first@negative_control_pairs, second@negative_control_pairs)
    expect_identical(
        as.data.frame(first@calibration_result),
        as.data.frame(second@calibration_result)
    )
    expect_identical(
        as.data.frame(first@discovery_result),
        as.data.frame(second@discovery_result)
    )
})

test_that("BPPARAM is checked and selects the back end", {
    expect_error(
        get_bpparam(2),
        "`BPPARAM` should be a `BiocParallelParam` object",
        fixed = TRUE
    )
    expect_null(get_bpparam(NULL))
    expect_null(get_bpparam(BiocParallel::SerialParam()))
    snow <- BiocParallel::SnowParam(2)
    expect_identical(get_bpparam(snow), snow)
    # parallel = TRUE is an error on Windows (it forks), checked in
    # test-parallel_windows.R
    skip_on_os("windows")
    expect_warning(
        expect_error(
            get_bpparam(BiocParallel::SerialParam(), TRUE, 1),
            "`n_processors` should be set to the string 'auto'",
            fixed = TRUE
        ),
        "`parallel` and `n_processors` are deprecated"
    )
    expect_warning(
        bpparam <- get_bpparam(BiocParallel::SerialParam(), TRUE, 2),
        "`parallel` and `n_processors` are deprecated"
    )
    expect_s4_class(bpparam, "MulticoreParam")
    expect_identical(BiocParallel::bpnworkers(bpparam), 2L)
    # a BPPARAM with several workers takes precedence
    expect_warning(
        expect_identical(get_bpparam(snow, TRUE, 4), snow),
        "deprecated"
    )
})

test_that("the deprecated parallel = TRUE gives the results of BPPARAM", {
    skip_on_os("windows")
    object <- refactor_import(
        refactor_make_inputs(),
        refactor_case_definitions()[["low_complement_perm"]]
    ) |>
        set_analysis_parameters(
            discovery_pairs = refactor_make_inputs()$discovery_pairs,
            control_group = "complement",
            resampling_mechanism = "permutations",
            seed = 4L
        ) |>
        assign_grnas(
            method = "thresholding",
            threshold = 40,
            print_progress = FALSE
        ) |>
        run_qc(
            response_n_umis_range = c(0, 1),
            response_n_nonzero_range = c(0, 1),
            n_nonzero_trt_thresh = 7L,
            n_nonzero_cntrl_thresh = 7L
        )
    expect_warning(
        deprecated <- suppressMessages(run_discovery_analysis(
            object,
            print_progress = FALSE,
            parallel = TRUE,
            n_processors = 2
        )),
        "`parallel` and `n_processors` are deprecated"
    )
    current <- suppressMessages(run_discovery_analysis(
        object,
        print_progress = FALSE,
        BPPARAM = BiocParallel::MulticoreParam(2)
    ))
    expect_identical(deprecated@discovery_result, current@discovery_result)
})

test_that("parallel results equal serial results", {
    skip_on_os("windows")
    for (case in c("low_complement_perm", "high_complement_crt")) {
        object <- refactor_import(
            refactor_make_inputs(),
            refactor_case_definitions()[[case]]
        ) |>
            set_analysis_parameters(
                discovery_pairs = refactor_make_inputs()$discovery_pairs,
                control_group = refactor_case_definitions()[[case]]$control,
                resampling_mechanism = refactor_case_definitions()[[case]]$mechanism,
                seed = 4L
            )
        serial <- suppressWarnings(suppressMessages(
            seed_test_analysis(object)
        ))
        parallel <- suppressWarnings(suppressMessages(seed_test_analysis(
            object,
            BPPARAM = BiocParallel::MulticoreParam(2)
        )))
        for (slot_name in c(
            "initial_grna_assignment_list",
            "negative_control_pairs",
            "calibration_result",
            "discovery_result",
            "response_precomputations"
        )) {
            expect_identical(
                methods::slot(parallel, slot_name),
                methods::slot(serial, slot_name),
                info = paste(case, slot_name)
            )
        }
    }
})

test_that("results on SnowParam workers equal serial results", {
    # SnowParam is the parallel back end on Windows; its workers are fresh R
    # sessions that load the installed package
    skip_on_cran()
    object <- refactor_import(
        refactor_make_inputs(),
        refactor_case_definitions()[["low_complement_perm"]]
    ) |>
        set_analysis_parameters(
            discovery_pairs = refactor_make_inputs()$discovery_pairs,
            control_group = "complement",
            resampling_mechanism = "permutations",
            seed = 4L
        )
    serial <- suppressWarnings(suppressMessages(
        seed_test_analysis(object)
    ))
    snow <- suppressWarnings(suppressMessages(seed_test_analysis(
        object,
        BPPARAM = BiocParallel::SnowParam(2)
    )))
    for (slot_name in c(
        "negative_control_pairs",
        "calibration_result",
        "discovery_result"
    )) {
        expect_identical(
            methods::slot(snow, slot_name),
            methods::slot(serial, slot_name),
            info = slot_name
        )
    }
})
