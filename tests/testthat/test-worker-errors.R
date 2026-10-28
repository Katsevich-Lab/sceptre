test_that("worker result checks report the original error and missing replies", {
    failed <- try(stop("Cannot read response gene_X"), silent = TRUE)
    expect_error(
        check_worker_results(list(list(result = 1), failed)),
        "Parallel worker 2 failed: Cannot read response gene_X",
        fixed = TRUE
    )
    expect_error(
        check_worker_results(list(list(result = 1), NULL)),
        "Parallel worker 2 returned no result",
        fixed = TRUE
    )
    ok <- list(list(result = 1), list(), data.table::data.table())
    expect_identical(check_worker_results(ok), ok)
})

# Five responses, so that two parallel workers each receive some of them
worker_error_object <- function(case = "low_complement_perm") {
    inputs <- refactor_make_inputs()
    inputs$discovery_pairs <- data.frame(
        grna_target = "A",
        response_id = c("r1", "r5", "r6", "r7", "r8")
    )
    object <- suppressMessages(refactor_prepare(
        inputs,
        refactor_case_definitions()[[case]]
    ))
    expect_identical(object@n_ok_discovery_pairs, 5L)
    object
}

test_that("a failing parallel worker reports its original error", {
    skip_on_os("windows")
    parent <- Sys.getpid()
    local_mocked_bindings(perform_response_precomputation = function(...) {
        if (Sys.getpid() != parent) {
            stop("Controlled worker failure")
        }
        stop("The test expected this call in a worker process")
    })
    # permutation test (low MOI) and CRT (high MOI)
    for (case in c("low_complement_perm", "high_complement_crt")) {
        object <- worker_error_object(case)
        expect_error(
            suppressWarnings(suppressMessages(run_discovery_analysis(
                object,
                parallel = TRUE,
                n_processors = 2,
                print_progress = FALSE
            ))),
            "Parallel worker 1 failed: Controlled worker failure",
            fixed = TRUE
        )
    }
})

test_that("a failing CRT association worker reports its original error", {
    skip_on_os("windows")
    # with singleton gRNA groups there are four groups (a1, a2, b1, b2), so
    # the association step runs in two workers
    parent <- Sys.getpid()
    local_mocked_bindings(crt_glm_factored_out = function(...) {
        if (Sys.getpid() != parent) {
            stop("Controlled association failure")
        }
        stop("The test expected this call in a worker process")
    })
    inputs <- refactor_make_inputs()
    case <- refactor_case_definitions()$high_complement_crt
    case$integration <- "singleton"
    object <- suppressMessages(refactor_prepare(inputs, case))
    expect_error(
        suppressWarnings(suppressMessages(run_discovery_analysis(
            object,
            parallel = TRUE,
            n_processors = 2,
            print_progress = FALSE
        ))),
        "Parallel worker 1 failed: Controlled association failure",
        fixed = TRUE
    )
})

test_that("a failing mixture assignment worker reports its original error", {
    skip_on_os("windows")
    parent <- Sys.getpid()
    local_mocked_bindings(obtain_em_assignments = function(...) {
        if (Sys.getpid() != parent) {
            stop("Controlled mixture failure")
        }
        stop("The test expected this call in a worker process")
    })
    inputs <- refactor_make_inputs()
    object <- refactor_import(inputs, refactor_case_definitions()$high_perm)
    object <- set_analysis_parameters(
        object,
        discovery_pairs = inputs$discovery_pairs,
        resampling_mechanism = "permutations"
    )
    expect_error(
        suppressWarnings(suppressMessages(assign_grnas(
            object,
            method = "mixture",
            formula_object = ~ log(response_n_umis),
            parallel = TRUE,
            n_processors = 2,
            print_progress = FALSE
        ))),
        "Parallel worker 1 failed: Controlled mixture failure",
        fixed = TRUE
    )
})

test_that("a parallel worker that exits without replying cannot drop results", {
    skip_on_os("windows")
    object <- worker_error_object()
    parent <- Sys.getpid()
    original <- perform_response_precomputation
    flag <- file.path(withr::local_tempdir(), "killed")
    local_mocked_bindings(perform_response_precomputation = function(...) {
        # Terminate one test-created child process (the first to create flag)
        if (Sys.getpid() != parent && dir.create(flag, showWarnings = FALSE)) {
            tools::pskill(Sys.getpid(), signal = 9L)
        }
        original(...)
    })
    expect_error(
        suppressWarnings(suppressMessages(run_discovery_analysis(
            object,
            parallel = TRUE,
            n_processors = 2,
            print_progress = FALSE
        ))),
        "Parallel worker [12] returned no result"
    )
})

test_that("parallel and serial discovery results agree", {
    skip_on_os("windows")
    object <- worker_error_object()
    serial <- suppressMessages(run_discovery_analysis(
        object,
        print_progress = FALSE
    ))
    parallel <- suppressMessages(run_discovery_analysis(
        object,
        parallel = TRUE,
        n_processors = 2,
        print_progress = FALSE
    ))
    expect_equal(parallel@discovery_result, serial@discovery_result)
})
