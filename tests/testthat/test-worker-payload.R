# A parallel back end that serializes the function it runs (SnowParam, the
# only back end on Windows, or a started cluster) sends that function and its
# environment to every worker. The worker functions used to carry the whole
# frame of the analysis, including several copies of the response matrix and
# the sceptre_object; they must carry only what they use, and of in-memory
# response counts only the rows they read. Here bplapply_checked() is replaced
# by a serial loop that records the size of the function it is given.

payload_test_object <- function(n_extra_responses = 3000L) {
    data(lowmoi_example_data, package = "sceptre", envir = environment())
    d <- lowmoi_example_data
    set.seed(1)
    n_cells <- ncol(d$response_matrix)
    extra <- Matrix::rsparsematrix(
        n_extra_responses,
        n_cells,
        density = 0.3,
        rand.x = function(n) stats::rpois(n, 2) + 1
    )
    rownames(extra) <- paste0("extra_response_", seq_len(n_extra_responses))
    colnames(extra) <- colnames(d$response_matrix)
    response_matrix <- rbind(
        methods::as(d$response_matrix, "CsparseMatrix"),
        methods::as(extra, "CsparseMatrix")
    )
    import_data(
        response_matrix = response_matrix,
        grna_matrix = d$grna_matrix,
        grna_target_data_frame = d$grna_target_data_frame,
        extra_covariates = d$extra_covariates,
        moi = "low"
    )
}


test_that("parallel workers receive only the objects and rows they use", {
    sceptre_object <- payload_test_object()
    response_size <- length(serialize(get_response_matrix(sceptre_object), NULL))
    sizes <- numeric()
    local_mocked_bindings(
        bplapply_checked = function(X, FUN, BPPARAM) {
            sizes <<- c(sizes, length(serialize(FUN, NULL)))
            lapply(X, FUN)
        }
    )
    analyze <- function(mechanism, control_group, BPPARAM) {
        suppressMessages(
            sceptre_object |>
                set_analysis_parameters(
                    discovery_pairs = construct_positive_control_pairs(sceptre_object),
                    resampling_mechanism = mechanism,
                    control_group = control_group,
                    seed = 4
                ) |>
                assign_grnas(print_progress = FALSE) |>
                run_qc() |>
                run_discovery_analysis(print_progress = FALSE, BPPARAM = BPPARAM)
        )
    }
    # the complement control group runs the response regressions in parallel
    # too (a separate worker function under the CRT)
    for (mechanism in c("permutations", "crt")) {
        for (control_group in c("nt_cells", "complement")) {
            sizes <- numeric()
            # SnowParam(2) is never started: bplapply_checked() is replaced
            parallel <- analyze(
                mechanism,
                control_group,
                BiocParallel::SnowParam(2)
            )
            expect_gt(length(sizes), 0L)
            # the discovery pairs involve a few of the 3,100 responses
            expect_lt(max(sizes), 0.1 * response_size)
            expect_equal(
                get_result(parallel, "run_discovery_analysis"),
                get_result(
                    analyze(mechanism, control_group, BiocParallel::SerialParam()),
                    "run_discovery_analysis"
                )
            )
        }
    }
})


test_that("mixture assignment workers do not receive the cell covariates", {
    data(highmoi_example_data, package = "sceptre", envir = environment())
    data(grna_target_data_frame_highmoi, package = "sceptre", envir = environment())
    sceptre_object <- import_data(
        response_matrix = highmoi_example_data$response_matrix,
        grna_matrix = highmoi_example_data$grna_matrix,
        grna_target_data_frame = grna_target_data_frame_highmoi,
        moi = "high",
        extra_covariates = highmoi_example_data$extra_covariates
    )
    shipped <- list()
    local_mocked_bindings(
        bplapply_checked = function(X, FUN, BPPARAM) {
            shipped[[length(shipped) + 1L]] <<- FUN
            lapply(X, FUN)
        }
    )
    parallel <- suppressMessages(
        sceptre_object |>
            set_analysis_parameters(
                discovery_pairs = construct_positive_control_pairs(sceptre_object),
                side = "left",
                seed = 4
            ) |>
            assign_grnas(
                method = "mixture",
                print_progress = FALSE,
                BPPARAM = BiocParallel::SnowParam(2)
            )
    )
    expect_length(shipped, 1L)
    worker <- environment(shipped[[1]])$analyze_given_grna_ids
    hyperparameters <- environment(worker)$grna_assignment_hyperparameters
    # the formula, which holds the cell covariates in its environment, stays
    expect_false("formula_object" %in% names(hyperparameters))
    expect_false(any(vapply(
        as.list(environment(worker)),
        function(x) inherits(x, "formula"),
        logical(1)
    )))
    serial <- suppressMessages(
        sceptre_object |>
            set_analysis_parameters(
                discovery_pairs = construct_positive_control_pairs(sceptre_object),
                side = "left",
                seed = 4
            ) |>
            assign_grnas(method = "mixture", print_progress = FALSE)
    )
    expect_identical(
        parallel@initial_grna_assignment_list,
        serial@initial_grna_assignment_list
    )
})


test_that("counts_for_workers() keeps the rows that are read", {
    counts <- methods::as(
        Matrix::rsparsematrix(6, 4, 0.5),
        "RsparseMatrix"
    )
    rownames(counts) <- letters[1:6]
    kept <- counts_for_workers(counts, c("e", "b"), BiocParallel::SnowParam(2))
    expect_s4_class(kept, "dgRMatrix")
    expect_identical(as.matrix(kept), as.matrix(counts[c("e", "b"), ]))
    expect_identical(counts_for_workers(counts, c("e", "b"), NULL), counts)
    # most rows read: no subset (it would only add a copy in the master)
    expect_identical(
        counts_for_workers(counts, letters[1:4], BiocParallel::SnowParam(2)),
        counts
    )
})


test_that("counts_for_workers() makes no copy for forked workers", {
    skip_on_os("windows")
    counts <- methods::as(Matrix::rsparsematrix(6, 4, 0.5), "RsparseMatrix")
    rownames(counts) <- letters[1:6]
    # forked workers share the master's memory; a subset would be held in
    # addition to the counts
    expect_identical(
        counts_for_workers(counts, c("e", "b"), BiocParallel::MulticoreParam(2)),
        counts
    )
})
