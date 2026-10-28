test_that("complement calibration bounds the permutation sample size by the cell count", {
    # Record the sampler request instead of drawing. Two overlapping NT gRNAs
    # of four cells each cover five cells; their sizes sum to eight.
    local_mocked_bindings(fisher_yates_sampler = function(n_tot, M, B) {
        list(n_tot = n_tot, M = M, B = B)
    })
    assignments <- list(indiv_nt_grna_idxs = list(nt1 = 1:4, nt2 = 2:5))
    request <- get_synthetic_permutation_idxs(
        assignments,
        7L,
        TRUE,
        TRUE,
        2L,
        5L
    )
    expect_identical(request, list(n_tot = 5L, M = 5L, B = 7L))
    # Without overlap, the request is unchanged
    disjoint <- list(indiv_nt_grna_idxs = list(nt1 = 1:2, nt2 = 3:5))
    request <- get_synthetic_permutation_idxs(disjoint, 7L, TRUE, TRUE, 2L, 10L)
    expect_identical(request, list(n_tot = 10L, M = 5L, B = 7L))
})

test_that("public high-MOI calibration accepts overlapping NT unions", {
    inputs <- refactor_make_inputs()
    case <- refactor_case_definitions()$high_perm
    nt <- inputs$grna_target_data_frame$grna_id[
        inputs$grna_target_data_frame$grna_target == "non-targeting"
    ]
    n <- ncol(inputs$high_grna_matrix)
    common <- seq_len(n - 6L)
    inputs$high_grna_matrix[nt, ] <- 0
    inputs$high_grna_matrix[nt, common] <- 50
    target <- setdiff(rownames(inputs$high_grna_matrix), nt)[[1L]]
    inputs$high_grna_matrix[target, setdiff(seq_len(n), common)] <- 50
    object <- suppressMessages(refactor_prepare(inputs, case))
    object <- suppressMessages(sceptre::run_qc(
        object,
        response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 0L,
        n_nonzero_cntrl_thresh = 0L,
        additional_cells_to_remove = inputs$remove_cells
    ))
    assignments <- object@grna_assignments$indiv_nt_grna_idxs
    expect_true(
        sum(sort(lengths(assignments), decreasing = TRUE)[1:2]) >
            length(object@cells_in_use)
    )
    result <- suppressMessages(sceptre::run_calibration_check(
        object,
        calibration_group_size = 2L,
        n_calibration_pairs = 1L,
        print_progress = FALSE
    ))
    expect_equal(nrow(result@calibration_result), 1L)
    pairs <- result@negative_control_pairs
    guides <- strsplit(pairs$grna_group[[1L]], "&", fixed = TRUE)[[1L]]
    selected <- unique(unlist(assignments[guides], use.names = FALSE))
    source <- object@cells_in_use[selected]
    complement <- setdiff(object@cells_in_use, source)
    response <- inputs$response_matrix[pairs$response_id[[1L]], ]
    expect_equal(pairs$n_nonzero_trt, sum(response[source] > 0))
    expect_equal(pairs$n_nonzero_cntrl, sum(response[complement] > 0))
})

test_that("the Fisher-Yates sampler rejects a sample larger than its population", {
    expect_error(
        fisher_yates_sampler(n_tot = 1L, M = 2L, B = 1L),
        "M must not exceed n_tot",
        fixed = TRUE
    )
})
