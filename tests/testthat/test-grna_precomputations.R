####################################################################
##            Caching of the gRNA precomputation                  ##
####################################################################

# a small, well-conditioned design on which the logistic regression converges
make_mock_grna_design <- function(n_cells = 200L, seed = 5L) {
    set.seed(seed)
    covariate_matrix <- cbind(
        "(Intercept)" = rep(1, n_cells),
        n_umis = rnorm(n_cells),
        p_mito = rnorm(n_cells)
    )
    list(
        covariate_matrix = covariate_matrix,
        trt_idxs = sort(sample.int(n_cells, size = 60L))
    )
}


test_that("the fitted values are recovered exactly from the coefficients", {
    design <- make_mock_grna_design()

    from_fit <- perform_grna_precomputation(
        trt_idxs = design$trt_idxs,
        covariate_matrix = design$covariate_matrix,
        return_fitted_values = TRUE
    )
    from_coefs <- reconstruct_grna_fitted_probabilities(
        fitted_coefs = perform_grna_precomputation(
            trt_idxs = design$trt_idxs,
            covariate_matrix = design$covariate_matrix,
            return_fitted_values = FALSE
        ),
        covariate_matrix = design$covariate_matrix
    )

    # `identical()`, not `equal()`: a one-ulp difference would flip a Bernoulli
    # draw in `crt_index_sampler_fast()` and change a p-value
    expect_identical(from_coefs, from_fit)
})


test_that("a cache hit and a cache miss produce the same probabilities", {
    design <- make_mock_grna_design()

    miss <- get_grna_fitted_probabilities(
        grna_precomputations = list(),
        curr_grna_group = "t1",
        trt_idxs = design$trt_idxs,
        covariate_matrix = design$covariate_matrix,
        fit_idxs = design$trt_idxs
    )
    expect_false(is.null(miss$to_cache))

    hit <- get_grna_fitted_probabilities(
        grna_precomputations = list(t1 = miss$to_cache),
        curr_grna_group = "t1",
        trt_idxs = design$trt_idxs,
        covariate_matrix = design$covariate_matrix,
        fit_idxs = design$trt_idxs
    )
    # nothing to add: the cache already holds this fit
    expect_null(hit$to_cache)
    expect_identical(hit$fitted_probabilities, miss$fitted_probabilities)

    # the cached fit is a handful of coefficients rather than one value per cell
    expect_length(
        miss$to_cache$fitted_coefs,
        ncol(design$covariate_matrix)
    )
})


test_that("a stale cache entry is refit rather than replayed", {
    design <- make_mock_grna_design()
    entry <- get_grna_fitted_probabilities(
        grna_precomputations = list(),
        curr_grna_group = "t1",
        trt_idxs = design$trt_idxs,
        covariate_matrix = design$covariate_matrix,
        fit_idxs = design$trt_idxs
    )$to_cache

    # a different set of treated cells, as after a rerun of `assign_grnas()`
    expect_false(grna_precomputation_is_current(
        grna_precomp = entry,
        covariate_matrix = design$covariate_matrix,
        fit_idxs = design$trt_idxs[-1]
    ))
    # a different set of covariates, as after a change to `formula_object`
    expect_false(grna_precomputation_is_current(
        grna_precomp = entry,
        covariate_matrix = design$covariate_matrix[, 1:2, drop = FALSE],
        fit_idxs = design$trt_idxs
    ))
    # an absent entry is a cache miss, not an error
    expect_false(grna_precomputation_is_current(
        grna_precomp = NULL,
        covariate_matrix = design$covariate_matrix,
        fit_idxs = design$trt_idxs
    ))
    expect_true(grna_precomputation_is_current(
        grna_precomp = entry,
        covariate_matrix = design$covariate_matrix,
        fit_idxs = design$trt_idxs
    ))

    # a stale entry yields the fit for the design actually in hand
    refit_idxs <- design$trt_idxs[-1]
    refit <- get_grna_fitted_probabilities(
        grna_precomputations = list(t1 = entry),
        curr_grna_group = "t1",
        trt_idxs = refit_idxs,
        covariate_matrix = design$covariate_matrix,
        fit_idxs = refit_idxs
    )
    expect_identical(
        refit$fitted_probabilities,
        perform_grna_precomputation(
            trt_idxs = refit_idxs,
            covariate_matrix = design$covariate_matrix,
            return_fitted_values = TRUE
        )
    )
    expect_identical(refit$to_cache$fit_idxs, as.integer(refit_idxs))
})


test_that("a rank-deficient fit is not cached", {
    design <- make_mock_grna_design()
    # an aliased column makes the fit rank deficient; the fitted values are
    # still valid, but the coefficients contain an NA and cannot be replayed
    covariate_matrix <- cbind(
        design$covariate_matrix,
        n_umis_copy = design$covariate_matrix[, "n_umis"]
    )

    out <- suppressWarnings(get_grna_fitted_probabilities(
        grna_precomputations = list(),
        curr_grna_group = "t1",
        trt_idxs = design$trt_idxs,
        covariate_matrix = covariate_matrix,
        fit_idxs = design$trt_idxs
    ))

    expect_null(out$to_cache)
    expect_identical(
        out$fitted_probabilities,
        suppressWarnings(perform_grna_precomputation(
            trt_idxs = design$trt_idxs,
            covariate_matrix = covariate_matrix,
            return_fitted_values = TRUE
        ))
    )
})


test_that("`get_grna_precomputations()` tolerates an object without the slot", {
    # a `sceptre_object` deserialized from a version predating the slot
    old_object <- structure(list(), class = "sceptre_object")
    expect_identical(get_grna_precomputations(old_object), list())
})


####################################################################
##                     End-to-end behavior                        ##
####################################################################

# a minimal high-MOI CRT setup; `control_group` selects which of the two CRT
# workhorse functions, and hence which row space, the analysis uses
make_mock_crt_object <- function(control_group, seed = 3L) {
    set.seed(seed)
    num_cells <- 150L
    num_responses <- 8L

    grna_target_data_frame <- make_mock_grna_target_data(c(2, 2, 2), 1, 1, 10)
    grna_matrix <- rpois(num_cells * nrow(grna_target_data_frame), 5) |>
        matrix(nrow = nrow(grna_target_data_frame), ncol = num_cells) |>
        `rownames<-`(grna_target_data_frame$grna_id)
    response_matrix <- rpois(num_cells * num_responses, 5) |>
        matrix(nrow = num_responses, ncol = num_cells) |>
        `rownames<-`(paste0("response_", seq_len(num_responses)))
    discovery_pairs <- expand.grid(
        grna_target = c("t1_c1_d1", "t1_c2_d1"),
        response_id = rownames(response_matrix),
        stringsAsFactors = FALSE
    )

    import_data(
        grna_matrix = grna_matrix,
        response_matrix = response_matrix,
        grna_target_data_frame = grna_target_data_frame,
        moi = "high"
    ) |>
        set_analysis_parameters(
            discovery_pairs = discovery_pairs,
            control_group = control_group,
            resampling_mechanism = "crt"
        ) |>
        assign_grnas(method = "thresholding", threshold = 5) |>
        run_qc(
            response_n_umis_range = c(0, 1),
            response_n_nonzero_range = c(0, 1),
            n_nonzero_trt_thresh = 0,
            n_nonzero_cntrl_thresh = 0
        )
}


for (control_group in c("complement", "nt_cells")) {
    test_that(
        paste0(
            "the gRNA cache is populated and reused (control_group = ",
            control_group,
            ")"
        ),
        {
            scep <- make_mock_crt_object(control_group)
            expect_length(scep@grna_precomputations, 0L)

            # the mock design is small enough that the logistic fit separates,
            # which is incidental here but see the warning test below
            first <- suppressWarnings(
                run_discovery_analysis(scep, print_progress = FALSE)
            )
            # one entry per gRNA target analyzed
            expect_setequal(
                names(first@grna_precomputations),
                c("t1_c1_d1", "t1_c2_d1")
            )

            # a second run serves every fit from the cache and must agree with
            # the first to the last bit
            second <- run_discovery_analysis(first, print_progress = FALSE)
            expect_identical(
                second@discovery_result,
                first@discovery_result
            )
            expect_identical(
                second@grna_precomputations,
                first@grna_precomputations
            )
        }
    )
}


test_that("a warning from the fit is raised only when the fit is performed", {
    scep <- make_mock_crt_object("complement")

    # this mock design separates, so the first run warns
    expect_warning(
        first <- run_discovery_analysis(scep, print_progress = FALSE),
        "fitted probabilities numerically 0 or 1"
    )
    # the cached run performs no logistic fit, so it has nothing to warn about.
    # Any warning the gRNA regression raises is therefore emitted on the run
    # that fits, not on the runs that reuse the fit.
    expect_length(
        testthat::capture_warnings(
            run_discovery_analysis(first, print_progress = FALSE)
        ),
        0L
    )
})


test_that("a cache built for one control group is not replayed in the other", {
    complement <- suppressWarnings(
        make_mock_crt_object("complement") |>
            run_discovery_analysis(print_progress = FALSE)
    )
    nt_cells <- make_mock_crt_object("nt_cells")

    # the two analyses regress on different row spaces, so transplanting the
    # cache must change nothing: every entry is stale and is refit
    poisoned <- nt_cells
    poisoned@grna_precomputations <- complement@grna_precomputations
    expect_identical(
        suppressWarnings(
            run_discovery_analysis(poisoned, print_progress = FALSE)
        )@discovery_result,
        suppressWarnings(
            run_discovery_analysis(nt_cells, print_progress = FALSE)
        )@discovery_result
    )
})


test_that("the calibration check leaves the gRNA cache untouched", {
    scep <- suppressWarnings(
        make_mock_crt_object("complement") |>
            run_discovery_analysis(print_progress = FALSE)
    )
    cache <- scep@grna_precomputations
    expect_gt(length(cache), 0L)

    # the calibration check's gRNA groups are randomly chosen sets of
    # non-targeting gRNAs, resampled on every call: nothing durable to cache
    calib <- run_calibration_check(
        scep,
        calibration_group_size = 1,
        print_progress = FALSE
    )
    expect_identical(calib@grna_precomputations, cache)
})


test_that("the gRNA cache is reset when its inputs change", {
    scep <- make_mock_crt_object("complement")
    scep@grna_precomputations <- list(t1_c1_d1 = "stale")

    # rerunning the gRNA assignment changes which cells are treated
    expect_length(
        assign_grnas(
            scep,
            method = "thresholding",
            threshold = 6
        )@grna_precomputations,
        0L
    )
    # rerunning QC can change the cells in use
    expect_length(
        run_qc(
            scep,
            response_n_umis_range = c(0.05, 0.95)
        )@grna_precomputations,
        0L
    )
})
