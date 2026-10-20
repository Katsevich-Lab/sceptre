# Actual independent count outcomes, distinct from raw native NaN injection.
observed_score_public_object <- function(side = "both", mechanism = "permutations",
                                         control = "complement", near = FALSE,
                                         genes = paste0("r", 1:4), seed = 23703L) {
    input <- withr::with_seed(seed, {
        n <- 160L
        assigned <- rep(1:4, length.out = n)
        cells <- paste0("cell", seq_len(n))
        response <- matrix(stats::rpois(4L * n, 5), 4L,
            dimnames = list(paste0("r", 1:4), cells))
        guide <- matrix(0, 4L, n, dimnames = list(c("a1", "a2", "nt1", "nt2"), cells))
        guide[cbind(assigned, seq_len(n))] <- 50
        conditioning <- as.numeric(assigned <= 2L)
        if (near) conditioning <- conditioning + 0.02 * sin(seq_len(n))
        list(response = response, guide = guide,
            covariates = data.frame(conditioning = conditioning, row.names = cells))
    })
    metadata <- data.frame(grna_id = rownames(input$guide),
        grna_target = c("A", "A", "non-targeting", "non-targeting"))
    object <- suppressMessages(import_data(input$response, input$guide, metadata, moi = "low",
        extra_covariates = input$covariates, response_names = rownames(input$response)))
    object <- suppressMessages(set_analysis_parameters(object,
        discovery_pairs = data.frame(grna_target = "A", response_id = genes),
        formula_object = stats::as.formula("~ conditioning", env = baseenv()),
        control_group = control, side = side, resampling_mechanism = mechanism,
        resampling_approximation = "no_approximation"))
    object <- suppressMessages(assign_grnas(object, method = "thresholding", threshold = 40,
        print_progress = FALSE))
    suppressMessages(run_qc(object, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 7L,
        n_nonzero_cntrl_thresh = 7L))
}

observed_score_capture <- function(code) {
    warnings <- character()
    result <- tryCatch(withCallingHandlers(force(code), warning = function(cause) {
        warnings <<- c(warnings, conditionMessage(cause))
        invokeRestart("muffleWarning")
    }), error = identity)
    # Exact treatment prediction can cause a documented glm convergence warning.
    expect_true(all(grepl("^glm.fit: algorithm did not converge$|^glm.fit: fitted probabilities numerically 0 or 1 occurred$",
                          warnings)))
    result
}

observed_score_check_failure <- function(failure) {
    expect_s3_class(failure, "sceptre_worker_error")
    expect_s3_class(failure$cause, "sceptre_pair_error")
    expect_true(failure$cause$response_id %in% paste0("r", 1:4))
    expect_identical(failure$cause$grna_group, "A")
    expect_match(conditionMessage(failure), "Observed association score is not finite", fixed = TRUE)
    expect_match(conditionMessage(failure), "required covariates", fixed = TRUE)
    expect_identical(failure$cause$parent, failure$cause$cause)
}

test_that("public undefined observed scores abort every tail and control path", {
    for (mechanism in c("permutations", "crt")) for (control in c("complement", "nt_cells")) {
        for (side in c("left", "right", "both")) {
            object <- observed_score_public_object(side, mechanism, control)
            before <- serialize(object, NULL)
            failure <- observed_score_capture(suppressMessages(run_discovery_analysis(object,
                print_progress = FALSE)))
            observed_score_check_failure(failure)
            expect_identical(serialize(object, NULL), before)
            expect_false(sceptre:::result_record_is_current(object, "run_discovery_analysis"))
        }
    }
})

test_that("distributed and fork failures retain pair identity and abort the family", {
    for (mechanism in c("permutations", "crt")) for (control in c("complement", "nt_cells")) {
        object <- observed_score_public_object(mechanism = mechanism, control = control)
        execution <- prepare_execution(object, "run_discovery_analysis", 1L)
        failure <- observed_score_capture(execute_work_unit(object, execution, execution$units[[1L]]$id))
        observed_score_check_failure(failure)
    }
    skip_on_os("windows")
    object <- observed_score_public_object()
    failure <- observed_score_capture(suppressMessages(run_discovery_analysis(object,
        parallel = TRUE, n_processors = 2L, print_progress = FALSE)))
    observed_score_check_failure(failure)
})

test_that("finite near-collinear treatment tests remain available", {
    for (mechanism in c("permutations", "crt")) for (control in c("complement", "nt_cells")) {
        object <- observed_score_public_object(mechanism = mechanism, control = control, near = TRUE)
        completed <- observed_score_capture(suppressMessages(run_discovery_analysis(object,
            print_progress = FALSE, output_amount = 3L)))
        expect_s4_class(completed, "sceptre_object")
        table <- get_result(completed, "run_discovery_analysis")
        expect_true(all(is.finite(table$z_orig)))
        expect_true(all(table$p_value >= 0 & table$p_value <= 1))
        expect_true(sceptre:::result_record_is_current(completed, "run_discovery_analysis"))
    }
})

test_that("native observed NaN/Inf reject before drawing and tiny finite variance works", {
    run <- sceptre:::run_association_test
    bridge <- refactor_native_bridge(test_path("fixtures", "refactor", "native-bridge.cpp"))
    bank <- bridge$refactor_native_pointer(list(0L, 1L), 2L)
    call <- function(a, w, side, pointer = bank) run(y = c(1, 2), mu = c(1, 1),
        a = a, w = w, D = matrix(0, 0, 2), trt_idxs = 1L, n_trt = 1L,
        use_all_cells = FALSE, synthetic_idxs = pointer, B1 = 2L, B2 = 0L, B3 = 0L,
        fit_parametric_curve = FALSE, return_resampling_dist = TRUE, side_code = side)
    set.seed(23704L)
    seed <- .Random.seed
    for (side in -1:1) {
        for (a in list(c(0, 0), c(1, 0), c(-1, 0))) {
            expect_error(call(a, c(0, 1), side), "Observed association score is not finite")
            expect_error(call(a, c(0, 1), side, NULL), "Observed association score is not finite")
        }
        expect_error(call(c(NA_real_, 0), c(1, 1), side), "Observed association score is not finite")
        result <- call(c(1e-12, -1e-12), c(1e-20, 1e-20), side)
        expect_equal(result$z_orig, 0.01)
        expect_equal(result$resampling_dist, c(0.01, -0.01))
        expect_identical(result$stage, 1L)
    }
    expect_identical(.Random.seed, seed)
})
