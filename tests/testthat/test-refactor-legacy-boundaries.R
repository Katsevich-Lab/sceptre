legacy_plot_object <- function() {
    pairs <- data.frame(response_id = paste0("r", 1:6),
        grna_group = "g", pass_qc = rep(c(TRUE, FALSE), 3L),
        n_nonzero_trt = 1:6, n_nonzero_cntrl = 6:1)
    object <- methods::new("sceptre_object",
        data = methods::new("ScreenData", identity = "plot-fixture",
            low_moi = FALSE,
            response_matrix = set_matrix_accessibility(matrix(1, 6L, 4L,
                dimnames = list(paste0("r", 1:6), NULL))),
            grna_matrix = set_matrix_accessibility(matrix(c(0, 1, 11, 15),
                1L, dimnames = list("g", NULL))),
            grna_target_data_frame = data.frame(grna_id = "g", grna_target = "g")),
        spec = methods::new("AnalysisSpec", n_nonzero_trt_thresh = 7L,
            n_nonzero_cntrl_thresh = 7L, grna_integration_strategy = "union",
            multiple_testing_alpha = 0.1, multiple_testing_method = "BH"),
        preparation = methods::new("PreparationRecord",
            guide_groups = "g", design_matrix = matrix(1, 4L, 1L),
            guide_assignments = list(g = 3:4),
            source_groups = list(grna_group_idxs = list(g = 3:4), indiv_nt_grna_idxs = list()),
            analysis_groups = list(grna_group_idxs = list(g = 3:4), indiv_nt_grna_idxs = list()),
            analysis_cells = 1:4, cell_removal_metrics = rep(0L, 6L),
            discovery_pairs = pairs[c("response_id", "grna_group")],
            discovery_pairs_with_info = pairs))
    for (stage in c("design", "assignment", "cell_qc", "pair_qc")) {
        object <- refactor_private("commit_preparation_stage")(object, stage)
    }
    object
}

test_that("high explicit guide thresholds do not request negative expansion bins", {
    object <- legacy_plot_object()
    before <- serialize(object, NULL)
    set.seed(184L)
    seed <- .Random.seed
    for (threshold in c(16, 40)) {
        panels <- sceptre::plot_grna_count_distributions(object, grnas_to_plot = "g",
            threshold = threshold, return_indiv_plots = TRUE)
        expect_length(panels, 1L)
        expect_identical(as.character(panels[[1L]]$data$grna_expressions_bin),
                         c("0", "1", "11", "15"))
        expect_identical(panels[[1L]]$data$bin_counts, rep(1L, 4L))
    }
    # Preserve the existing exponential-bin labels in the default branch.
    default <- sceptre::plot_grna_count_distributions(object, grnas_to_plot = "g",
                                                      return_indiv_plots = TRUE)
    expect_identical(as.character(default[[1L]]$data$grna_expressions_bin),
                     c("0", "1", "10", "11"))
    expect_identical(default[[1L]]$data$bin_counts, rep(1L, 4L))
    expect_identical(.Random.seed, seed)
    expect_identical(serialize(object, NULL), before)
})

test_that("volcano colors use existing flags including significant threshold ties", {
    result <- data.frame(p_value = c(.01, .02, .02, .03),
        log_2_fold_change = c(-2, 0, .5, 2),
        significant = c(TRUE, TRUE, FALSE, FALSE))
    before <- serialize(result, NULL)
    plot <- refactor_private("make_volcano_plot")(result, p_thresh = .02)
    expect_identical(plot$data$reject, result$significant)
    expect_identical(plot$data$p_value, result$p_value)
    expect_identical(plot$data$log_2_fold_change, c(-1.5, 0, .5, 1.5))
    expect_identical(plot$layers[[2L]]$data$yintercept, .02)
    expect_identical(serialize(result, NULL), before)
})

test_that("QC forwards nondefaults while retaining default and single-plot returns", {
    object <- legacy_plot_object()
    before <- serialize(object, NULL)
    set.seed(829L)
    seed <- .Random.seed
    panels <- sceptre::plot_run_qc(object, downsample_pairs = 2L, point_size = 3,
                                  transparency = .2, return_indiv_plots = TRUE)
    expect_length(panels, 2L)
    expected_rows <- withr::with_seed(4, sample.int(6L, 2L))
    expect_identical(panels[[2L]]$data$response_id,
                     object@preparation@discovery_pairs_with_info$response_id[expected_rows])
    expect_identical(panels[[2L]]$layers[[1L]]$aes_params$size, 3)
    expect_identical(panels[[2L]]$layers[[1L]]$aes_params$alpha, .2)
    default <- sceptre::plot_run_qc(object, return_indiv_plots = TRUE)[[2L]]
    direct <- refactor_private("plot_pairwise_qc")(object)
    expect_identical(default$data, direct$data)
    expect_identical(default$layers[[1L]]$aes_params, direct$layers[[1L]]$aes_params)
    expect_identical(nrow(default$data), 6L)
    object@preparation@discovery_pairs_with_info <- object@preparation@discovery_pairs_with_info[FALSE, ]
    object@preparation@discovery_pairs <- object@preparation@discovery_pairs[FALSE, ]
    for (stage in c("design", "assignment", "cell_qc", "pair_qc")) {
        object <- refactor_private("commit_preparation_stage")(object, stage)
    }
    single <- sceptre::plot_run_qc(object, return_indiv_plots = TRUE)
    expect_s3_class(single, "ggplot")
    expect_identical(.Random.seed, seed)
    object <- unserialize(before)
    expect_identical(serialize(object, NULL), before)
})

test_that("QQ bands keep one group separately in every panel", {
    values <- c(.05, .2, .4, .6, .9)
    data <- data.frame(panel = rep(c("A", "B"), each = 10L),
        group = rep(c("g1", "g2", "g3", "g4"), each = 5L),
        p = rep(values, 4L))
    before <- serialize(data, NULL)
    set.seed(114L)
    seed <- .Random.seed
    for (panels in list(data, data[data$panel == "A", ])) {
        plot <- ggplot2::ggplot(panels, ggplot2::aes(y = p, group = group)) +
            refactor_private("stat_qq_band")() + ggplot2::facet_wrap(~panel)
        built <- ggplot2::ggplot_build(plot)
        band <- built$data[[1L]]
        expect_setequal(as.character(unique(band$PANEL)),
                        as.character(built$layout$layout$PANEL))
        # Independent uniform order-statistic bands; no new ranking, tie or
        # sampling method is selected by the per-panel group repair.
        for (panel in unique(band$PANEL)) {
            actual <- band[band$PANEL == panel, ]
            actual <- actual[order(actual$x), ]
            expect_equal(actual$x, stats::ppoints(5L), tolerance = 1e-12)
            expect_equal(actual$ymin, stats::qbeta(.025, 1:5, 5:1), tolerance = 1e-12)
            expect_equal(actual$ymax, stats::qbeta(.975, 1:5, 5:1), tolerance = 1e-12)
        }
    }
    malformed <- data.frame(y = 1:5, PANEL = 1L, group = c(1L, 1L, 2L, 2L, 2L))
    expect_error(refactor_private("StatQQBand")$setup_data(malformed, list()),
                 "Within each panel, you must have the same number of points per QQ plot!", fixed = TRUE)
    expect_identical(.Random.seed, seed)
    expect_identical(serialize(data, NULL), before)
})

legacy_mixture_parameters <- function() list(n_em_rep = 1,
    pi_guess_range = c(.01, .1), g_pert_guess_range = c(0, 20),
    n_nonzero_cells_cutoff = 1, backup_threshold = 1,
    probability_threshold = .8, formula_object = ~ 1)

test_that("mixture perturbation bounds use the documented perturbation range", {
    object <- legacy_plot_object()
    check <- refactor_private("check_assign_grna_inputs")
    for (upper in c(19.9, 20)) {
        parameters <- legacy_mixture_parameters()
        parameters$g_pert_guess_range[2L] <- upper
        expect_null(check(object, "mixture", parameters, "auto"))
    }
    parameters$g_pert_guess_range[2L] <- 20.1
    expect_error(check(object, "mixture", parameters, "auto"),
                 "`g_pert_guess_range` should be a numeric vector of length 2", fixed = TRUE)
    parameters$pi_guess_range[2L] <- 2
    expect_error(check(object, "mixture", parameters, 1L),
                 "`pi_guess_range` should be a numeric vector of length 2", fixed = TRUE)
})

test_that("shared processor validation keeps existing predicates and error order", {
    check <- refactor_private("check_processor_count")
    for (value in list("auto", 2L, 2.5, Inf)) expect_null(check(value))
    for (value in list("2", 1L)) expect_error(check(value),
        "`n_processors` should be set to the string 'auto' or an integer greater than or equal to 2.", fixed = TRUE)
    expect_error(check(NA_real_), "missing value where TRUE/FALSE needed", fixed = TRUE)
    object <- legacy_plot_object()
    expect_error(refactor_private("check_assign_grna_inputs")(object, "invalid",
        list(), 1L), "`assignment_method` must be", fixed = TRUE)
    object@spec@control_group_complement <- FALSE
    object@preparation@analysis_groups <- list(indiv_nt_grna_idxs = list(nt1 = 1L, nt2 = 2L))
    expect_error(refactor_private("check_calibration_check_inputs")(object,
        n_calibration_pairs = 0L, n_processors = 1L, calibration_group_size = 1L),
        "Cannot run a calibration check on zero negative control pairs.", fixed = TRUE)
    messages <- character()
    outcome <- withCallingHandlers(tryCatch(
        refactor_private("check_discovery_analysis_inputs")(
            response_grna_group_pairs = data.frame(response_id = "r", grna_group = "A"),
            control_group_complement = TRUE, grna_target_data_frame = data.frame(),
            pc_analysis = FALSE, calibration_result = data.frame(), n_ok_pairs = 0L,
            n_processors = 1L), error = conditionMessage),
        message = function(cnd) {
            messages <<- c(messages, conditionMessage(cnd)); invokeRestart("muffleMessage")
        })
    expect_match(messages, "The calibration check (`run_calibration_check()`) should be run before the discovery analysis.", fixed = TRUE)
    expect_identical(outcome, "Zero pairs pass pairwise QC. Cannot run analysis.")
})


test_that("maximum assignment needs only aligned import summaries and retains order", {
    ids <- c("b", "a", "b", "other", "a", "b")
    fractions <- c(.8, .81, .4, 1, .9, .95)
    totals <- c(5, 4, 100, 8, 10, 5)
    requested <- c("b", "a", "b")
    before <- serialize(list(ids, fractions, totals, requested), NULL)
    set.seed(1831L)
    seed <- .Random.seed
    actual <- refactor_private("assign_grnas_to_cells_maximum")(
        max_grna = ids, max_grna_frac_umis = fractions,
        umi_fraction_threshold = .8, grna_n_umis = totals,
        min_grna_n_umis_threshold = 5, grna_ids = requested)
    expect_identical(actual$initial_assignment_list,
        stats::setNames(list(c(1L, 3L, 6L), c(2L, 5L), c(1L, 3L, 6L)), requested))
    # Existing union order is ambiguous cells first, then new low-total cells.
    expect_identical(actual$cells_w_zero_or_twoplus_grnas, c(1L, 3L, 2L))
    expect_identical(serialize(list(ids, fractions, totals, requested), NULL), before)
    expect_identical(.Random.seed, seed)
})
