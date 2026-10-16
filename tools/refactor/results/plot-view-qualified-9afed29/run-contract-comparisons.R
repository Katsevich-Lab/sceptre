out <- '/private/tmp/sceptre-plot-view-9afed29'
.libPaths(c('/private/tmp/sceptre-cl2-dense-final-909aa3a/candidate/library',
            '/private/tmp/sceptre-rf0-ondisc-library', .libPaths()))
namespace <- loadNamespace('sceptre')
checks <- logical()
check <- function(name, actual, expected) {
    passed <- identical(actual, expected)
    checks[[name]] <<- passed
    if (!passed) stop('Contract mismatch: ', name)
    invisible(actual)
}
source_functions <- function(path) {
    env <- new.env(parent = namespace)
    eval(parse(file = path, keep.source = FALSE), envir = env)
    env
}
old_score <- parse(file.path(out, 'reference/score_preparation.R'), keep.source = FALSE)
new_score <- parse(file.path(out, 'source/R/score_preparation.R'), keep.source = FALSE)
check('score module parsed expressions unchanged', new_score, old_score)
old_context <- source_functions(file.path(out, 'reference/analysis_context.R'))
new_context <- source_functions(file.path(out, 'source/R/analysis_context.R'))
expected_body <- body(old_context$make_fit_view)
last <- length(expected_body)
stopifnot(identical(expected_body[[last]][[1L]], as.name('list')),
          'control_fit_pos' %in% names(as.list(expected_body[[last]])))
expected_body[[last]][['control_fit_pos']] <- NULL
check('fit-view body changes only by deleting unused field',
      body(new_context$make_fit_view), expected_body)
check('fit-view formals unchanged', formals(new_context$make_fit_view),
      formals(old_context$make_fit_view))
for (name in setdiff(ls(old_context), 'make_fit_view')) {
    check(paste('context body unchanged', name), body(new_context[[name]]), body(old_context[[name]]))
    check(paste('context formals unchanged', name), formals(new_context[[name]]), formals(old_context[[name]]))
}
base_plot <- source_functions(file.path(out, 'reference/base-plotting_functions.R'))
new_plot <- source_functions(file.path(out, 'source/R/plotting_functions.R'))
for (name in setdiff(ls(base_plot), 'plot_run_discovery_analysis')) {
    check(paste('other plot body unchanged', name), body(new_plot[[name]]), body(base_plot[[name]]))
    check(paste('other plot formals unchanged', name), formals(new_plot[[name]]), formals(base_plot[[name]]))
}
check('discovery plot public formals unchanged', formals(new_plot$plot_run_discovery_analysis),
      formals(base_plot$plot_run_discovery_analysis))

set.seed(935L)
seed <- .Random.seed
index <- list(analysis_to_source = c(2L, 4L, 5L, 8L, 11L),
    groups_analysis = list(A = c(2L, 3L)),
    nt_to_analysis = c(4L, 1L, 5L),
    nt_groups_analysis = list(ntA = 4L, ntB = c(1L, 5L)),
    low_moi = TRUE, calibration_check = FALSE, control_group_complement = FALSE)
index_before <- serialize(index, NULL)
view <- function(i, group, expected, expected_control_source) {
    actual <- new_context$make_fit_view(i, group)
    check(paste('literal fit mapping', length(views) + 1L), actual, expected)
    controls <- setdiff(seq_along(actual$source_cells), actual$trt_fit_pos)
    check(paste('literal control sources', length(views) + 1L),
          actual$source_cells[controls], expected_control_source)
    old <- old_context$make_fit_view(i, group)
    old$control_fit_pos <- NULL
    check(paste('retained view fields unchanged', length(views) + 1L), actual, old)
    views[[length(views) + 1L]] <<- actual
}
views <- list()
view(index, 'A', list(fit_to_analysis = c(2L, 3L, 4L, 1L, 5L),
    source_cells = c(4L, 5L, 8L, 2L, 11L), trt_fit_pos = 1:2, n_trt = 2L), c(8L, 2L, 11L))
view(index, NULL, list(fit_to_analysis = 1:5,
    source_cells = c(2L, 4L, 5L, 8L, 11L), trt_fit_pos = integer(), n_trt = 0L), c(2L, 4L, 5L, 8L, 11L))
calibration <- index; calibration$calibration_check <- TRUE
view(calibration, 'ntB', list(fit_to_analysis = c(4L, 1L, 5L),
    source_cells = c(8L, 2L, 11L), trt_fit_pos = c(2L, 3L), n_trt = 2L), 8L)
view(calibration, NULL, list(fit_to_analysis = c(4L, 1L, 5L),
    source_cells = c(8L, 2L, 11L), trt_fit_pos = integer(), n_trt = 0L), c(8L, 2L, 11L))
complement <- index; complement$control_group_complement <- TRUE
view(complement, 'A', list(fit_to_analysis = 1:5,
    source_cells = c(2L, 4L, 5L, 8L, 11L), trt_fit_pos = c(2L, 3L), n_trt = 2L), c(2L, 8L, 11L))
calibration$control_group_complement <- TRUE
view(calibration, 'ntB', list(fit_to_analysis = 1:5,
    source_cells = c(2L, 4L, 5L, 8L, 11L), trt_fit_pos = c(1L, 5L), n_trt = 2L), c(4L, 5L, 8L))
check('mapping input unchanged', serialize(index, NULL), index_before)
check('mapping RNG unchanged', .Random.seed, seed)

# Full amount1 parity with the immediate pre-fix integration and qualified
# original21f parity, separating the already accepted B4 volcano-color repair.
sys.source(file.path(out, 'reference/base-qq_plot_helpers.R'), envir = base_plot)
sys.source(file.path(out, 'source/R/qq_plot_helpers.R'), envir = new_plot)
original_plot <- source_functions(file.path(out, 'reference/plotting_functions.R'))
sys.source(file.path(out, 'reference/qq_plot_helpers.R'), envir = original_plot)
core <- data.frame(response_id = c('response_a', 'response_b'),
    grna_target = c('target_a', 'target_b'),
    n_nonzero_trt = c(15L, 18L), n_nonzero_cntrl = c(20L, 21L),
    pass_qc = TRUE, p_value = c(0.0001, 0.5), fold_change = c(0.5, 1.1),
    se_fold_change = c(0.1, 0.2), log_2_fold_change = log2(c(0.5, 1.1)),
    significant = c(TRUE, FALSE))
comparisons <- list()
for (table_class in c('data.frame', 'data.table')) {
    object <- methods::new('sceptre_object')
    object@functs_called <- c(run_discovery_analysis = TRUE)
    object@multiple_testing_alpha <- 0.1
    object@discovery_result <- if (table_class == 'data.table') data.table::as.data.table(core) else core
    object@calibration_result <- data.frame(p_value = c(0.1, 0.8))
    before <- serialize(object, NULL)
    panels <- lapply(list(original21f = original_plot, base3d9fafc = base_plot, candidate9afed29 = new_plot),
        function(env) env$plot_run_discovery_analysis(object, return_indiv_plots = TRUE))
    built <- lapply(panels, function(p) lapply(p, function(g) ggplot2::ggplot_build(g)$data))
    check(paste(table_class, 'complete built amount1 data unchanged from immediate base'),
          built$candidate9afed29, built$base3d9fafc)
    for (panel in c(1L, 2L, 4L)) {
        check(paste(table_class, 'original21f QQ/text data unchanged', panel),
              built$candidate9afed29[[panel]], built$original21f[[panel]])
    }
    primary <- c('response_id', 'p_value', 'log_2_fold_change', 'significant', 'bin')
    check(paste(table_class, 'original21f volcano primary rows/downsampling unchanged'),
          as.data.frame(panels$candidate9afed29[[3L]]$data)[, primary],
          as.data.frame(panels$original21f[[3L]]$data)[, primary])
    check(paste(table_class, 'original21f displayed volcano coordinates unchanged'),
          built$candidate9afed29[[3L]][[1L]][, c('x', 'y')],
          built$original21f[[3L]][[1L]][, c('x', 'y')])
    check(paste(table_class, 'known original21f strict-threshold color defect characterized'),
          panels$original21f[[3L]]$data$reject, c(FALSE, FALSE))
    check(paste(table_class, 'current significance flag repair retained'),
          panels$candidate9afed29[[3L]]$data$reject, c(TRUE, FALSE))
    check(paste(table_class, 'amount1 input unchanged'), serialize(object, NULL), before)
    check(paste(table_class, 'amount1 RNG unchanged'), .Random.seed, seed)
    comparisons[[table_class]] <- list(panels_data = lapply(panels, function(p) lapply(p, function(g) g$data)), built_data = built)
}
utils::write.csv(data.frame(check = names(checks), passed = unname(checks)),
                 file.path(out, 'contract-comparisons.csv'), row.names = FALSE)
saveRDS(list(checks = checks, views = views, amount1_comparisons = comparisons),
        file.path(out, 'contract-comparisons.rds'))
cat(length(checks), 'contract checks passed;6 literal fit views; exact base amount1 parity; qualified original21f parity.\n')
print(sessionInfo())
