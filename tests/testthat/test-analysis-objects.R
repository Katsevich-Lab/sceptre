# The Screen/Analysis objects and their files. These tests cover sceptre's own
# mapping; SingleCellExperiment, MuData, anndataR and arrow test their formats.

analysis_workflow <- function() {
    data(lowmoi_example_data, package = "sceptre", envir = environment())
    object <- import_data(
        response_matrix = lowmoi_example_data$response_matrix,
        grna_matrix = lowmoi_example_data$grna_matrix,
        extra_covariates = lowmoi_example_data$extra_covariates,
        grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
        moi = "low"
    )
    positive <- construct_positive_control_pairs(object)
    discovery <- construct_trans_pairs(
        object,
        positive_control_pairs = positive,
        pairs_to_exclude = "pc_pairs"
    )[1:30, ]
    steps <- list(
        parameters = function(o) {
            set_analysis_parameters(
                o,
                discovery_pairs = discovery,
                positive_control_pairs = positive
            )
        },
        assign = function(o) assign_grnas(o, method = "thresholding"),
        qc = function(o) run_qc(o),
        discovery = function(o) {
            run_discovery_analysis(o, print_progress = FALSE)
        }
    )
    list(start = object, steps = steps)
}

run_steps <- function(object, steps) {
    for (step in steps) {
        object <- step(object)
    }
    object
}

test_that("every sceptre_object slot belongs to exactly one part", {
    mapped <- unname(unlist(analysis_layout))
    expect_setequal(mapped, methods::slotNames("sceptre_object"))
    expect_false(anyDuplicated(mapped) > 0L)
})

test_that("converting to an Analysis and back loses nothing at any stage", {
    workflow <- analysis_workflow()
    object <- workflow$start
    expect_identical(as_sceptre_object(as_analysis(object)), object)
    for (step in workflow$steps) {
        object <- suppressMessages(step(object))
        expect_identical(as_sceptre_object(as_analysis(object)), object)
    }
    screen <- as_analysis(object)@screen
    expect_s4_class(screen, "SingleCellExperiment")
    expect_identical(SingleCellExperiment::altExpNames(screen), "gRNA")
})

skip_without_file_packages <- function() {
    for (package in c("anndataR", "arrow", "jsonlite", "MuData")) {
        skip_if_not_installed(package)
    }
}

# Values that files cannot hold: a formula's environment, data.table's internal
# pointer, and whether a table was a tibble.
without_r_internals <- function(x) {
    if (inherits(x, "tbl_df")) {
        x <- as.data.frame(x)
    }
    if (is.data.frame(x)) {
        attr(x, ".internal.selfref") <- NULL
    }
    if (inherits(x, "formula")) {
        environment(x) <- globalenv()
    }
    if (is.list(x) && !is.data.frame(x)) {
        x[] <- lapply(x, without_r_internals)
    }
    x
}

expect_same_object <- function(actual, expected) {
    for (name in methods::slotNames(expected)) {
        expect_identical(
            without_r_internals(methods::slot(actual, name)),
            without_r_internals(methods::slot(expected, name)),
            label = name
        )
    }
}

test_that("an analysis read back from files is the same at every stage", {
    skip_without_file_packages()
    workflow <- analysis_workflow()
    object <- workflow$start
    for (step in c(list(identity), workflow$steps)) {
        object <- suppressMessages(step(object))
        directory <- withr::local_tempfile()
        write_analysis(as_analysis(object), directory)
        expect_same_object(as_sceptre_object(read_analysis(directory)), object)
    }
    # a changed file is detected before anything is read
    arrow::write_parquet(
        data.frame(changed = 1L),
        file.path(directory, "results", "discovery.parquet")
    )
    expect_error(read_analysis(directory), "results/discovery.parquet")
    expect_error(
        write_analysis(as_analysis(object), directory),
        "must be empty"
    )
})

test_that("an analysis saved before discovery finishes with the same results", {
    skip_without_file_packages()
    workflow <- analysis_workflow()
    steps <- workflow$steps
    direct <- suppressMessages(run_steps(workflow$start, steps))
    before <- suppressMessages(run_steps(workflow$start, steps[1:3]))
    directory <- withr::local_tempdir()
    write_analysis(as_analysis(before), directory)
    resumed <- suppressMessages(
        steps$discovery(as_sceptre_object(read_analysis(directory)))
    )
    expect_identical(resumed@discovery_result, direct@discovery_result)
})

test_that("formulas inside settings survive the state file", {
    skip_without_file_packages()
    object <- analysis_workflow()$start
    object@grna_assignment_hyperparameters <- list(
        n_em_rep = 5L,
        formula_object = ~ log(grna_n_umis) + batch
    )
    directory <- withr::local_tempdir()
    write_analysis(as_analysis(object), directory)
    expect_same_object(as_sceptre_object(read_analysis(directory)), object)
})

test_that("missing gRNA locations survive the screen file", {
    skip_without_file_packages()
    data(highmoi_example_data, package = "sceptre", envir = environment())
    data(
        grna_target_data_frame_highmoi,
        package = "sceptre",
        envir = environment()
    )
    object <- import_data(
        response_matrix = highmoi_example_data$response_matrix,
        grna_matrix = highmoi_example_data$grna_matrix,
        extra_covariates = highmoi_example_data$extra_covariates,
        grna_target_data_frame = grna_target_data_frame_highmoi,
        moi = "high"
    )
    expect_true(anyNA(object@grna_target_data_frame$chr))
    directory <- withr::local_tempdir()
    write_analysis(as_analysis(object), directory)
    read_back <- as_sceptre_object(read_analysis(directory))
    expect_identical(
        read_back@grna_target_data_frame,
        object@grna_target_data_frame
    )
    expect_identical(
        read_back@covariate_data_frame,
        object@covariate_data_frame
    )
})
