# Defensively unwind only resources acquired by this fixture, even when an
# assertion or the tested implementation fails.
with_output_sinks <- function(callback) {
    depth <- sink.number(type = "output")
    message_sink <- sink.number(type = "message")
    output_path <- tempfile()
    message_path <- tempfile()
    output <- file(output_path, "wt")
    messages <- file(message_path, "wt")
    on.exit({
        while (sink.number(type = "output") > depth) sink(type = "output")
        if (message_sink == 2L) sink(type = "message") else {
            sink(getConnection(message_sink), type = "message")
        }
        close(output)
        close(messages)
        unlink(c(output_path, message_path))
    })
    sink(output)
    sink(messages, type = "message")
    callback(depth + 1L, sink.number(type = "message"))
    cat("caller output restored\n")
    message("caller messages restored")
    flush(output)
    flush(messages)
    expect_identical(readLines(output_path), "caller output restored")
    expect_identical(readLines(message_path), "caller messages restored")
}

test_that("summary resources restore caller sinks on success, errors and interrupts", {
    run <- refactor_private("write_analysis_summary")
    directory <- tempfile()
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    path <- file.path(directory, "summary.txt")
    with_output_sinks(function(depth, message_sink) {
        run("normal output", path)
        expect_identical(readLines(path), '[1] "normal output"')
        expect_identical(sink.number(type = "output"), depth)
        expect_identical(sink.number(type = "message"), message_sink)
    })
    for (kind in c("error", "interrupt")) {
        original <- structure(list(message = paste("fixture", kind), call = NULL,
            fixture = 42L), class = c(paste0("fixture_", kind), kind, "condition"))
        testthat::local_mocked_bindings(print = function(...) stop(original), .package = "sceptre")
        with_output_sinks(function(depth, message_sink) {
            observed <- tryCatch(write_outputs_to_directory(new("sceptre_object"), directory),
                error = identity, interrupt = identity)
            expect_identical(observed, original)
            expect_identical(sink.number(type = "output"), depth)
            expect_identical(sink.number(type = "message"), message_sink)
        })
    }
})

test_that("successful public export restores existing caller sinks", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_nt_perm
    object <- refactor_capture(refactor_prepare(inputs, case))$value
    directory <- tempfile()
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    with_output_sinks(function(depth, message_sink) {
        expect_no_error(write_outputs_to_directory(object, directory))
        expect_identical(sink.number(type = "output"), depth)
        expect_identical(sink.number(type = "message"), message_sink)
        expect_true(file.exists(file.path(directory, "analysis_summary.txt")))
    })
})

test_that("public export of a saved closed backing cannot leak a sink", {
    skip_if_not_installed("ondisc")
    inputs <- readRDS(file.path(refactor_fixture_dir(), "backed-inputs.rds"))
    case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_nt_perm
    directory <- tempfile()
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    object <- refactor_capture(refactor_import(inputs, case, TRUE, directory))$value
    write_ondisc_backed_sceptre_object(object, directory)
    closed <- readRDS(file.path(directory, "sceptre_object.rds"))
    with_output_sinks(function(depth, message_sink) {
        expect_error(write_outputs_to_directory(closed, file.path(directory, "export")), "closed")
        expect_identical(sink.number(type = "output"), depth)
        expect_identical(sink.number(type = "message"), message_sink)
    })
})
