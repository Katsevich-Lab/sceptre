# Acquire each resource with cleanup already registered. In particular, errors
# and interrupts from printing a closed/unreadable store must not leave stdout
# redirected. Existing caller output/message sinks are left in place.
write_analysis_summary <- function(object, path) {
    connection <- file(path, open = "wt")
    on.exit(close(connection), add = TRUE)
    depth <- sink.number(type = "output")
    sink(connection, type = "output")
    on.exit(
        {
            while (sink.number(type = "output") > depth) {
                sink(type = "output")
            }
        },
        add = TRUE,
        after = FALSE
    )
    print(object)
    invisible(NULL)
}
