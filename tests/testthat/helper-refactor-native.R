refactor_native_function <- function(name) {
    getFromNamespace(name, "sceptre")
}

refactor_native_bridge <- local({
    compiled <- NULL
    function(path) {
        if (is.null(compiled)) {
            compiled <- new.env(parent = baseenv())
            Rcpp::sourceCpp(path, env = compiled, rebuild = FALSE, showOutput = FALSE)
        }
        compiled
    }
})

refactor_native_platform <- function() {
    compiler <- paste(
        system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX"), stdout = TRUE),
        collapse = " "
    )
    compiler_command <- strsplit(compiler, " ", fixed = TRUE)[[1L]][[1L]]
    list(
        platform = R.version$platform,
        system = unname(Sys.info()[["sysname"]]),
        R = as.character(getRversion()),
        BH = as.character(utils::packageVersion("BH")),
        Rcpp = as.character(utils::packageVersion("Rcpp")),
        compiler = compiler,
        compiler_version = paste(
            system2(compiler_command, "--version", stdout = TRUE),
            collapse = " "
        )
    )
}

refactor_native_run_case <- function(case, bridge) {
    stopifnot(
        is.integer(case$trt_idxs),
        all(case$trt_idxs >= 1L & case$trt_idxs <= length(case$a)),
        length(case$draws) == case$B1 + case$B2 + case$B3,
        all(lengths(case$draws) >= case$n_trt)
    )
    args <- case[setdiff(names(case), c("draws", "expected_stage", "expected_sn"))]
    args$synthetic_idxs <- bridge$refactor_native_pointer(case$draws, length(case$a))
    do.call(refactor_native_function("run_low_level_test_full_v4"), args)
}

refactor_native_sample <- function(config, bridge) {
    pointer <- do.call(refactor_native_function(config$function_name), config$args)
    # Materialize the entire requested bank before taking a diagnostic prefix.
    bridge$refactor_native_materialize(pointer)
}
