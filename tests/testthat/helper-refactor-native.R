refactor_native_function <- function(name) {
    getFromNamespace(name, "sceptre")
}

refactor_native_bridge <- local({
    compiled <- NULL
    function(path) {
        if (is.null(compiled)) {
            environment <- new.env(parent = baseenv())
            Rcpp::sourceCpp(path, env = environment, rebuild = FALSE, showOutput = FALSE)
            compiled <<- environment
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
    do.call(refactor_native_function("run_association_test"), args)
}

# Frozen configurations are historical inputs, never a candidate name registry.
# Keep this literal map separate so the retired slow CRT oracle cannot silently
# re-enter active coverage or be replaced by a reference implementation.
refactor_native_candidate_samplers <- function() {
    c(fisher_yates = "fisher_yates_sampler",
      hybrid = "hybrid_fisher_iwor_sampler",
      crt_fast_production = "sample_crt_indices")
}

refactor_native_sample <- function(config, candidate_name, bridge) {
    path <- getNamespaceInfo(asNamespace("sceptre"), "path")
    if (!file.exists(file.path(path, "Meta", "package.rds"))) {
        stop("Full-bank candidate comparisons require an installed sceptre namespace.")
    }
    pointer <- do.call(refactor_native_function(candidate_name), config$args)
    # Actual always comes from the installed candidate. The independent frozen
    # source oracle below supplies expected; materialize the entire B bank.
    bridge$refactor_native_materialize(pointer)
}

refactor_native_reference <- function() {
    list(
        commit = "21f9ea098b69c444a884a49bc254de711968e3b5",
        source = "src/generate_samples_functions.cpp",
        git_blob = "5a839ba737db9d6ef0e5e031adb7486d3692ff91",
        sha256 = "deb2dee4338f5779588408f4d80caa9a42d64e1df5e2a5c519c8361dd1c06cfd"
    )
}
