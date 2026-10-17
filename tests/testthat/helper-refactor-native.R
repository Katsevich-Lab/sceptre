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

# The three currently supported production samplers; no retired implementation.
refactor_native_candidate_samplers <- function() {
    c(fisher_yates = "fisher_yates_sampler",
      hybrid = "hybrid_fisher_iwor_sampler",
      crt_fast_production = "sample_crt_indices")
}

refactor_native_sample <- function(config, candidate_name, bridge) {
    pointer <- do.call(refactor_native_function(candidate_name), config$args)
    # Materialize the entire bank from the current package namespace.
    bridge$refactor_native_materialize(pointer)
}
