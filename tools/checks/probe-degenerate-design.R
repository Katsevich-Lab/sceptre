#!/usr/bin/env Rscript
# Diagnostic, not an assertion of calibration: actual independent Poisson(5)
# outcomes, with the target indicator exactly included in the response design.
# Save every warning/error and actual result; do not alter default inference.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
library(sceptre)
set.seed(23703L)
n <- 160L
cell <- paste0("cell", seq_len(n))
gene <- paste0("r", seq_len(4L))
assigned <- rep(seq_len(4L), length.out = n)
guide <- matrix(0, 4L, n, dimnames = list(c("a1", "a2", "nt1", "nt2"), cell))
guide[cbind(assigned, seq_len(n))] <- 50
response <- matrix(stats::rpois(4L * n, 5), 4L, dimnames = list(gene, cell))
metadata <- data.frame(grna_id = rownames(guide), grna_target = c("A", "A", "non-targeting", "non-targeting"))
covariates <- data.frame(treatment = as.integer(assigned <= 2L), row.names = cell)
results <- list()
for (mechanism in c("permutations", "crt")) {
    warnings <- character()
    result <- tryCatch(withCallingHandlers({
        object <- suppressMessages(import_data(response, guide, metadata, moi = "low",
            extra_covariates = covariates, response_names = gene))
        object <- suppressMessages(set_analysis_parameters(object,
            discovery_pairs = data.frame(grna_target = "A", response_id = gene),
            formula_object = ~ treatment, control_group = "complement",
            side = "both", resampling_mechanism = mechanism,
            resampling_approximation = "no_approximation"))
        object <- suppressMessages(assign_grnas(object, method = "thresholding", threshold = 40,
            print_progress = FALSE))
        object <- suppressMessages(run_qc(object, response_n_umis_range = c(0, 1),
            response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 7L,
            n_nonzero_cntrl_thresh = 7L))
        completed <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE,
            output_amount = 3L))
        list(status = "completed", table = get_result(completed, "run_discovery_analysis"))
    }, warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
    }), error = function(condition) list(status = "error", message = conditionMessage(condition),
        classes = class(condition)))
    result$warnings <- warnings
    results[[mechanism]] <- result
}
saveRDS(list(seed = 23703L, distribution = "Independent Poisson(5), all genes null",
    design = "intercept and exact target membership", replication_unit = "one dataset, four genes",
    counts = response, guide_counts = guide, covariates = covariates,
    settings = list(side = "both", approximation = "no_approximation", control = "complement"),
    results = results, session = sessionInfo()), args[[1L]])
print(results)
