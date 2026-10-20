#!/usr/bin/env Rscript
# Diagnostic of the deliberately limited finite-observed-score guard, not a
# calibration experiment or general estimability certificate.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
library(sceptre)
eval(parse(file.path(args[[1L]], "tests/testthat/test-observed-score-validity.R"))[[1L]])
records <- list()
for (seed in c(23703L, 23705L, 23711L)) for (mechanism in c("permutations", "crt")) {
    for (control in c("complement", "nt_cells")) for (side in c("left", "right", "both")) for (gene in paste0("r", 1:4)) {
        object <- observed_score_public_object(side = side, mechanism = mechanism, control = control,
            genes = gene, seed = seed)
        warnings <- character()
        result <- tryCatch(withCallingHandlers({
            completed <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE,
                output_amount = 3L))
            table <- as.data.frame(get_result(completed, "run_discovery_analysis"))
            list(status = "completed", table = table[, intersect(c("response_id", "pass_qc",
                "p_value", "z_orig", "significant"), names(table)), drop = FALSE])
        }, warning = function(cause) {
            warnings <<- c(warnings, conditionMessage(cause))
            invokeRestart("muffleWarning")
        }), error = function(cause) list(status = "error", message = conditionMessage(cause),
            classes = class(cause)))
        records[[length(records) + 1L]] <- c(list(seed = seed, mechanism = mechanism,
            control = control, side = side, gene = gene, warnings = warnings), result)
    }
}
summary <- do.call(rbind, lapply(records, function(x) data.frame(seed = x$seed,
    mechanism = x$mechanism, control = x$control, side = x$side, gene = x$gene, status = x$status,
    z = if (is.null(x$table)) NA_real_ else x$table$z_orig,
    p = if (is.null(x$table)) NA_real_ else x$table$p_value,
    significant = if (is.null(x$table)) NA else x$table$significant)))
saveRDS(list(scope = "Exact-treatment collinearity singleton diagnostic; three datasets, not calibration",
    records = records, summary = summary, session = sessionInfo()), args[[2L]])
print(summary)
cat("Completed finite significant singletons:", sum(summary$status == "completed" & summary$significant, na.rm = TRUE), "\n")
