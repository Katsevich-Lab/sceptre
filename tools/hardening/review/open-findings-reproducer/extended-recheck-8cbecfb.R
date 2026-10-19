# Run with the same installed LIBRARY OUTPUT_DIR fixed arguments as reproduce.R.
# This wrapper lives beside reproduce.R so its fixture resolution stays explicit.
arg <- grep("^--file=", commandArgs(), value = TRUE)[[1L]]
directory <- dirname(normalizePath(sub("^--file=", "", arg)))
source(file.path(directory, "reproduce.R"))
stopifnot(identical(mode, "installed"), result$all_fixed)
checks <- list(original_conflict_rejected = rejected,
    original_scope_vectors_absent = result$compact_proposal_certificates)
roundtrip <- function(x) unserialize(serialize(x, NULL))
recertify <- function(proposal) {
    proposal$identity <- NULL
    proposal$identity <- e$execution_digest(proposal)
    proposal
}
before <- serialize(warm, NULL, version = 2L)
for (kind in c("diagnostics", "remove")) {
    changed <- roundtrip(proposals)
    if (kind == "diagnostics") {
        changed[[unit_index]]$perturbation_fits[[group]]$diagnostics$iterations <-
            changed[[unit_index]]$perturbation_fits[[group]]$diagnostics$iterations + 1L
    } else changed[[unit_index]]$perturbation_fits[group] <- list(NULL)
    changed[[unit_index]] <- recertify(changed[[unit_index]])
    outcome <- tryCatch(e$collect_execution(warm, manifest, changed), error = identity)
    checks[[paste0(kind, "_conflict_rejected")]] <- inherits(outcome, "error") &&
        grepl("Conflicting perturbation fit", conditionMessage(outcome), fixed = TRUE)
    checks[[paste0(kind, "_input_unchanged")]] <- identical(before, serialize(warm, NULL, version = 2L))
}
checks$matching_warm_fit_preserved <- identical(control@fits@perturbation_fits,
    warm@fits@perturbation_fits)
checks$full_context_identities_restored <- identical(control@fits@response_identity,
    warm@fits@response_identity) && identical(control@fits@perturbation_identity,
    warm@fits@perturbation_identity)

# A stale per-entry treatment identity must miss and be replaced, rather than
# being treated as a conflicting valid warm entry. No input count file is edited.
stale <- warm
stale@fits@perturbation_fits[[group]]$identity$treatment <- strrep("0", 64L)
stale_manifest <- e$prepare_execution(stale, "run_discovery_analysis", n_work_units = 2L)
stale_proposals <- lapply(vapply(stale_manifest$units, function(unit) unit$id, character(1)),
    function(id) e$execute_work_unit(stale, stale_manifest, id))
fresh <- e$collect_execution(stale, stale_manifest, rev(stale_proposals))
checks$stale_entry_replaced <- identical(fresh@fits@perturbation_fits[[group]],
    warm@fits@perturbation_fits[[group]])
checks$stale_replacement_current <- e$result_record_is_current(fresh, "run_discovery_analysis")

# Retained public proposal fit attributes must still be isolated after admission.
stored_names <- names(control@fits@perturbation_fits[[group]]$fitted_coefs)
coefs <- proposals[[unit_index]]$perturbation_fits[[group]]$fitted_coefs
data.table::setattr(coefs, "names", rep("caller-only-probe", length(coefs)))
checks$perturbation_proposal_attribute_isolated <- identical(stored_names,
    names(control@fits@perturbation_fits[[group]]$fitted_coefs))
response <- names(proposals[[unit_index]]$response_fits)[[1L]]
stored_names <- names(control@fits@response_fits[[response]]$fitted_coefs)
coefs <- proposals[[unit_index]]$response_fits[[response]]$fitted_coefs
data.table::setattr(coefs, "names", rep("caller-only-probe", length(coefs)))
checks$response_proposal_attribute_isolated <- identical(stored_names,
    names(control@fits@response_fits[[response]]$fitted_coefs))
print(checks)
saveRDS(list(checks = checks, all_passed = all(unlist(checks)),
    installed_package = find.package("sceptre"), runtime = result$loaded_runtime_certificate),
    file.path(output_dir, "extended-results.rds"))
writeLines(capture.output(dput(checks)), file.path(output_dir, "extended-results.txt"))
stopifnot(all(unlist(checks)))
