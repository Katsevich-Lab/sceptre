# Run from the repository root. Compare available numerical payloads, not the
# deliberately changed object schema or empty not-yet-run table placeholders.
validation_dir <- Sys.getenv("SCEPTRE_PARITY_OUTPUT",
    file.path("tools", "hardening", "validation"))
current <- readRDS(file.path(validation_dir, "parity-current.rds"))
reference <- readRDS(file.path("tools", "components", "schema3", "validation", "reference-matched-budgets.rds"))
field_map <- c(cells_in_use = "analysis_cells", grna_assignments = "analysis_groups",
    cell_removal_metrics = "cell_removal_metrics", M_matrix = "nt_nonzero_matrix",
    n_nonzero_tot_vector = "control_nonzero_totals",
    discovery_pairs_with_info = "discovery_pairs_with_info",
    positive_control_pairs_with_info = "positive_control_pairs_with_info",
    side_code = "side_code", response_precomputations = "response_precomputations")
checks <- lapply(names(current), function(name) {
    old <- reference[[name]]; new <- current[[name]]
    compare <- list(warnings = identical(old$warnings, new$warnings),
        rng_after = identical(old$value$rng_after, new$value$rng_after),
        calibrated_calibration = all.equal(old$value$calibrated$calibration_result,
            new$value$calibrated$calibration_result, tolerance = 0))
    for (table in c("calibration_result", "power_result", "discovery_result")) {
        compare[[paste0("final/", table)]] <- all.equal(old$value$final[[table]],
            new$value$final[[table]], tolerance = 0)
    }
    for (stage in c("prepared", "calibrated", "final")) {
        for (field in names(field_map)) {
            compare[[paste(stage, field_map[[field]], sep = "/")]] <- all.equal(
                old$value[[stage]][[field]], new$value[[stage]][[field_map[[field]]]],
                tolerance = 0)
        }
    }
    cat(name, ":", if (all(vapply(compare, isTRUE, logical(1)))) "EXACT MATCH" else "DIFFERENCE", "\n")
    if (!all(vapply(compare, isTRUE, logical(1)))) print(compare)
    compare
})
names(checks) <- names(current)
saveRDS(checks, file.path(validation_dir, "parity-checks.rds"))
jsonlite::write_json(list(tolerance = 0, cases = length(checks),
    comparisons_per_case = lengths(checks),
    all_equal = all(vapply(checks, function(x) all(vapply(x, isTRUE, logical(1))), logical(1))),
    checks = checks, preparation_field_map = as.list(field_map)),
    file.path(validation_dir, "parity-summary.json"), pretty = TRUE, auto_unbox = TRUE)
stopifnot(all(vapply(checks, function(x) all(vapply(x, isTRUE, logical(1))), logical(1))))
