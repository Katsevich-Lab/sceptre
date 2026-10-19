# Usage: Rscript --vanilla reproduce.R overlay|installed LIBRARY OUTPUT_DIR [fixed]
# Overlay mode pins the historical f9f8142 R snapshot. Installed mode runs only
# the explicitly selected candidate namespace and is intended for final rechecks.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L, args[[1L]] %in% c("overlay", "installed"))
mode <- args[[1L]]
library_path <- normalizePath(args[[2L]], mustWork = TRUE)
output_dir <- args[[3L]]
expect_fixed <- length(args) >= 4L && identical(args[[4L]], "fixed")
.libPaths(c(library_path, .libPaths()))
library(sceptre)
stopifnot(identical(normalizePath(find.package("sceptre")),
                    normalizePath(file.path(library_path, "sceptre"))))
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[[1L]]
here <- dirname(normalizePath(sub("^--file=", "", script_arg)))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
namespace <- asNamespace("sceptre")
e <- namespace
if (identical(mode, "overlay")) {
    e <- new.env(parent = namespace)
    source_dir <- file.path(here, "source-f9f8142", "R")
    sys.source(file.path(source_dir, "scientific_identifiers.R"), e)
    excluded <- c("s4_class.R", "RcppExports.R", "sceptre.R", "scientific_identifiers.R")
    for (file in list.files(source_dir, full.names = TRUE)) {
        if (!basename(file) %in% excluded) sys.source(file, e)
    }
}
stopifnot(all(vapply(c("prepare_execution", "execute_work_unit", "collect_execution"),
                       function(name) is.function(get(name, e)), logical(1))))
inputs <- readRDS(file.path(here, "inputs.rds"))
case <- readRDS(file.path(here, "cases.rds"))$high_complement_crt
set.seed(20261019)
object <- suppressMessages({
    value <- e$import_data(response_matrix = inputs$response_matrix,
        grna_matrix = inputs$high_grna_matrix,
        grna_target_data_frame = inputs$grna_target_data_frame,
        extra_covariates = inputs$extra_covariates, moi = "high")
    value <- e$set_analysis_parameters(value,
        discovery_pairs = inputs$discovery_pairs,
        positive_control_pairs = inputs$positive_control_pairs,
        formula_object = stats::as.formula("~ batch + z", env = baseenv()),
        side = case$side, control_group = "complement",
        resampling_mechanism = "crt", resampling_approximation = "skew_normal",
        grna_integration_strategy = "union")
    value <- e$assign_grnas(value, method = "thresholding", threshold = 40,
        print_progress = FALSE)
    e$run_qc(value, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 7L,
        n_nonzero_cntrl_thresh = 7L, additional_cells_to_remove = inputs$remove_cells)
})
warm <- suppressMessages(e$run_discovery_analysis(object, output_amount = 1L,
    print_progress = FALSE))
manifest <- e$prepare_execution(warm, "run_discovery_analysis", n_work_units = 2L,
    output_amount = 1L)
proposals <- lapply(vapply(manifest$units, function(unit) unit$id, character(1)),
    function(id) e$execute_work_unit(warm, manifest, id))
proposals <- unserialize(serialize(proposals, NULL))
# A matching control must collect before any fault injection is interpreted.
control <- e$collect_execution(warm, manifest, rev(proposals))
stopifnot(e$result_record_is_current(control, "run_discovery_analysis"))

# Find full ordered-cell vectors structurally rather than relying on field names.
ordered_cells <- warm@preparation@analysis_cells
cell_paths <- function(value, path) {
    found <- character()
    walk <- function(x, where) {
        if (is.environment(x) || is.function(x) || typeof(x) == "externalptr") return(invisible(NULL))
        if (is.integer(x) && identical(x, ordered_cells)) found <<- c(found, where)
        if (isS4(x)) {
            for (name in methods::slotNames(x)) walk(methods::slot(x, name), paste0(where, "@", name))
        } else if (is.list(x)) {
            for (i in seq_along(x)) walk(x[[i]], paste0(where, "$",
                if (is.null(names(x))) i else names(x)[[i]]))
        }
        attrs <- attributes(x)
        skip <- c("class", if (isS4(x)) methods::slotNames(x))
        for (name in setdiff(names(attrs), skip)) walk(attrs[[name]], paste0(where, "<", name, ">"))
    }
    walk(value, path)
    found
}
manifest_cells <- cell_paths(manifest, "manifest")
proposal_cells <- lapply(seq_along(proposals), function(i)
    cell_paths(proposals[[i]], paste0("proposal[[", i, "]]")))

# Deliberate semantic fault injection, not an adversarial-forgery claim:
# recompute the ordinary payload checksum so the fit-conflict guard is tested.
bad <- unserialize(serialize(proposals, NULL))
unit_index <- which(vapply(bad, function(proposal)
    any(vapply(proposal$perturbation_fits, function(fit) !is.null(fit), logical(1))), logical(1)))[[1L]]
group <- names(Filter(Negate(is.null), bad[[unit_index]]$perturbation_fits))[[1L]]
coefficients <- bad[[unit_index]]$perturbation_fits[[group]]$fitted_coefs
position <- which(is.finite(coefficients))[[1L]]
original <- coefficients[[position]]
bad[[unit_index]]$perturbation_fits[[group]]$fitted_coefs[[position]] <- original + 1
stopifnot(identical(bad[[unit_index]]$perturbation_fits[[group]]$identity,
                    warm@fits@perturbation_fits[[group]]$identity))
bad[[unit_index]]$identity <- NULL
bad[[unit_index]]$identity <- e$execution_digest(bad[[unit_index]])
before <- serialize(warm, NULL, version = 2L)
collected <- tryCatch(e$collect_execution(warm, manifest, bad), error = identity)
rejected <- inherits(collected, "error")
unchanged <- identical(before, serialize(warm, NULL, version = 2L))
stopifnot(unchanged)
result <- list(mode = mode, package = find.package("sceptre"),
    package_version = as.character(utils::packageVersion("sceptre")),
    R = R.version.string,
    overlay_source = if (identical(mode, "overlay")) "f9f8142bd404c6da1dc39eca2c993f78cb6e5eb1" else NULL,
    loaded_runtime_certificate = e$execution_runtime_certificate(warm),
    manifest_cell_paths = manifest_cells,
    proposal_cell_paths = proposal_cells,
    matching_control_collected_current = TRUE,
    conflict_rejected = rejected,
    conflict_error = if (rejected) conditionMessage(collected) else NULL,
    original_coefficient = original,
    injected_coefficient = original + 1,
    collected_coefficient = if (rejected) NULL else collected@fits@perturbation_fits[[group]]$fitted_coefs[[position]],
    accepted_result_current = if (rejected) NULL else e$result_record_is_current(collected, "run_discovery_analysis"),
    caller_unchanged = unchanged,
    compact_proposal_certificates = !any(lengths(proposal_cells)),
    group = group, unit_id = bad[[unit_index]]$unit_id,
    session = capture.output(sessionInfo()))
result$all_fixed <- result$conflict_rejected && result$compact_proposal_certificates
saveRDS(result, file.path(output_dir, "reproduction-results.rds"))
print(result[setdiff(names(result), c("session", "loaded_runtime_certificate"))])
writeLines(capture.output(dput(result[setdiff(names(result), "session")])),
    file.path(output_dir, "reproduction-results.txt"))
if (expect_fixed) stopifnot(result$all_fixed)
