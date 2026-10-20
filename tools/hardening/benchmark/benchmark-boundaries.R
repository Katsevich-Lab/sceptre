# Run from the repository root:
# Rscript tools/hardening/benchmark/benchmark-boundaries.R
# Optional SCEPTRE_BENCH_LIBRARY selects a freshly installed package library.
# These bounded timings characterize this local workload, not large screens.
main <- function() {
root <- normalizePath('.', mustWork = TRUE)
output <- file.path(root, 'tools', 'hardening', 'benchmark', 'results')
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source_paths <- sort(c(list.files(file.path(root, 'R'), pattern = "\\.R$", full.names = TRUE),
    file.path(root, 'tools', 'hardening', 'benchmark', 'benchmark-boundaries.R')))
source_hashes <- tools::md5sum(source_paths)
names(source_hashes) <- substring(names(source_hashes), nchar(root) + 2L)

selected_library <- Sys.getenv('SCEPTRE_BENCH_LIBRARY')
if (nzchar(selected_library)) {
    .libPaths(c(normalizePath(selected_library), .libPaths()))
    library(sceptre)
} else pkgload::load_all(root, recompile = FALSE)
private <- function(name) getFromNamespace(name, 'sceptre')
source(file.path(root, 'tests', 'testthat', 'helper-refactor.R'))
options(sceptre.refactor.fixture_dir = file.path(root, 'tests', 'testthat', 'fixtures', 'refactor'))

fixture <- readRDS(file.path(root, 'tests', 'testthat', 'fixtures', 'saved-objects',
    'sceptre-schema3-list-plans-2a0938e.rds'))
object <- private('migrate_sceptre_object')(fixture)
record <- object@results@records$run_discovery_analysis
stopifnot(private('result_record_is_current')(object, 'run_discovery_analysis'))
# The bounded synthetic graph is explicitly not a scientific result. It uses
# 500 rows x 128 numeric draws (64,000 values, about 512 KiB before R overhead),
# a nested data.table attribute and a frozen genuine RunPlan to exercise ownership.
synthetic <- withr::with_seed(19283L, {
    table <- data.table::data.table(response_id = paste0('example_', seq_len(500L)),
        p_value = seq_len(500L) / 501,
        resamples = lapply(seq_len(500L), function(i) rnorm(128L)))
    data.table::setattr(table, 'diagnostics', list(table = data.table::data.table(n = 500L)))
    list(table = table, provenance = list(plan = record$provenance$plan), history = list())
})
inputs <- readRDS(file.path(refactor_fixture_dir(), 'inputs.rds'))
case <- readRDS(file.path(refactor_fixture_dir(), 'cases.rds'))$high_complement_crt
prepared <- refactor_capture(refactor_prepare(inputs, case))$value
run <- function(value) withr::with_seed(172839L,
    refactor_capture(sceptre::run_discovery_analysis(value, print_progress = FALSE, output_amount = 3L)))
cold <- run(prepared)
warm <- run(cold$value)
stopifnot(identical(refactor_plain(private('current_result')(cold$value, 'run_discovery_analysis')),
                    refactor_plain(private('current_result')(warm$value, 'run_discovery_analysis'))),
          identical(cold$warnings, warm$warnings),
          identical(cold$value@fits, warm$value@fits))

values <- list(actual_five_row_table = record$table,
    actual_record_with_history = record,
    actual_crt_output3_record = cold$value@results@records$run_discovery_analysis,
    bounded_nested_draw_graph = synthetic)

profile_copy <- function(value) {
    profile <- tempfile('sceptre-copy-allocations-')
    on.exit({ Rprofmem(NULL); unlink(profile) }, add = TRUE)
    Rprofmem(profile)
    copy <- private('owned_result_value')(value)
    Rprofmem(NULL)
    lines <- readLines(profile, warn = FALSE)
    bytes <- suppressWarnings(as.numeric(sub(' .*$', '', lines)))
    list(total_allocated_bytes = sum(bytes, na.rm = TRUE), allocations = sum(!is.na(bytes)),
        copy_size_bytes = as.numeric(object.size(copy)))
}
rows <- list()
for (name in names(values)) {
    value <- values[[name]]
    invisible(private('owned_result_value')(value))
    times <- vapply(seq_len(20L), function(i) {
        started <- proc.time()[['elapsed']]
        for (j in seq_len(10L)) copy <- private('owned_result_value')(value)
        (proc.time()[['elapsed']] - started) / 10
    }, numeric(1))
    memory <- profile_copy(value)
    rows[[name]] <- data.frame(workload = name, input_size_bytes = as.numeric(object.size(value)),
        median_seconds = median(times), min_seconds = min(times), max_seconds = max(times),
        allocated_bytes = memory$total_allocated_bytes, allocations = memory$allocations,
        copy_size_bytes = memory$copy_size_bytes, samples = 20L, copies_per_sample = 10L)
}
copy_results <- do.call(rbind, rows)
rownames(copy_results) <- NULL
write.csv(copy_results, file.path(output, 'copy-boundaries.csv'), row.names = FALSE)

# Time the uninstrumented public workflow. The cold object is always the same
# prepared input; the warm object is the corresponding completed cached screen.
cache_times <- do.call(rbind, lapply(seq_len(3L), function(i) {
    cold_time <- system.time(invisible(run(prepared)))[['elapsed']]
    warm_time <- system.time(invisible(run(cold$value)))[['elapsed']]
    data.frame(repetition = i, cold_seconds = cold_time, warm_seconds = warm_time)
}))
write.csv(cache_times, file.path(output, 'crt-cache-times.csv'), row.names = FALSE)

# Count actual fit calls in a separate isolated probe; tracing is excluded from
# the timings above and performs no random operations.
counts <- new.env(parent = emptyenv())
counts$response <- counts$perturbation <- 0L
namespace <- asNamespace('sceptre')
on.exit(untrace('fit_response', where = namespace), add = TRUE)
trace('fit_response', where = namespace, print = FALSE,
    tracer = function() counts$response <- counts$response + 1L)
on.exit(untrace('fit_perturbation_glm', where = namespace), add = TRUE)
trace('fit_perturbation_glm', where = namespace, print = FALSE,
    tracer = function() counts$perturbation <- counts$perturbation + 1L)
traced_cold <- run(prepared)
cold_counts <- c(response = counts$response, perturbation = counts$perturbation)
counts$response <- counts$perturbation <- 0L
traced_warm <- run(traced_cold$value)
warm_counts <- c(response = counts$response, perturbation = counts$perturbation)
# The main() exit handler removes both traces even if a later check fails.
stopifnot(identical(refactor_plain(private('current_result')(traced_cold$value, 'run_discovery_analysis')),
                    refactor_plain(private('current_result')(traced_warm$value, 'run_discovery_analysis'))),
          all(warm_counts == 0L))
fit_counts <- data.frame(path = c('cold', 'warm'),
    response_fit_calls = c(cold_counts[['response']], warm_counts[['response']]),
    perturbation_fit_calls = c(cold_counts[['perturbation']], warm_counts[['perturbation']]))
write.csv(fit_counts, file.path(output, 'crt-fit-counts.csv'), row.names = FALSE)
plan <- cold$value@results@records$run_discovery_analysis$provenance$plan
stopifnot(identical(unname(source_hashes), unname(tools::md5sum(source_paths))))
write.csv(data.frame(path = names(source_hashes), md5 = unname(source_hashes)),
    file.path(output, 'source-hashes.csv'), row.names = FALSE)
metadata <- list(source_hash_algorithm = 'MD5', source_hashes = source_hashes,
    timestamp_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE),
    revision = system2('git', c('rev-parse', 'HEAD'), stdout = TRUE),
    package_path = find.package('sceptre'), source_or_library = if (nzchar(selected_library)) 'installed' else 'source pkgload',
    R = R.version.string, platform = R.version$platform,
    data_table = as.character(packageVersion('data.table')),
    copy_sampling = '20 batches of 10 copies; allocation profile is one separate copy; no baseline allocation subtracted',
    allocation_limit = 'Rprofmem total allocations, not peak resident memory; profiler metadata and timing noise may contribute',
    synthetic = list(rows = 500L, draws_per_row = 128L, draw_values = 64000L, scientific_result = FALSE),
    crt = list(source_cells = ncol(inputs$response_matrix), prepared_cells = length(prepared@preparation@analysis_cells),
        requested_pairs = nrow(inputs$discovery_pairs), B1 = plan@B1, B2 = plan@B2, B3 = plan@B3,
        exact_cold_warm_results = TRUE, exact_cold_warm_fits = TRUE,
        exact_cold_warm_warnings = TRUE, seeded_run = 172839L),
    limitations = 'Single local small fixture, three timing repetitions, warmed process/filesystem, no cross-platform or large-screen inference.')
saveRDS(list(metadata = metadata, copy = copy_results, cache = cache_times, fits = fit_counts),
    file.path(output, 'benchmark-results.rds'))
dput(metadata, file.path(output, 'metadata.R'))
print(copy_results)
print(cache_times)
print(fit_counts)

}
main()
