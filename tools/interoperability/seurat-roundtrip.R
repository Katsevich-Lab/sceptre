# Usage: Rscript seurat-roundtrip.R PIPELINE_CHECKOUT NEW_OUTPUT_DIRECTORY
# Install candidate sceptre and optional SeuratObject >=5; SCE is not needed.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, !dir.exists(args[2L]))
pipeline <- normalizePath(args[1L])
dir.create(args[2L], recursive = TRUE)
out <- normalizePath(args[2L])
library(sceptre)
library(SeuratObject)
data("lowmoi_example_data", package = "sceptre")
x <- lowmoi_example_data
metadata <- x$grna_target_data_frame
# Seurat constructors replace '_' in feature IDs. Make the mapping explicit
# BEFORE constructing the object; the adapter itself never guesses this mapping.
metadata$original_grna_id <- metadata$grna_id
ids <- paste0("guide", seq_len(nrow(x$grna_matrix)))
rownames(x$grna_matrix) <- metadata$grna_id <- ids
s <- CreateSeuratObject(counts = Matrix::Matrix(x$response_matrix, sparse = TRUE),
                        min.cells = 0, min.features = 0,
                        meta.data = x$extra_covariates)
s[["gRNA"]] <- CreateAssay5Object(counts = Matrix::Matrix(x$grna_matrix, sparse = TRUE),
                                  min.cells = 0, min.features = 0)
object <- import_from_seurat(s, metadata, "low", covariates = "batch")
targets <- unique(metadata$grna_target[metadata$grna_target != "non-targeting"])
pairs <- data.frame(grna_target = head(targets, 2L),
                    response_id = rep(rownames(x$response_matrix)[1L], 2L))
object <- set_analysis_parameters(object, discovery_pairs = pairs,
    formula_object = ~ log(response_n_umis), control_group = "complement",
    side = "left", resampling_approximation = "no_approximation")
object <- set_grna_assignments(object, x$grna_matrix >= 5L,
    "Seurat worked example raw UMI >=5; not a universal caller recommendation")
config <- list(qc = list(response_n_umis_range = c(0, 1),
    response_n_nonzero_range = c(0, 1), p_mito_threshold = 1,
    n_nonzero_trt_thresh = 3L, n_nonzero_cntrl_thresh = 3L),
    execution = list(output_amount = 1L))
saveRDS(s, file.path(out, "source-seurat.rds"))
saveRDS(object, file.path(out, "assigned-object.rds"))
saveRDS(config, file.path(out, "config.rds"))
saveRDS(metadata, file.path(out, "explicit-guide-id-map.rds"))
source(file.path(pipeline, "bin", "pipeline.R"))
lock <- package_identity()
saveRDS(lock, file.path(out, "lock.rds"))
create_input(object, file.path(out, "input"), config, lock)
previous <- file.path(out, "input")
for (stage in c("parameters", "assignment", "qc", "discovery")) {
    destination <- file.path(out, stage)
    status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(
        file.path(pipeline, "bin", "stage.R"), file.path(out, "config.rds"),
        previous, destination, stage, "20261019", file.path(out, "lock.rds"), "")),
        stdout = file.path(out, paste0(stage, ".log")),
        stderr = file.path(out, paste0(stage, ".stderr.log")))
    stopifnot(status == 0L)
    previous <- destination
}
final <- read_bundle(previous, config, lock)$object
stopifnot(identical(get_cell_covariates(final), get_cell_covariates(object)),
          identical(get_grna_assignments(final), get_grna_assignments(object)))
direct <- do.call(run_qc, c(list(object), config$qc))
set.seed(20261019)
direct <- run_discovery_analysis(direct, print_progress = FALSE)
stopifnot(isTRUE(all.equal(get_result(final, "run_discovery_analysis"),
                          get_result(direct, "run_discovery_analysis"), tolerance = 0)))
before <- serialize(s, NULL)
returned <- export_results_to_seurat(final, s)
stopifnot(identical(serialize(s, NULL), before),
          identical(returned[["RNA"]], s[["RNA"]]),
          identical(returned[["gRNA"]], s[["gRNA"]]),
          identical(returned[[]], s[[]]))
saveRDS(returned, file.path(out, "results-seurat.rds"))
writeLines(c(capture.output(sessionInfo()),
    paste("pipeline_commit", system2("git", c("-C", shQuote(pipeline), "rev-parse", "HEAD"), stdout = TRUE)),
    paste("installed_package_identity", lock$installed_sha256)), file.path(out, "session.txt"))
cat("PASS: raw Seurat -> external calls -> fresh-process pipeline -> owned Seurat results; direct equality tolerance 0\n")
