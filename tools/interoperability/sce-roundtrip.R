# Execute in a fresh Rscript with the candidate package and optional SCE on R_LIBS.
# Usage: Rscript sce-roundtrip.R PIPELINE_CHECKOUT NEW_OUTPUT_DIRECTORY
# Pipeline checkout must contain W05's preserve-external-assignment companion.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, !dir.exists(args[2L]))
pipeline <- normalizePath(args[1L])
out <- args[2L]
dir.create(out, recursive = TRUE)
out <- normalizePath(out)
library(sceptre)
library(SingleCellExperiment)
data("lowmoi_example_data", package = "sceptre")
x <- lowmoi_example_data
# Use the complete bundled count library. Only supply missing cell IDs once,
# before constructing the SCE, using the same source ordering in both assays.
cells <- colnames(x$response_matrix)
if (is.null(cells)) cells <- paste0("cell_", seq_len(ncol(x$response_matrix)))
colnames(x$response_matrix) <- colnames(x$grna_matrix) <- cells
rownames(x$extra_covariates) <- cells
sce <- SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(x$response_matrix, sparse = TRUE)),
    colData = S4Vectors::DataFrame(x$extra_covariates))
altExp(sce, "gRNA") <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = Matrix::Matrix(x$grna_matrix, sparse = TRUE)))
metadata <- x$grna_target_data_frame
metadata <- metadata[match(rownames(x$grna_matrix), metadata$grna_id), ]
object <- import_from_sce(sce, metadata, "low",
                         covariates = names(x$extra_covariates))
targets <- unique(metadata$grna_target[metadata$grna_target != "non-targeting"])
pairs <- data.frame(grna_target = head(targets, 2L),
                    response_id = rep(rownames(sce)[1L], min(2L, length(targets))))
object <- set_analysis_parameters(object, discovery_pairs = pairs,
    formula_object = ~ log(response_n_umis), control_group = "complement",
    side = "left", resampling_approximation = "no_approximation")
# A transparent example external caller; replace with your independently
# supplied binary calls. Raw guide counts stay intact in the count library.
calls <- x$grna_matrix >= 5L
object <- set_grna_assignments(object, calls,
    "worked example threshold >=5; not a recommended universal caller")
config <- list(qc = list(response_n_umis_range = c(0, 1),
    response_n_nonzero_range = c(0, 1), p_mito_threshold = 1,
    n_nonzero_trt_thresh = 3L, n_nonzero_cntrl_thresh = 3L),
    execution = list(output_amount = 1L))
saveRDS(sce, file.path(out, "source-sce.rds"))
saveRDS(object, file.path(out, "assigned-object.rds"))
saveRDS(config, file.path(out, "config.rds"))
source(file.path(pipeline, "bin", "pipeline.R"))
lock <- package_identity()
saveRDS(lock, file.path(out, "lock.rds"))
create_input(object, file.path(out, "input"), config, lock)
# Launch each boundary in a fresh R process using the pipeline's actual CLI.
previous <- file.path(out, "input")
for (stage in c("parameters", "assignment", "qc", "discovery")) {
    destination <- file.path(out, stage)
    status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(
        file.path(pipeline, "bin", "stage.R"), file.path(out, "config.rds"),
        previous, destination, stage, "20261019", file.path(out, "lock.rds"), "")),
        stdout = file.path(out, paste0(stage, ".log")), stderr = TRUE)
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
returned <- export_results_to_sce(final, sce)
saveRDS(returned, file.path(out, "results-sce.rds"))
stopifnot(identical(assays(returned), assays(sce)),
          is.null(S4Vectors::metadata(sce)$sceptre))
writeLines(c(capture.output(sessionInfo()),
    paste("pipeline_commit", system2("git", c("-C", shQuote(pipeline), "rev-parse", "HEAD"), stdout = TRUE)),
    paste("package_archive_identity", lock$installed_sha256)), file.path(out, "session.txt"))
cat("PASS: full-library SCE -> external assignments -> fresh-process pipeline -> owned SCE results; direct equality tolerance 0\n")
