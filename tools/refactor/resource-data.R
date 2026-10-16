#!/usr/bin/env Rscript
# Deliberate benchmark input generation; no expected scientific outputs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: resource-data.R <fixtures> <output-dir>")
fixtures <- normalizePath(args[[1L]])
directory <- args[[2L]]
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
library(Matrix)
source(file.path(dirname(dirname(fixtures)), "helper-refactor.R"))
set.seed(20261016)

sparse_responses <- function(n_cells, n_responses, n_active, difficult) {
    batches <- rep(seq_len(4L), length.out = n_cells)
    z <- sin(seq_len(n_cells) / 317)
    columns <- vector("list", n_responses)
    values <- vector("list", n_responses)
    for (r in seq_len(n_responses)) {
        if (r <= n_active) {
            baseline <- c(0.15, 0.4, 1.5, 3)[(r - 1L) %% 4L + 1L]
            size <- c(0.1, 0.15, 0.3, 1)[(r - 1L) %% 4L + 1L]
            mu <- baseline * exp(0.7 * z + c(-0.8, 0.2, 0.4, 0)[batches])
            if (difficult && r %% 16L == 0L) mu[batches != 1L] <- 0
            y <- rnbinom(n_cells, mu = mu, size = size)
            columns[[r]] <- which(y != 0)
            values[[r]] <- y[columns[[r]]]
        } else {
            columns[[r]] <- sort(sample.int(n_cells, max(20L, round(n_cells * 0.01))))
            values[[r]] <- 1L + rgeom(length(columns[[r]]), 0.6)
        }
    }
    sparseMatrix(i = rep(seq_len(n_responses), lengths(columns)),
                 j = unlist(columns), x = unlist(values),
                 dims = c(n_responses, n_cells),
                 dimnames = list(paste0("r", seq_len(n_responses)), paste0("c", seq_len(n_cells))),
                 repr = "R")
}

synthetic_case <- function(n_cells, n_responses, n_active, mechanism, backed, difficult) {
    response <- sparse_responses(n_cells, n_responses, n_active, difficult)
    n_guides <- if (mechanism == "crt") 24L else if (backed) 256L else 120L
    guide_ids <- paste0("g", seq_len(n_guides))
    if (mechanism == "crt") {
        membership <- lapply(seq_len(n_guides), function(g) {
            z <- sin(seq_len(n_cells) / 317)
            which(runif(n_cells) < plogis(-3.9 + 0.5 * z))
        })
        guides <- sparseMatrix(i = rep(seq_len(n_guides), lengths(membership)),
                               j = unlist(membership), x = 50,
                               dims = c(n_guides, n_cells), repr = "R")
    } else {
        # Only eight requested target groups; the rest are NT. This keeps the
        # production draw-bank payload bounded without reducing its B budget.
        guides <- sparseMatrix(i = rep(seq_len(n_guides), length.out = n_cells),
                               j = seq_len(n_cells), x = 50,
                               dims = c(n_guides, n_cells), repr = "R")
    }
    dimnames(guides) <- list(guide_ids, colnames(response))
    covariates <- data.frame(batch = factor(rep(seq_len(4L), length.out = n_cells)),
                             z = sin(seq_len(n_cells) / 317))
    if (difficult) covariates$z_near <- covariates$z + 0.01 * cos(seq_len(n_cells) / 211)
    formula <- if (difficult) ~ batch + z + z_near else ~ batch + z
    list(
        data = list(response_matrix = response, grna_matrix = guides,
                    grna_target_data_frame = data.frame(grna_id = guide_ids,
                        grna_target = c(paste0("T", seq_len(8L)), rep("non-targeting", n_guides - 8L))),
                    extra_covariates = covariates),
        parameters = list(moi = if (mechanism == "crt") "high" else "low",
                          mechanism = mechanism, control = "complement", backed = backed,
                          formula = formula, assignment_threshold = 40,
                          discovery_pairs = expand.grid(grna_target = paste0("T", seq_len(8L)),
                              response_id = rownames(response)[seq_len(n_active)], stringsAsFactors = FALSE)),
        recipe = list(kind = "deterministic synthetic; not actual Replogle data", seed = 20261016L,
                      cells = n_cells, responses = n_responses, active_responses = n_active,
                      response_nonzero = length(response@x), response_density = length(response@x) / prod(dim(response)),
                      dense_response_bytes = prod(as.double(dim(response))) * 8,
                      declared_dense_assay_budget_bytes = if (backed) 256 * 1024^2 else NULL,
                      budget_is_os_enforced = FALSE,
                      dispersion_sizes = c(0.1, 0.15, 0.3, 1),
                      difficult = if (difficult) "Every16th active response expressed only in batch1; near-collinear z/z_near covariates" else "overdispersed negative-binomial counts",
                      design_condition_number = kappa(model.matrix(formula, covariates)),
                      guide_assignment = if (mechanism == "crt") "overlapping logistic memberships" else "balanced singleton memberships")
    )
}

write_case <- function(name, object) {
    environment(object$parameters$formula) <- baseenv()
    path <- file.path(directory, paste0(name, ".rds"))
    saveRDS(object, path, compress = "gzip")
    list(name = name, path = basename(path), sha256 = digest::digest(file = path, algo = "sha256"),
         recipe = object$recipe, parameters = c(object$parameters[c("moi", "mechanism", "control", "backed", "assignment_threshold")],
             list(formula = paste(deparse(object$parameters$formula), collapse = ""),
                  requested_pairs = nrow(object$parameters$discovery_pairs))))
}

inputs <- readRDS(file.path(fixtures, "backed-inputs.rds"))
compact <- list(
    data = list(response_matrix = inputs$response_matrix, grna_matrix = inputs$low_grna_matrix,
                grna_target_data_frame = inputs$grna_target_data_frame, extra_covariates = inputs$extra_covariates),
    parameters = list(moi = "low", mechanism = "permutations", control = "complement", backed = FALSE,
                      formula = ~ batch + z, assignment_threshold = 40, discovery_pairs = inputs$discovery_pairs),
    recipe = list(kind = "RF0 compact synthetic", cells = ncol(inputs$response_matrix), responses = nrow(inputs$response_matrix),
                  source = "backed-inputs.rds", source_sha256 = digest::digest(file = file.path(fixtures, "backed-inputs.rds"), algo = "sha256"))
)
manifest <- list(write_case("compact_perm", compact))
rm(compact, inputs)

data("lowmoi_example_data", package = "sceptredata")
actual <- lowmoi_example_data
targets <- setdiff(unique(actual$grna_target_data_frame$grna_target), "non-targeting")[seq_len(6L)]
real_case <- list(
    data = actual,
    parameters = list(moi = "low", mechanism = "permutations", control = "complement", backed = FALSE,
        formula = ~ log(response_n_umis) + log(response_n_nonzero) + bio_rep,
        assignment_threshold = 5, discovery_pairs = expand.grid(grna_target = targets,
            response_id = rownames(actual$response_matrix)[seq_len(48L)], stringsAsFactors = FALSE)),
    recipe = list(kind = "actual Papalexi2021 example counts from installed sceptredata; downsampled genes",
                  cells = ncol(actual$response_matrix), responses = nrow(actual$response_matrix),
                  sceptredata_version = as.character(packageVersion("sceptredata")),
                  response_nonzero = length(actual$response_matrix@x),
                  source = "sceptredata::lowmoi_example_data")
)
manifest[[length(manifest) + 1L]] <- write_case("papalexi_perm", real_case)
rm(actual, real_case, lowmoi_example_data)

manifest[[length(manifest) + 1L]] <- write_case("sparse_difficult_perm",
    synthetic_case(12000L, 768L, 96L, "permutations", FALSE, TRUE))
manifest[[length(manifest) + 1L]] <- write_case("shared_crt",
    synthetic_case(12000L, 192L, 64L, "crt", FALSE, FALSE))
manifest[[length(manifest) + 1L]] <- write_case("backed_perm_512mib",
    synthetic_case(16384L, 4096L, 128L, "permutations", TRUE, FALSE))
jsonlite::write_json(list(schema_version = 1L, generator = "resource-data.R", R = R.version.string,
                          cases = manifest), file.path(directory, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
