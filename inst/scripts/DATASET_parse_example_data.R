# Generates the example Parse Biosciences data in
# inst/extdata/parse_example/. The data are a small deterministic subset
# of Parse's public "Technical Performance of CRISPR Detect in Cell
# Lines" demo (split-pipe v1.0.6p), written in Parse's native split-pipe
# DGE format so they exercise import_data_from_parse().
#
# Source (CC BY 4.0, attribution to Parse Biosciences):
#   https://www.parsebiosciences.com/datasets/
#   technical-performance-of-crispr-detect-in-cell-lines/

library(Matrix)

# pin the sampler so the subset is reproducible across R versions
# (sample()'s algorithm changed in R 3.6); "Rejection" is the R >= 3.6
# default, so this matches the committed files.
RNGkind(sample.kind = "Rejection")
set.seed(4)
n_cells <- 5000L
n_genes <- 266L
base_url <- "https://cdn.parsebiosciences.com/crispr_detect/v1_0_6p"
out_dir <- file.path("inst", "extdata", "parse_example")

# ---- download the relevant demo files --------------------------------
dl <- function(path) {
    dest <- file.path(tempdir(), gsub("/", "_", path))
    download.file(file.path(base_url, path), dest, mode = "wb")
    dest
}
wt_mtx_fp <- dl("wt_library/DGE.mtx")
wt_genes_fp <- dl("wt_library/all_genes.csv")
crispr_mtx_fp <- dl("crispr_library/DGE.mtx")
crispr_genes_fp <- dl("crispr_library/all_genes.csv")

# ---- read (Parse stores matrices as cells x features) ----------------
wt_mat <- readMM(wt_mtx_fp)             # 51242 cells x 58395 genes
crispr_mat <- readMM(crispr_mtx_fp)     # 51242 cells x 3 guides
wt_genes <- read.csv(wt_genes_fp, colClasses = "character")
crispr_genes <- read.csv(crispr_genes_fp, colClasses = "character")
# the wt and crispr libraries share the same cells in the same row order,
# and each all_genes.csv is row-aligned with its matrix columns
stopifnot(
    nrow(wt_mat) == nrow(crispr_mat),
    ncol(wt_mat) == nrow(wt_genes),
    ncol(crispr_mat) == nrow(crispr_genes)
)

# ---- choose a deterministic subset of cells and genes ----------------
# one cell-index vector subsets both matrices, keeping them aligned
cell_idx <- sort(sample.int(nrow(wt_mat), n_cells))
wt_sub <- wt_mat[cell_idx, ]
# keep genes detected in the chosen cells, then sample n_genes of them
detected <- which(colSums(wt_sub) > 0)
stopifnot(length(detected) >= n_genes)
gene_idx <- sort(sample(detected, n_genes))

gene_mat <- wt_sub[, gene_idx]
grna_mat <- crispr_mat[cell_idx, ]
genes_out <- wt_genes[gene_idx, ]

# ---- write the four files in Parse's native format -------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_parse_mtx <- function(mat, fp) {
    trip <- as(drop0(mat), "TsparseMatrix")
    o <- order(trip@i, trip@j)
    con <- file(fp, "w")
    on.exit(close(con))
    writeLines("%%MatrixMarket matrix coordinate integer general", con)
    writeLines(sprintf("%%Rows=cells (%d), Cols=genes (%d)",
                       nrow(mat), ncol(mat)), con)
    writeLines(sprintf("%d %d %d", nrow(mat), ncol(mat),
                       length(trip@x)), con)
    writeLines(sprintf("%d %d %d", trip@i[o] + 1L, trip@j[o] + 1L,
                       as.integer(trip@x[o])), con)
}

write_parse_mtx(gene_mat, file.path(out_dir, "gene_mat.mtx"))
write_parse_mtx(grna_mat, file.path(out_dir, "grna_mat.mtx"))
write.csv(genes_out, file.path(out_dir, "all_genes.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(crispr_genes, file.path(out_dir, "all_grnas.csv"),
          row.names = FALSE, quote = FALSE)
