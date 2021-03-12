load_all(helpers = FALSE)
library(magrittr)
library(fst)
library(bigstatsr)
library(dplyr)
library(Matrix)
library(readr)
library(ondisc)

load_fbm <- function(fbm_metadata) {
  out <- FBM(nrow = fbm_metadata$nrow, ncol = fbm_metadata$ncol, type = fbm_metadata$type, backingfile = fbm_metadata$backingfile, is_read_only = TRUE, create_bk = FALSE)
  return(out)
}

args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE-manuscript/SCEPTRE/" else args[1]
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_gasp/file_paths_to_dirs.R"))
source(paste0(code_dir, "sceptre_paper/analysis_drivers/analysis_drivers_gasp/sceptre_function_args.R"))

### First: an example gene and gRNA, stored in rda format for easy access

gene_id <- "ENSG00000164713"
gRNA_id <- "chr7.3255_top_two"
# expression data
expressions <- cell_gene_expression_matrix[,which(ordered_gene_ids == gene_id)]
gRNA_indics <- read_fst(path = gRNA_indicator_matrix_fp, gRNA_id) %>% pull() %>% as.integer()
# subset according to cell id
expressions <- expressions[cell_subset]
gRNA_indicators <- gRNA_indics[cell_subset]
covariate_matrix <- covariate_matrix[cell_subset,]
covariate_matrix$prep_batch <- factor(covariate_matrix$prep_batch)
# save the data via use_this
usethis::use_data(expressions, gRNA_indicators, covariate_matrix, overwrite = TRUE)

### Next: example ondisc_matrix objects representing several genes and gRNAs

expression_nos <- c(2, 3, 5, 8, 10)
expressions <- cell_gene_expression_matrix[,expression_nos]
sparse_expressions <- t(Matrix(as.matrix(expressions), sparse = TRUE))
expression_ids <- ordered_gene_ids[expression_nos]

gRNA_ids <- c("chr11.2465_top_two", "chr10.1153_top_two", "chr1.10869_top_two", "chr1.12987_top_two")
gRNA_indicators <- read_fst(path = gRNA_indicator_matrix_fp, columns = gRNA_ids)
sparse_gRNA_indicators <- t(Matrix(as.matrix(gRNA_indicators), sparse = TRUE))

# copy cell barcodes
system(command = paste0("cp ", raw_data_dir, "/GSE120861_at_scale_screen.cells.txt ", tempdir()))
barcodes_fp <- paste0(tempdir(), "/GSE120861_at_scale_screen.cells.txt")
# write .mtx files for expressions and perturbations
exp_fp <- paste0(tempdir(), "/expressions.mtx")
writeMM(obj = sparse_expressions[,cell_subset], file = exp_fp)
pert_fp <- paste0(tempdir(), "/perturbations.mtx")
writeMM(obj = sparse_gRNA_indicators[,cell_subset], file = pert_fp)
# write tsv files for grna and gene ids
grnas_fp <- paste0(tempdir(), "/grna_ids.tsv")
write_tsv(x = data.frame(gRNA_ids), col_names = FALSE, file = grnas_fp)
genes_fp <- paste0(tempdir(), "/gene_ids.tsv")
write_tsv(x = data.frame(expression_ids), col_names = FALSE, file = genes_fp)
# set the extdata dir
extdata_dir <- "./inst/extdata"

expression_odm <- create_ondisc_matrix_from_mtx(mtx_fp = exp_fp, barcodes_fp = barcodes_fp, features_fp = genes_fp, on_disk_dir = extdata_dir, file_name = "expressions")[["ondisc_matrix"]]
perturbation_odm <- create_ondisc_matrix_from_mtx(mtx_fp = pert_fp, barcodes_fp = barcodes_fp, features_fp = grnas_fp, on_disk_dir = extdata_dir, file_name = "perturbations")[["ondisc_matrix"]]

# finally, save the covariate matrix and gene-gRNA pairs to extdata
saveRDS(object = covariate_matrix, file = paste0(extdata_dir, "/covariate_matrix.rds"))
gene_gRNA_pairs <- expand.grid(gene_id = expression_ids, gRNA_id = gRNA_ids)
saveRDS(object = gene_gRNA_pairs, file = paste0(extdata_dir, "/gene_gRNA_pairs.rds"))
