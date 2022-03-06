# example data from Xie dataset
library(ondisc)
library(magrittr)

# load the Xie odms corresponding to genes and gRNAs
xie_fp <- paste0(.get_config_path("LOCAL_XIE_2019_DATA_DIR"), "processed/")
glmeiv_fp <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/xie/data/")

# Load gene and gRNA ODMs, alongside covariate matrix and gene-gRNA pair list
gene_odm <- read_odm(odm_fp = paste0(xie_fp, "gene/expression_matrix.odm"),
                     metadata_fp = paste0(glmeiv_fp, "gene_metadata.rds"))
gRNA_odm <- read_odm(odm_fp = paste0(xie_fp, "gRNA/raw_grouped.odm"),
                     metadata_fp = paste0(glmeiv_fp, "gRNA_metadata.rds"))
raw_covariate_matrix <- readRDS(paste0(glmeiv_fp, "covariate_matrix.rds"))
lg_gRNA_lib_size <- readRDS(paste0(glmeiv_fp, "g_offset.rds"))
lg_gene_lib_size <- readRDS(paste0(glmeiv_fp, "m_offset.rds"))
covariate_matrix <- dplyr::mutate(raw_covariate_matrix, lg_gene_lib_size = lg_gene_lib_size, lg_gRNA_lib_size = lg_gRNA_lib_size)

# determine pairs (randomy sample 5 genes and 25 gRNAs)
set.seed(4)
gene_ids <- get_feature_ids(gene_odm)
gRNA_ids <- get_feature_ids(gRNA_odm)
my_genes <- sample(x = gene_ids, size = 5, replace = FALSE)
my_gRNAs <- sample(x = gRNA_ids, size = 25, replace = FALSE)

# load matrices corresponding to sampled genes and gRNAs
gene_matrix <- as.matrix(gene_odm[[my_genes,]])
gRNA_matrix <- gRNA_odm[[my_gRNAs,]]

# create the list of gene-gRNA pairs to analyze
gene_gRNA_pairs <- expand.grid(gene_id = my_genes, gRNA_id = my_gRNAs) %>%
  dplyr::sample_n(100)

# save the expression matrix, perturbation matrix, covariate matrix, and gene-gRNA pairs matrix
usethis::use_data(gene_matrix, gRNA_matrix, covariate_matrix, gene_gRNA_pairs, overwrite = TRUE)

