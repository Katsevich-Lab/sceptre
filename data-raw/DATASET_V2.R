# example data from Xie dataset
library(ondisc)
library(magrittr)

# load the Xie odms corresponding to genes and gRNAs
xie_fp <- paste0(.get_config_path("LOCAL_XIE_2019_DATA_DIR"), "processed/")
glmeiv_fp <- paste0(.get_config_path("LOCAL_GLMEIV_DATA_DIR"), "public/xie/data/")

# Load gene and gRNA ODMs, alongside covariate matrix and gene-gRNA pair list
gene_odm <- read_odm(odm_fp = paste0(xie_fp, "gene/expression_matrix.odm"),
                     metadata_fp = paste0(xie_fp, "gene/metadata.rds"))
gRNA_odm_ungrouped <- read_odm(odm_fp = paste0(xie_fp, "gRNA/raw_ungrouped.odm"),
                               metadata_fp = paste0(xie_fp, "gRNA/raw_ungrouped_metadata.rds"))

# compute covariate matrix
covariate_matrix <- gene_odm %>% get_cell_covariates() %>% dplyr::select(n_umis, p_mito, batch) %>%
  dplyr::mutate(n_gRNA_umis = (gRNA_odm_ungrouped %>% get_cell_covariates() %>% dplyr::pull(n_umis))) %>%
  dplyr::summarize(lg_gene_lib_size = log(n_umis), p_mito, batch, lg_gRNA_lib_size = log(n_gRNA_umis)) %>%
  dplyr::select(lg_gene_lib_size, lg_gRNA_lib_size, batch, p_mito)
row.names(covariate_matrix) <- get_cell_barcodes(gene_odm)

# Determine gRNAs to analyze, as well as groupings
set.seed(4)
gRNA_groups_df <- readRDS(paste0(xie_fp, "../intermediate/guide_seqs.rds"))
my_enh_regions <- sample(x = unique(gRNA_groups_df$hg38_enh_region), size = 5, replace = FALSE)
site_table <- gRNA_groups_df %>% dplyr::filter(hg38_enh_region %in% my_enh_regions) %>%
  dplyr::rename("site" = "hg38_enh_region", "gRNA_id" = "spacer_seq")
gRNAs_to_analyze <- gRNA_grps$gRNA_id %>% unique()

# determine genes and gRNAs to analyze
genes_to_analyze <- get_feature_covariates(gene_odm) %>% dplyr::filter(mean_expression > 1) %>%
  dplyr::sample_n(5) %>% row.names()

# load matrices corresponding to sampled genes and gRNAs
gene_matrix <- as.matrix(gene_odm[[genes_to_analyze,]])
row.names(gene_matrix) <- genes_to_analyze; colnames(gene_matrix) <- get_cell_barcodes(gene_odm)

gRNA_matrix <- gRNA_odm_ungrouped[[gRNAs_to_analyze,]]
row.names(gRNA_matrix) <- gRNAs_to_analyze; colnames(gRNA_matrix) <- get_cell_barcodes(gRNA_odm_ungrouped)

# create the data frame of gene-gRNA pairs to analyze
gene_gRNA_pairs <- expand.grid(gene_id = genes_to_analyze, gRNA_id = my_enh_regions) %>%
  dplyr::sample_n(23)

# save the expression matrix, perturbation matrix, covariate matrix, and gene-gRNA pairs matrix
usethis::use_data(gene_matrix, gRNA_matrix, covariate_matrix, gene_gRNA_pairs, site_table, overwrite = TRUE)
