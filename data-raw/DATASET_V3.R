# example data from Xie dataset
library(dplyr)
library(ondisc)
library(magrittr)

# set the Gasperini directory
gasp_fp <- paste0(.get_config_path("LOCAL_GASPERINI_2019_DATA_DIR"), "at-scale/")
pairs_ungrouped <- readRDS(paste0(gasp_fp, "processed/gRNA_ungrouped/pairs_ungrouped.rds"))
n_cells <- 40000

#########################################
# 1. Construct gene-gRNA pairs data frame
#########################################
set.seed(62)
gene_gRNA_pairs_raw <- pairs_ungrouped %>%
  select(gene_id, gRNA_group, site_type) %>%
  distinct()
gene_gRNA_pairs_pc_candidate <- rbind(
  gene_gRNA_pairs_raw %>%
    filter(site_type == "selfTSS") %>%
    sample_n(5),
  gene_gRNA_pairs_raw %>%
    filter(site_type == "DHS") %>%
    sample_n(15)
)
my_gene_ids <- gene_gRNA_pairs_pc_candidate$gene_id %>%
  as.character() %>% unique()
neg_control_gRNAs <- gene_gRNA_pairs_raw %>%
  filter(site_type == "NTC") %>%
  pull(gRNA_group) %>% unique() %>%
  sample(., size = 5) %>% as.character()
gene_gRNA_pairs_NTC <- expand.grid(gene_id = my_gene_ids, gRNA_group = neg_control_gRNAs) %>%
  mutate(site_type = "NTC")
gene_gRNA_pairs <- rbind(gene_gRNA_pairs_pc_candidate,
                         gene_gRNA_pairs_NTC) %>%
  mutate_all(as.character) %>%
  mutate(pair_type = factor(site_type, c("selfTSS", "DHS", "NTC"),
                       c("positive_control", "candidate", "negative_control")),
         site_type = NULL)
my_gRNA_groups <- unique(as.character(gene_gRNA_pairs$gRNA_group))

#####################################
# 2. Construct gRNA groups data frame
#####################################
gRNA_groups_table <- pairs_ungrouped %>% select(gRNA_group, barcode, site_type) %>%
  distinct() %>% filter(gRNA_group %in% my_gRNA_groups) %>% filter(site_type != "TSS") %>%
  rename(gRNA_id = barcode, gRNA_type = site_type) %>%
  mutate(gRNA_type = factor(gRNA_type, c("DHS", "selfTSS", "NTC"), c("enh_target", "tss_target", "non_target"))) %>%
  mutate_all(.tbl = ., .funs = as.character) %>%
  arrange(gRNA_type) %>%
  select(gRNA_id, gRNA_group, gRNA_type)
my_gRNA_ids <- as.character(gRNA_groups_table$gRNA_id)

###########
# Load ODMs
###########
gRNA_odm <- read_odm(odm_fp = paste0(gasp_fp, "processed/gRNA_ungrouped/gasp_scale_gRNA_counts_ungrouped.odm"),
                     metadata_fp = paste0(gasp_fp, "processed/gRNA_ungrouped/gasp_scale_gRNA_metadata_ungrouped.rds"))

gene_odm <- read_odm(odm_fp = paste0(gasp_fp, "processed/gene/gasp_scale_gene_expressions.odm"),
                     metadata_fp = paste0(gasp_fp, "processed/gene/gasp_scale_gene_metadata.rds"))

###################################
# multimodal data; downsample cells
###################################
multimodal_odm <- multimodal_ondisc_matrix(covariate_ondisc_matrix_list = list(gRNA = gRNA_odm,
                                                                                gene = gene_odm))
multimodal_odm <- multimodal_odm[,(multimodal_odm %>% get_cell_covariates() %>% pull(gRNA_n_umis)) != 0]
cell_ids <- sample(seq(1, ncol(multimodal_odm)), n_cells)
multimodal_odm_downsample <- multimodal_odm[,cell_ids]

############################
# Simplify the cell barcodes
############################
cell_barcodes <- get_cell_barcodes(multimodal_odm_downsample)
my_cell_barcodes <- substr(x = cell_barcodes, start = 0, stop = 23)

##############################
# Load gRNA and gene data data
##############################
gRNA_sub <- get_modality(multimodal_odm_downsample, "gRNA")
gene_sub <- get_modality(multimodal_odm_downsample, "gene")

gene_matrix <- gene_sub[[my_gene_ids,]]
row.names(gene_matrix) <- my_gene_ids
colnames(gene_matrix) <- my_cell_barcodes

gRNA_matrix <- gRNA_sub[[my_gRNA_ids,]]
row.names(gRNA_matrix) <- my_gRNA_ids
colnames(gRNA_matrix) <- my_cell_barcodes

##############################
# Compute the covariate matrix
##############################
covariate_matrix <- multimodal_odm_downsample %>% get_cell_covariates() %>%
  summarize(lg_gRNA_lib_size = log(gRNA_n_umis),
         lg_gene_lib_size = log(gene_n_umis),
         p_mito = gene_p_mito,
         batch = gene_batch)
row.names(covariate_matrix) <- my_cell_barcodes

###############################
# Simplify the gRNA group names
###############################
simplify_grna_names <- function(vect) {
  vect <- gsub(x = vect, pattern = "_top_two", replacement = "")
  ntc_orig_names <- sort(unique(grep(pattern = "random*|scrambled*", x = vect, value = TRUE)))
  ntc_new_names <- paste0("NTC_", seq(1, length(ntc_orig_names)))
  for (i in seq(1, length(ntc_orig_names))) {
    orig_name <- ntc_orig_names[i]
    vect[orig_name ==  vect] <- ntc_new_names[i]
  }
  return(vect)
}

gene_gRNA_pairs$gRNA_group <- simplify_grna_names(gene_gRNA_pairs$gRNA_group )
gRNA_groups_table$gRNA_group <- simplify_grna_names(gRNA_groups_table$gRNA_group)
gene_gRNA_group_pairs <- gene_gRNA_pairs

usethis::use_data(gene_matrix, gRNA_matrix, covariate_matrix, gene_gRNA_group_pairs, gRNA_groups_table, overwrite = TRUE)
