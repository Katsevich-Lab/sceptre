gasperini_dir <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/processed/")

# 0. set the Gasperini directory load the group, gene, and grna data
gasp_fp <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/processed/")
pairs_grouped <- readRDS(paste0(gasp_fp, "/pairs_grouped.rds")) |> dplyr::distinct()

gene_odm_fp <- paste0(gasp_fp, "gene/matrix.odm")
gene_metadata_fp <- paste0(gasp_fp, "gene/metadata.rds")
gene_odm <- ondisc::read_odm(odm_fp = gene_odm_fp, metadata_fp = gene_metadata_fp)

grna_odm_fp <- paste0(gasp_fp, "grna_expression/matrix.odm")
grna_metadata_fp <-  paste0(gasp_fp, "grna_expression/metadata.rds")
grna_odm <- ondisc::read_odm(odm_fp = grna_odm_fp, metadata_fp = grna_metadata_fp)
grna_feature_df <- grna_odm |>
  ondisc::get_feature_covariates()
grna_feature_df <- grna_feature_df |>
  dplyr::mutate(grna_id = rownames(grna_feature_df)) |>
  `rownames<-`(NULL)

# 1. Construct the discovery pairs data frame, which contains both cis and positive control pairs
set.seed(4)
my_cis_grna_groups <- c(pairs_grouped |>
  dplyr::filter(site_type == "cis") |>
  dplyr::pull(grna_group) |>
  sample(12), "chr8.847_top_two", "chr9.1594_top_two", "chr9.2869_top_two", "chr9.3633_top_two", "chr9.871_top_two")
cis_pairs <- pairs_grouped |>
  dplyr::filter(grna_group %in% my_cis_grna_groups)
my_pc_grna_groups <- pairs_grouped |>
  dplyr::filter(site_type == "pos_cntrl") |>
  dplyr::pull(grna_group) |>
  sample(10)
pc_pairs <- pairs_grouped |>
  dplyr::filter(grna_group %in% my_pc_grna_groups, site_type == "pos_cntrl")
discovery_pairs <- rbind(cis_pairs, pc_pairs) |>
  dplyr::rename(type = site_type, response_id = gene_id)

# 2. sample some set of negative control pairs
nt_grnas <- grna_feature_df |>
  dplyr::filter(target_type == "non-targeting") |>
  dplyr::sample_n(25) |> dplyr::pull(grna_id)

# 3. construct the grna group table
grna_group_data_frame_highmoi <- rbind(grna_feature_df |>
  dplyr::filter(grna_group %in% discovery_pairs$grna_group) |>
  dplyr::select(grna_id, grna_group),
  data.frame(grna_id = nt_grnas, grna_group = "non-targeting")) |>
  dplyr::arrange(grna_group)

# 4. construct the respomse and grna matrices; downsample cells
multimodal_odm <- ondisc::multimodal_ondisc_matrix(covariate_ondisc_matrix_list = list(grna = grna_odm, gene = gene_odm))
cell_ids <- sample(seq(1, dim(multimodal_odm)[[1]][2]), 40002)
multimodal_odm_downsample <- multimodal_odm[,cell_ids]

# 5. get the gene matrix
gene_sub <- ondisc::get_modality(multimodal_odm_downsample, "gene")
my_gene_ids <- unique(discovery_pairs$response_id)
gene_matrix <- gene_sub[[my_gene_ids,]]
rownames(gene_matrix) <- my_gene_ids

# 6. get the grna matrix
grna_sub <- ondisc::get_modality(multimodal_odm_downsample, "grna")
my_grna_ids <- unique(grna_group_data_frame_highmoi$grna_id)
grna_matrix <- grna_sub[[my_grna_ids,]]
rownames(grna_matrix) <- my_grna_ids

# 7. Compute the covariate matrix
covariate_matrix <- multimodal_odm_downsample |>
  ondisc::get_cell_covariates() |>
  dplyr::select(batch = gene_batch) |>
  `rownames<-`(NULL)

# 8. sort according to batch; also remove cells containing no gRNAs
cell_order <- order(covariate_matrix$batch)
response_matrix_highmoi <- gene_matrix[,cell_order]
grna_matrix_highmoi <- grna_matrix[,cell_order]
covariate_data_frame_highmoi <- covariate_matrix[cell_order,,drop=FALSE]

# 9. rename the data objects
discovery_pairs_highmoi <- discovery_pairs |> dplyr::filter(type == "cis")
pc_pairs_highmoi <- discovery_pairs |> dplyr::filter(type == "pos_cntrl")

response_matrix_highmoi_experimental <- response_matrix_highmoi
grna_matrix_highmoi_experimental <- grna_matrix_highmoi
extra_covariates_highmoi_experimental <- covariate_data_frame_highmoi
grna_group_data_frame_highmoi_experimental <- grna_group_data_frame_highmoi
discovery_pairs_highmoi_experimental <- discovery_pairs_highmoi
pc_pairs_highmoi_experimental <- pc_pairs_highmoi

usethis::use_data(response_matrix_highmoi_experimental, grna_matrix_highmoi_experimental,
                  extra_covariates_highmoi_experimental, grna_group_data_frame_highmoi_experimental,
                  discovery_pairs_highmoi_experimental, pc_pairs_highmoi_experimental, overwrite = TRUE)
