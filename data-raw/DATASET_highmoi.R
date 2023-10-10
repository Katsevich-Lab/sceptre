gasperini_dir <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/processed/")
gasperini_dir_raw <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/raw/")

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
set.seed(6) # 4
my_cis_grna_groups <- c(pairs_grouped |>
  dplyr::filter(site_type == "cis") |>
  dplyr::pull(grna_group) |> # 12
  sample(20), "chr8.847_top_two", "chr9.1594_top_two", "chr9.2869_top_two", "chr9.3633_top_two", "chr9.871_top_two")
cis_pairs <- pairs_grouped |>
  dplyr::filter(grna_group %in% my_cis_grna_groups, gene_id %in% gene_table$gene_id)
my_pc_grna_groups <- pairs_grouped |>
  dplyr::filter(site_type == "pos_cntrl") |>
  dplyr::pull(grna_group) |>
  sample(10)
pc_pairs <- pairs_grouped |>
  dplyr::filter(grna_group %in% my_pc_grna_groups, site_type == "pos_cntrl",
                gene_id %in% gene_table$gene_id)
targeting_pairs <- rbind(cis_pairs, pc_pairs) |>
  dplyr::rename(type = site_type, response_id = gene_id)

# 2. sample some set of negative control pairs
nt_grnas <- grna_feature_df |>
  dplyr::filter(target_type == "non-targeting") |>
  dplyr::sample_n(25) |> dplyr::pull(grna_id)

# 3. construct the grna group table
grna_group_data_frame_highmoi <- rbind(grna_feature_df |>
  dplyr::filter(grna_group %in% targeting_pairs$grna_group) |>
  dplyr::select(grna_id, grna_group),
  data.frame(grna_id = nt_grnas, grna_group = "non-targeting")) |>
  dplyr::arrange(grna_group)

# 4. add chromosomal locations to the grna group table
grna_loc_info <- readr::read_tsv(file = paste0(gasperini_dir_raw, "/GSE120861_gene_gRNAgroup_pair_table.at_scale.txt"),
                                 col_names = c("chr", "start", "end", "grna_group", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8"),
                                 skip = 1) |>
  dplyr::select(chr, start, end, grna_group) |>
  dplyr::distinct()
targeting_grna_grps <- grna_group_data_frame_highmoi |> dplyr::filter(grna_group != "non-targeting") |> dplyr::pull()
target_group_loc_info <- grna_loc_info |> dplyr::filter(grna_group %in% targeting_grna_grps) |>
  dplyr::mutate(start = as.integer(start), end = as.integer(end))
grna_group_data_frame_highmoi <- dplyr::left_join(x = grna_group_data_frame_highmoi,
                 y = target_group_loc_info, by = c("grna_group"))

# 5. rename some of the grna groups
renamed_grna_group <- pc_pairs$gene_id[match(x = grna_group_data_frame_highmoi$grna_group, table = pc_pairs$grna_group)]
grna_group_data_frame_highmoi$grna_group[!is.na(renamed_grna_group)] <- renamed_grna_group[!is.na(renamed_grna_group)]
enh_group_idxs <- grep(pattern = "chr", x = grna_group_data_frame_highmoi$grna_group)
enh_groups <- grna_group_data_frame_highmoi$grna_group[enh_group_idxs]
grna_group_data_frame_highmoi$grna_group[enh_group_idxs] <- factor(x = enh_groups, levels = unique(enh_groups),
                                                                   labels = paste0("candidate_enh_", seq_along(unique(enh_groups)))) |>
  as.character()

# 6. reset the start and end positions
grna_group_data_frame_highmoi <- dplyr::group_by(grna_group_data_frame_highmoi, grna_group) |>
  dplyr::group_modify(.f = function(tbl, key) {
    if (!all(is.na(tbl$start))) {
      group_start <- min(tbl$start)
      group_end <- max(tbl$end)
      if (group_end - group_start == 1L) group_end <- group_end + 30
      posits <- as.integer(floor(seq(group_start, group_end, length.out = nrow(tbl) + 1)))
      tbl$start <- posits[seq(1, nrow(tbl))]
      tbl$end <- posits[seq(2, nrow(tbl) + 1)]
    }
    return(tbl)
  }) |> dplyr::relocate(grna_id)

# 7. determine cells that contain at least one gRNA
my_grna_ids <- unique(grna_group_data_frame_highmoi$grna_id)
grna_matrix <- grna_odm[[my_grna_ids,]]
cell_ids <- which(apply(X = as.matrix(grna_matrix >= 5), 2, any))

# 8. construct the response and grna matrices; downsample cells
multimodal_odm <- ondisc::multimodal_ondisc_matrix(covariate_ondisc_matrix_list = list(grna = grna_odm, gene = gene_odm))
multimodal_odm_downsample <- multimodal_odm[,cell_ids]

# 9. get the gene matrix, adding a couple MT genes for good measure
gene_sub <- ondisc::get_modality(multimodal_odm_downsample, "gene")
all_gene_ids <- ondisc::get_feature_ids(gene_sub)
mt_genes <- gene_table |> dplyr::filter(chr == "chrM", gene_id %in% all_gene_ids) |> dplyr::pull(gene_id)
my_gene_ids <- c(unique(targeting_pairs$response_id), sample(mt_genes, 2))

gene_matrix <- gene_sub[[my_gene_ids,]]
rownames(gene_matrix) <- my_gene_ids

# 10. get the grna matrix
grna_sub <- ondisc::get_modality(multimodal_odm_downsample, "grna")
my_grna_ids <- unique(grna_group_data_frame_highmoi$grna_id)
grna_matrix <- grna_sub[[my_grna_ids,]]
rownames(grna_matrix) <- my_grna_ids

# 11. Compute the covariate matrix
covariate_matrix <- multimodal_odm_downsample |>
  ondisc::get_cell_covariates() |>
  dplyr::select(batch = gene_batch) |>
  `rownames<-`(NULL)

# 12. sort according to batch; also remove cells containing no gRNAs
cell_order <- order(covariate_matrix$batch)
response_matrix_highmoi <- gene_matrix[,cell_order]
grna_matrix_highmoi <- grna_matrix[,cell_order]
covariate_data_frame_highmoi <- covariate_matrix[cell_order,,drop = FALSE]

# 13. update gRNA ids
curr_grna_ids <- grna_group_data_frame_highmoi$grna_id
new_grna_ids <- paste0("grna_", substr(curr_grna_ids, 0, 7))
grna_group_data_frame_highmoi$grna_id <- new_grna_ids
rownames(grna_matrix_highmoi) <- new_grna_ids[match(x = rownames(grna_matrix_highmoi), table = curr_grna_ids)]

# 14. rename the data objects
extra_covariates_highmoi <- covariate_data_frame_highmoi |>
  dplyr::mutate(batch = factor(batch, levels = c("prep_batch_1", "prep_batch_2"),
                               labels = c("b1", "b2")))
gene_names_highmoi <- gene_table$gene_name[match(rownames(response_matrix_highmoi), gene_table$gene_id)]
grna_target_data_frame_highmoi <- grna_group_data_frame_highmoi |> dplyr::rename("grna_target" = "grna_group") |> as.data.frame()

# save
highmoi_example_data <- list(response_matrix = response_matrix_highmoi,
                             gene_names = gene_names_highmoi,
                             grna_matrix = grna_matrix_highmoi,
                             extra_covariates = extra_covariates_highmoi)

usethis::use_data(highmoi_example_data, grna_target_data_frame_highmoi, overwrite = TRUE)
