data(response_matrix_highmoi)
data(grna_matrix_highmoi)
data(extra_covariates_highmoi)
data(grna_target_data_frame_highmoi)
data(gene_names_highmoi)

sceptre_object <- import_data(
  response_matrix = response_matrix_highmoi,
  grna_matrix = grna_matrix_highmoi,
  grna_target_data_frame = grna_target_data_frame_highmoi,
  moi = "high",
  extra_covariates = extra_covariates_highmoi,
  response_names = gene_names_highmoi)

write_sceptre_object_to_cellranger_format(sceptre_object = sceptre_object,
                                          directory = "~/research_code/sceptre/inst/extdata/highmoi_example")
