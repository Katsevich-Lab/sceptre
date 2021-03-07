args <- commandArgs(trailingOnly = TRUE)
code_dir <- if (is.na(args[1])) "/Users/timbarry/Box/SCEPTRE-manuscript/SCEPTRE/" else args[1]
source(paste0(code_dir, "/sceptre_paper/analysis_drivers/analysis_drivers_gasp/file_paths_to_dirs.R"))
source(paste0(code_dir, "sceptre_paper/analysis_drivers/analysis_drivers_gasp/sceptre_function_args.R"))

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
