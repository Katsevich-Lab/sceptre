###############################################################################
# HIGH-MOI CIS-REGULATORY PERTURB-SEQ SIMULATION
# - Gene-targeting gRNAs use *real ENSG IDs* in grna_target_data_frame
# - Enhancers remain "candidate_enh_...", but store local ENSG IDs
###############################################################################

# 0. LOAD DATA AND LIBRARIES --------------------------------------------------
# This assumes you have 'gene_position_data_frame_grch38' in your environment.
# e.g., data(gene_position_data_frame_grch38)
library(Matrix)
library(dplyr)
library(sceptre)
data(gene_position_data_frame_grch38)

###############################################################################
# 1. SIMULATION PARAMETERS
###############################################################################

RANDOM_SEED <- 1
set.seed(RANDOM_SEED)

# Genes measured (rows) in expression matrix
NUM_GENES <- 100

# Of these genes, how many are *directly targeted* with gene gRNAs?
NUM_GENE_TARGETS <- 5

# Number of candidate enhancers across the gene targets
NUM_ENHANCER_TARGETS <- 20

# Number of genes with nearby enhancers
NUM_ENHANCER_GENES <- 10

# Constraints on genes with nearby enhancers
ENHANCER_GENE_CHR <- "chr22"
ENHANCER_GENE_MIN_POS <- 2e7
ENHANCER_GENE_MAX_POS <- 3e7

# Fraction of enhancers that truly regulate their local gene
FRACTION_OF_ENHANCERS_REGULATING <- 1

# gRNA design details
GRNAS_PER_GENE <- 2
GRNAS_PER_ENH  <- 2
NUM_NT_GRNAS   <- 10

# Windows around gene TSS for gene-based gRNAs
GENE_WINDOW_SHIFT <- 50
GENE_WINDOW_SIZE  <- 50

# Enhancer location parameters
ENHANCER_CIS_WINDOW <- 5e5
ENHANCER_LENGTH     <- 200
ENHANCER_GRNA_SHIFT <- 50

# Average MOI => expected # of gRNAs per cell
AVG_MOI <- 10

# Number of cells
NUM_CELLS <- 500

# Gene expression parameters
MEAN_LOG_GENE_EXPRESSION <- log(1)
SD_LOG_GENE_EXPRESSION   <- 1.25
GENE_DISPERSION          <- 1

# gRNA expression parameters
MEAN_LOG_BACKGROUND_GRNA_EXPRESSION <- log(0.05)
SD_LOG_BACKGROUND_GRNA_EXPRESSION   <- 0.1
MEAN_LOG_GRNA_EXPRESSION           <- log(20)
SD_LOG_GRNA_EXPRESSION             <- 0.5

# On-target effect (log-fold change)
ON_TARGET_LOG_FOLD_CHANGE <- log(0.25)

# Enhancer effects on genes (log-fold change)
ENHANCER_GENE_LOG_FOLD_CHANGE <- log(0.5)

# Batch effects
NUM_BATCHES           <- 2
BATCH_LOG_FOLD_CHANGE <- log(0.8)

###############################################################################
# 2. SELECT GENES FOR EXPRESSION + FOR DIRECT GENE-TARGETING + FOR ENHANCER-TARGETING
###############################################################################
# 'gene_position_data_frame_grch38' has columns (response_id, response_name, chr, position, ...)
# We'll rename them for clarity, but our final row names in 'gene_expression_matrix' should be
# the real ENSG IDs: e.g., "ENSG00000123456".

selected_expression_genes <- gene_position_data_frame_grch38 %>%
  filter(chr == ENHANCER_GENE_CHR,
         position >= ENHANCER_GENE_MIN_POS,
         position <= ENHANCER_GENE_MAX_POS) %>%
  rename(
    gene_id   = response_id,   # e.g. "ENSG00000..."
    gene_name = response_name,
    gene_chr  = chr,
    gene_tss  = position
  ) %>%
  # randomly pick 100 to simulate
  sample_n(NUM_GENES)

# Out of those genes, pick subset for direct perturbation
selected_perturbed_genes <- selected_expression_genes %>%
  sample_n(NUM_GENE_TARGETS)
# Each row has: gene_id, gene_name, gene_chr, gene_tss

# We'll store these ENSG IDs for convenience
perturbed_gene_ids <- selected_perturbed_genes$gene_id

# Pick genes with enhancers
selected_enhancer_genes <- selected_expression_genes %>%
  sample_n(NUM_ENHANCER_GENES)

# Store these ENSG IDs for convenience
enhancer_gene_ids <- selected_enhancer_genes$gene_id

###############################################################################
# 3. GENE-TARGETING GNRAS (USING REAL ENSG IDS)
###############################################################################
# For each of the 3 perturbed genes, place GRNAS_PER_GENE around TSS.
# The key: grna_target == the real 'gene_id'.

gene_target_list <- list()

for (i in seq_len(nrow(selected_perturbed_genes))) {
  gene_id_i  <- selected_perturbed_genes$gene_id[i]   # e.g. "ENSG0000012345"
  gene_chr_i <- selected_perturbed_genes$gene_chr[i]
  gene_tss_i <- selected_perturbed_genes$gene_tss[i]

  starts_i <- gene_tss_i + seq.int(0, by = GENE_WINDOW_SHIFT, length.out = GRNAS_PER_GENE)
  ends_i   <- starts_i + (GENE_WINDOW_SIZE - 1)

  df_i <- data.frame(
    grna_id     = paste0(gene_id_i, "_grna", seq_len(GRNAS_PER_GENE)),
    grna_target = gene_id_i,  # REAL ENSG ID
    chr         = gene_chr_i,
    start       = starts_i,
    end         = ends_i,
    stringsAsFactors = FALSE
  )
  gene_target_list[[i]] <- df_i
}

gene_target_data_frame <- do.call(rbind, gene_target_list)

###############################################################################
# 4. ENHANCER-TARGETING GNRAS (CIS, ±1 MB), "candidate_enh_..."
###############################################################################
enhancers_per_gene <- rep(
  floor(NUM_ENHANCER_TARGETS / NUM_ENHANCER_GENES),
  NUM_ENHANCER_GENES
)
remainder <- NUM_ENHANCER_TARGETS - sum(enhancers_per_gene)
if (remainder > 0) {
  enhancers_per_gene[length(enhancers_per_gene)] <-
    enhancers_per_gene[length(enhancers_per_gene)] + remainder
}

enh_rows_list <- list()
enh_index <- 1

for (i in seq_len(nrow(selected_enhancer_genes))) {
  gene_id_i  <- selected_enhancer_genes$gene_id[i]
  gene_chr_i <- selected_enhancer_genes$gene_chr[i]
  gene_tss_i <- selected_enhancer_genes$gene_tss[i]

  num_enh_i <- enhancers_per_gene[i]
  if (num_enh_i == 0) next

  # fraction that truly regulate
  num_regulating <- floor(FRACTION_OF_ENHANCERS_REGULATING * num_enh_i)
  is_regulating_vec <- c(
    rep(TRUE,  num_regulating),
    rep(FALSE, num_enh_i - num_regulating)
  )
  is_regulating_vec <- sample(is_regulating_vec)

  # random centers ±ENHANCER_CIS_WINDOW
  region_min <- gene_tss_i - ENHANCER_CIS_WINDOW
  region_max <- gene_tss_i + ENHANCER_CIS_WINDOW

  enh_centers <- sample(seq.int(region_min, region_max),
                        size = num_enh_i, replace = TRUE)

  enh_starts <- enh_centers
  enh_ends   <- enh_centers + (ENHANCER_LENGTH - 1)

  # build enhancer labels: "candidate_enh_1", "candidate_enh_2", ...
  candidate_enh_labels <- paste0("candidate_enh_", enh_index + seq_len(num_enh_i) - 1)

  # each enhancer has GRNAS_PER_ENH gRNAs
  all_enh_df <- list()
  for (j in seq_len(num_enh_i)) {
    starts_j <- enh_starts[j] + seq.int(0, by = ENHANCER_GRNA_SHIFT, length.out = GRNAS_PER_ENH)
    ends_j   <- starts_j + (ENHANCER_GRNA_SHIFT - 1)

    df_j <- data.frame(
      grna_id     = paste0(candidate_enh_labels[j], "_grna", seq_len(GRNAS_PER_ENH)),
      grna_target = candidate_enh_labels[j],  # "candidate_enh_..."
      chr         = gene_chr_i,
      start       = starts_j,
      end         = ends_j,
      is_regulating = is_regulating_vec[j],
      # store the real gene ID that this enhancer is near
      local_gene    = gene_id_i,
      stringsAsFactors = FALSE
    )
    all_enh_df[[j]] <- df_j
  }

  enh_rows_list[[i]] <- do.call(rbind, all_enh_df)
  enh_index <- enh_index + num_enh_i
}

enhancer_target_data_frame <- do.call(rbind, enh_rows_list)

###############################################################################
# 5. NON-TARGETING GNRAS (NA POSITIONS)
###############################################################################
nt_rows <- data.frame(
  grna_id     = paste0("non-targeting_grna", seq_len(NUM_NT_GRNAS)),
  grna_target = "non-targeting",
  chr         = NA_character_,
  start       = NA,
  end         = NA,
  stringsAsFactors = FALSE
)

###############################################################################
# 6. COMBINE GENE, ENHANCER, AND NON-TARGETING GNRAS
###############################################################################
# STEP A: Keep a *full copy* of the enhancer table so we can build the map without merges:
enhancer_target_data_frame_full <- enhancer_target_data_frame

# STEP B: Build 'enhancer_regulatory_map' directly:
enhancer_regulatory_map <- enhancer_target_data_frame_full %>%
  filter(is_regulating) %>%              # i.e., TRUE
  distinct(grna_target, local_gene)

# STEP C: Now remove 'is_regulating'/'local_gene' only when binding:
enhancer_target_data_frame <- enhancer_target_data_frame_full %>%
  select(-is_regulating, -local_gene)

# STEP D: Combine gene-targeting gRNAs, enhancer-targeting gRNAs, and non-targeting gRNAs
grna_target_data_frame <- dplyr::bind_rows(
  gene_target_data_frame,
  enhancer_target_data_frame,
  nt_rows
) %>%
  # reorder columns if desired
  select(grna_id, grna_target, chr, start, end)

NUM_GRNAS <- nrow(grna_target_data_frame)

###############################################################################
# 7. ASSIGN MULTIPLE GNRAS PER CELL (HIGH MOI)
###############################################################################
gRNAs_per_cell_count <- rpois(n = NUM_CELLS, lambda = AVG_MOI)
assigned_gRNAs_list <- vector("list", length = NUM_CELLS)

for (cell_idx in seq_len(NUM_CELLS)) {
  k <- gRNAs_per_cell_count[cell_idx]
  if (k == 0) {
    assigned_gRNAs_list[[cell_idx]] <- integer(0)
  } else {
    assigned_gRNAs_list[[cell_idx]] <- sample.int(NUM_GRNAS, size = k, replace = FALSE)
  }
}

###############################################################################
# 8. GENERATE GRNA EXPRESSION MATRIX
###############################################################################
grna_matrix <- matrix(
  rpois(
    n = NUM_GRNAS * NUM_CELLS,
    lambda = exp(rnorm(
      n = NUM_GRNAS * NUM_CELLS,
      mean = MEAN_LOG_BACKGROUND_GRNA_EXPRESSION,
      sd   = SD_LOG_BACKGROUND_GRNA_EXPRESSION
    ))
  ),
  nrow = NUM_GRNAS,
  ncol = NUM_CELLS,
  dimnames = list(grna_target_data_frame$grna_id, paste0("Cell_", seq_len(NUM_CELLS)))
)

# override expression for gRNAs actually present in a cell
for (cell_idx in seq_len(NUM_CELLS)) {
  g_in_cell <- assigned_gRNAs_list[[cell_idx]]
  if (length(g_in_cell) > 0) {
    lam_vec <- exp(rnorm(
      n    = length(g_in_cell),
      mean = MEAN_LOG_GRNA_EXPRESSION,
      sd   = SD_LOG_GRNA_EXPRESSION
    ))
    grna_matrix[g_in_cell, cell_idx] <- rpois(length(g_in_cell), lam_vec)
  }
}

###############################################################################
# 9. BUILD LOG-FOLD CHANGE MATRIX (ROWS=GNRAS, COLS=GENES)
###############################################################################
# We'll label columns by the real ENSG IDs from 'selected_expression_genes$gene_id'.
all_gene_ids <- selected_expression_genes$gene_id  # length=NUM_GENES

grna_log_fold_change_matrix <- matrix(
  0,
  nrow = NUM_GRNAS,
  ncol = NUM_GENES,
  dimnames = list(grna_target_data_frame$grna_id, all_gene_ids)
)

# (A) On-target for gene-targeting gRNAs: each row has grna_target == real ENSG ID
for (ensg_id in perturbed_gene_ids) {
  row_idxs <- which(grna_target_data_frame$grna_target == ensg_id)
  col_j    <- match(ensg_id, all_gene_ids)

  if (!is.na(col_j) && length(row_idxs) > 0) {
    grna_log_fold_change_matrix[row_idxs, col_j] <- ON_TARGET_LOG_FOLD_CHANGE
  }
}

# (B) On-target for enhancers that truly regulate a local gene (ENSG)
# 'enhancer_regulatory_map' has (grna_target="candidate_enh_...", local_gene=<ENSG>).
for (k in seq_len(nrow(enhancer_regulatory_map))) {
  enh_label  <- enhancer_regulatory_map$grna_target[k] # "candidate_enh_..."
  local_ensg <- enhancer_regulatory_map$local_gene[k]  # e.g. "ENSG0000012345"

  row_idxs <- which(grna_target_data_frame$grna_target == enh_label)
  col_j    <- match(local_ensg, all_gene_ids)

  if (!is.na(col_j) && length(row_idxs) > 0) {
    grna_log_fold_change_matrix[row_idxs, col_j] <- ENHANCER_GENE_LOG_FOLD_CHANGE
  }
}

###############################################################################
# 10. BATCH COVARIATES
###############################################################################
batch_indicator <- sample(seq_len(NUM_BATCHES), size = NUM_CELLS, replace = TRUE)
extra_covariates <- data.frame(batch = factor(paste0("batch_", batch_indicator)))

batch_effect <- ifelse(batch_indicator == 1, 0, BATCH_LOG_FOLD_CHANGE)

###############################################################################
# 11. SIMULATE GENE EXPRESSION
###############################################################################
# (A) Baseline log-normal means for each gene
mu_gene <- exp(rnorm(
  NUM_GENES,
  mean = MEAN_LOG_GENE_EXPRESSION,
  sd   = SD_LOG_GENE_EXPRESSION
))
names(mu_gene) <- selected_expression_genes$gene_id
# set expression of targeted genes to a higher value
mu_gene[perturbed_gene_ids] <- exp(MEAN_LOG_GENE_EXPRESSION + SD_LOG_GENE_EXPRESSION)

# (B) Sum up all gRNA effects per cell => (genes x cells)
grna_effect_matrix <- matrix(0, nrow = NUM_GENES, ncol = NUM_CELLS)

for (cell_idx in seq_len(NUM_CELLS)) {
  gRNA_rows_in_cell <- assigned_gRNAs_list[[cell_idx]]
  if (length(gRNA_rows_in_cell) > 0) {
    grna_effect_matrix[, cell_idx] <-
      colSums(grna_log_fold_change_matrix[gRNA_rows_in_cell, , drop = FALSE])
  }
}

# (C) Add batch effect
batch_effect_matrix <- matrix(batch_effect,
                              nrow = NUM_GENES,
                              ncol = NUM_CELLS, byrow = TRUE)

# (D) Final log-scale means
log_mu_matrix <- matrix(log(mu_gene),
                        nrow = NUM_GENES, ncol = NUM_CELLS, byrow = FALSE
) + grna_effect_matrix + batch_effect_matrix

# (E) Convert to linear scale
mu_expression <- exp(log_mu_matrix)

# (F) Negative binomial sampling
size_parameter <- 1 / GENE_DISPERSION

gene_expression_vector <- rnbinom(
  n   = NUM_GENES * NUM_CELLS,
  mu  = mu_expression,
  size = size_parameter
)

# (G) Build final expression matrix; rownames = real ENSG IDs
gene_expression_matrix <- matrix(
  gene_expression_vector,
  nrow = NUM_GENES,
  ncol = NUM_CELLS,
  dimnames = list(all_gene_ids, paste0("Cell_", seq_len(NUM_CELLS)))
)

###############################################################################
# 12. PACKAGE AND SAVE
###############################################################################

# get the gene names
all_gene_names <- gene_position_data_frame_grch38 |> pull(response_name)
names(all_gene_names) <- gene_position_data_frame_grch38 |> pull(response_id)
gene_names <- all_gene_names[selected_expression_genes$gene_id] |> unname()

# package and save high-MOI example data
highmoi_example_data <- list(
  response_matrix = gene_expression_matrix,
  grna_matrix = grna_matrix,
  extra_covariates = extra_covariates,
  gene_names = gene_names
)
usethis::use_data(highmoi_example_data, overwrite = TRUE)

# save grna_target_data_frame
grna_target_data_frame_highmoi <- grna_target_data_frame
usethis::use_data(grna_target_data_frame_highmoi, overwrite = TRUE)

# write to 10x Cell Ranger format
ondisc:::write_sceptre_object_to_cellranger_format_v2(
  mats = list(gene = Matrix(gene_expression_matrix, sparse = TRUE),
              grna = Matrix(grna_matrix, sparse = TRUE)),
  gene_names <- gene_names,
  directory = file.path("inst", "extdata", "highmoi_example"),
  batch = paste0("gem_group_", batch_indicator)
)

###############################################################################
# 13. TEST
###############################################################################

# sceptre_object <- import_data(
#   response_matrix = gene_expression_matrix,
#   grna_matrix = grna_matrix,
#   extra_covariates = extra_covariates,
#   response_names = gene_names,
#   grna_target_data_frame = grna_target_data_frame,
#   moi = "high"
# )
#
# positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
# discovery_pairs <- construct_cis_pairs(
#   sceptre_object = sceptre_object,
#   positive_control_pairs = positive_control_pairs
# )
#
# sceptre_object <- set_analysis_parameters(
#   sceptre_object = sceptre_object,
#   discovery_pairs = discovery_pairs,
#   positive_control_pairs = positive_control_pairs,
#   side = "left"
# )
# print(sceptre_object)
#
# # 3. assign grnas
# plot_grna_count_distributions(sceptre_object)
# sceptre_object <- sceptre_object |> assign_grnas()
# plot(sceptre_object)
# print(sceptre_object)
#
# # 4. run qc
# plot_covariates(sceptre_object, p_mito_threshold = 0.075)
# sceptre_object <- sceptre_object |> run_qc(p_mito_threshold = 0.075)
# plot(sceptre_object)
# print(sceptre_object)
#
# # 5. run the calibration check
# sceptre_object <- run_calibration_check(sceptre_object, parallel = TRUE, n_processors = 2)
# plot(sceptre_object)
# print(sceptre_object)
#
# # 6. run power check
# sceptre_object <- run_power_check(sceptre_object, parallel = TRUE, n_processors = 2)
# plot(sceptre_object)
# print(sceptre_object)
#
# # 7. run discovery analysis
# sceptre_object <- run_discovery_analysis(sceptre_object, parallel = TRUE, n_processors = 2)
# plot(sceptre_object)
# print(sceptre_object)
