library(sceptre)
data(gene_position_data_frame_grch38)

# dataset parameters -----------------------------------------------------------
num_cells <- 1000
num_genes <- 100
num_targets <- 20
grnas_per_target <- 2
num_nt_grnas <- 10
mean_log_gene_expression <- log(1)
sd_log_gene_expression <- 0.5
gene_dispersion <- 1
mean_log_background_grna_expression <- log(0.05)
sd_log_background_grna_expression <- 0.1
mean_log_grna_expression <- log(20)
sd_log_grna_expression <- 0.5
nonnull_fraction <- 0.1
mean_log_fold_change <- 0 # to have positive and negative effects
sd_log_fold_change <- 1
on_target_log_fold_change <- log(0.25)
num_batches <- 2
batch_log_fold_change <- log(1.2)

# generate the data -----------------------------------------------------------

# Set seed for reproducibility
set.seed(123)

# generate the gRNA target data frame
gene_names <- sample(gene_position_data_frame_grch38$response_id, num_genes)
grna_target_data_frame <- dplyr::tibble(
  grna_target = c(
    rep(sample(gene_names, num_targets), each = grnas_per_target),
    rep("nt", num_nt_grnas)
  ),
  grna_id = paste0(grna_target, "_grna", c(
    rep(1:grnas_per_target, num_targets),
    1:num_nt_grnas
  ))
) |>
  dplyr::mutate(grna_target = ifelse(grna_target == "nt",
                                     "non-targeting",
                                     grna_target)) |>
  dplyr::relocate(grna_id)

# generate which cells received which gRNAs via multinomial sampling
num_grnas <- nrow(grna_target_data_frame)
grna_per_cell <- rmultinom(
  n = num_cells,
  size = 1,
  prob = rep(1 / num_grnas, num_grnas)
) |>
  apply(MARGIN = 2, FUN = function(col) which(col == 1))

# generate gRNA expression matrix based on the vector grna_per_cell; gRNAs not in a
# cell get a background expression level
grna_matrix <- matrix(
  rpois(
    n = num_grnas * num_cells,
    lambda = exp(rnorm(num_grnas * num_cells,
      mean = mean_log_background_grna_expression,
      sd = sd_log_background_grna_expression
    ))
  ),
  nrow = num_grnas,
  ncol = num_cells,
  dimnames = list(grna_target_data_frame$grna_id, paste0("Cell_", 1:num_cells))
)
grna_specific_values <- rpois(
  n = num_cells,
  lambda = exp(rnorm(
    n = num_cells,
    mean = mean_log_grna_expression,
    sd = sd_log_grna_expression
  ))
)
grna_matrix[cbind(grna_per_cell, 1:num_cells)] <- grna_specific_values


# matrix of gRNA effects on genes

# Initialize the log-fold change matrix with zeros
grna_log_fold_change_matrix <- matrix(
  0,
  nrow = num_grnas,
  ncol = num_genes,
  dimnames = list(grna_target_data_frame$grna_id, gene_names)
)

# Assign On-Target Log-Fold Changes
# Vectorized approach for efficiency
# Create a vector of column indices for target genes
target_gene_indices <- match(grna_target_data_frame$grna_target, gene_names)

# Assign the on-target log-fold change
grna_log_fold_change_matrix[cbind(1:num_grnas, target_gene_indices)] <- on_target_log_fold_change

# Identify Off-Target Pairs
# Create a logical matrix indicating on-target pairs
on_target_matrix <- matrix(FALSE, nrow = num_grnas, ncol = num_genes)
on_target_matrix[cbind(1:num_grnas, target_gene_indices)] <- TRUE

# All possible pairs are rows (gRNAs) x columns (genes)
# Off-target pairs are where on_target_matrix is FALSE
off_target_indices <- which(!on_target_matrix, arr.ind = TRUE)
# Exclude NTs from off-target indices
nt_grna_rows <- which(grna_target_data_frame$grna_target == "non-targeting")
off_target_indices <- off_target_indices[!off_target_indices[, 1] %in% nt_grna_rows, ]

# Total number of off-target pairs
num_off_target <- nrow(off_target_indices)

# Number of non-null off-target effects
num_nonnull <- floor(nonnull_fraction * num_off_target)

# Randomly select indices for non-null off-target effects
selected_off_target_indices <- off_target_indices[sample(1:num_off_target, num_nonnull), ]

# Assign random log-fold changes to selected off-target pairs
grna_log_fold_change_matrix[cbind(selected_off_target_indices[, 1], selected_off_target_indices[, 2])] <-
  rnorm(num_nonnull, mean = mean_log_fold_change, sd = sd_log_fold_change)

# vector of batch indicators
batch_indicator <- sample(1:num_batches, size = num_cells, replace = TRUE)
extra_covariates <- data.frame(batch = factor(paste0("batch_", batch_indicator)))

# matrix of gene expressions

# ----- Define Batch Effects -----
# Assign log-fold change based on batch
# Batch 1: 0 (no effect), Batch 2: log(1.2)
batch_effect <- ifelse(batch_indicator == 1, 0, batch_log_fold_change)

# ----- Generate Baseline Gene Expression -----
# For each gene, generate a baseline mean expression from a log-normal distribution
mu_gene <- exp(rnorm(num_genes, mean = mean_log_gene_expression, sd = sd_log_gene_expression))

# ----- Incorporate gRNA Effects -----
# grna_log_fold_change_matrix: num_grnas x num_genes
# grna_per_cell: vector of length num_cells, values from 1 to num_grnas

# Extract gRNA-specific log-fold changes for each cell
# This results in a num_cells x num_genes matrix
grna_effect_matrix <- grna_log_fold_change_matrix[grna_per_cell, ]

# Transpose to get a num_genes x num_cells matrix
grna_effect_matrix <- t(grna_effect_matrix) # Now rows = genes, columns = cells

# ----- Incorporate Batch Effects -----
# Create a num_genes x num_cells matrix where each column has the batch effect added to all genes
batch_effect_matrix <- matrix(batch_effect, nrow = num_genes, ncol = num_cells, byrow = TRUE)

# ----- Calculate Mean Expression Matrix -----
# For each gene and cell, compute the mean expression incorporating baseline, gRNA, and batch effects
# log_mu_matrix = log(mu_gene) + grna_effect + batch_effect
log_mu_matrix <- matrix(log(mu_gene), nrow = num_genes, ncol = num_cells, byrow = FALSE) +
  grna_effect_matrix +
  batch_effect_matrix

# Convert log-scale to linear scale
mu_expression <- exp(log_mu_matrix)

# ----- Simulate Gene Expression Counts -----
# Negative Binomial parameters:
# - mu = mu_expression
# - size = 1 / gene_dispersion
# Note: In R's rnbinom, 'size' is the dispersion parameter where variance = mu + mu^2 / size
size_parameter <- 1 / gene_dispersion # Here, size_parameter = 1

# Simulate counts for all gene-cell pairs
gene_expression_vector <- rnbinom(n = num_genes * num_cells,
                                  mu = mu_expression,
                                  size = size_parameter)

# Reshape the vector into a matrix
gene_expression_matrix <- matrix(
  gene_expression_vector,
  nrow = num_genes,
  ncol = num_cells,
  dimnames = list(gene_names, paste0("Cell_", 1:num_cells))
)

# ----- Package and save -----

lowmoi_example_data <- list(
  response_matrix = gene_expression_matrix,
  grna_matrix = grna_matrix,
  grna_target_data_frame = grna_target_data_frame,
  extra_covariates = extra_covariates
)
usethis::use_data(lowmoi_example_data, overwrite = TRUE)

# ----- Test --------

# library(sceptre)
#
# sceptre_object <- import_data(
#   response_matrix = gene_expression_matrix,
#   grna_matrix = grna_matrix,
#   extra_covariates = extra_covariates,
#   response_names = gene_names,
#   grna_target_data_frame = grna_target_data_frame,
#   moi = "low"
# )
#
# positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
# discovery_pairs <- construct_trans_pairs(
#   sceptre_object = sceptre_object,
#   positive_control_pairs = positive_control_pairs,
#   pairs_to_exclude = "pc_pairs"
# )
#
# sceptre_object <- set_analysis_parameters(
#   sceptre_object = sceptre_object,
#   discovery_pairs = discovery_pairs,
#   positive_control_pairs = positive_control_pairs
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
