#!/bin/bash
param_file_name="param_file.R"
param_file=$PWD"/"$param_file_name

# 0. Print parameter file being used
printf "Running SCEPTRE on parameter file $param_file.\n"

# 1. Obtain some basic parameters from param file
n_proc=$(Rscript -e "source('$param_file');\
cat(param_funct('n_processors'));")

# 2. create the offsite directory structure
printf "Creating storage directory structure.\n"
read -r gene_dir gRNA_dir res_dir log_dir <<<$(Rscript -e \
"library(sceptre); source('$param_file');\
storage_dir <- param_funct('storage_dir');\
dir_locs <- initialize_directories(storage_dir);\
cat(dir_locs);")

# 3. Create the file dictionaries
printf "Creating dictionaries.\n"
read -r n_gene_pods n_gRNA_pods n_pair_pods <<<$(Rscript -e \
"library(sceptre); source('$param_file');\
gene_gRNA_pairs <- param_funct('gene_gRNA_pairs');\
pod_sizes <- param_funct('pod_sizes');\
dicts <- create_and_store_dictionaries(gene_gRNA_pairs, '$gene_dir', '$gRNA_dir', '$res_dir', pod_sizes);\
cat(dicts)")

# 4. Run first round of gene precomputations
printf "Run first round of gene precomputations.\n"
for i in $(seq 1 $n_gene_pods)
  do
    Rscript -e \
    "library(sceptre); source('$param_file');\
    expression_matrix <- param_funct('expression_matrix');\
    covariate_matrix <- param_funct('covariate_matrix');\
    regularization_amount <- param_funct('regularization_amount');\
    run_gene_precomputation_at_scale_round_1($i, '$gene_dir', expression_matrix, covariate_matrix, regularization_amount, '$log_dir');" &
  done
wait

# 5. Regularize thetas
printf "Regularizing thetas.\n"
Rscript -e "library(sceptre); source('$param_file');\
regularization_amount <- param_funct('regularization_amount');\
regularize_gene_sizes_at_scale('$gene_dir', regularization_amount, '$log_dir');"

# 6. Run second round of gene precomputations
printf "Running second round of gene precomputations.\n"
for i in $(seq 1 $n_gene_pods)
  do
    Rscript -e \
    "library(sceptre); source('$param_file');\
    expression_matrix <- param_funct('expression_matrix');\
    covariate_matrix <- param_funct('covariate_matrix');\
    regularization_amount <- param_funct('regularization_amount');\
    run_gene_precomputation_at_scale_round_2($i, '$gene_dir', expression_matrix, covariate_matrix, regularization_amount, '$log_dir')" &
  done

# 7. Run gRNA precomputations.
printf "Running gRNA precomputations.\n"
for i in $(seq 1 $n_gRNA_pods)
  do
    Rscript -e \
    "library(sceptre); source('$param_file');\
    perturbation_matrix <- param_funct('perturbation_matrix');\
    covariate_matrix <- param_funct('covariate_matrix');\
    run_gRNA_precomputation_at_scale($i, '$gRNA_dir', perturbation_matrix, covariate_matrix, '$log_dir');" &
  done
wait

# 8. Run gene-gRNA pair analysis
printf "Computing p-values for gene-gRNA pairs.\n"
for i in $(seq 1 $n_pair_pods)
  do
    Rscript -e \
    "
    library(sceptre); source('$param_file');\
    expression_matrix <- param_funct('expression_matrix');\
    perturbation_matrix <- param_funct('perturbation_matrix');\
    covariate_matrix <- param_funct('covariate_matrix');\
    regularization_amount <- param_funct('regularization_amount');\
    side <- param_funct('side');\
    seed <- param_funct('seed');\
    B <- param_funct('B');\
    run_gRNA_gene_pair_analysis_at_scale($i, '$gene_dir', '$gRNA_dir', '$res_dir', '$log_dir', expression_matrix, perturbation_matrix, covariate_matrix, regularization_amount, side, seed, B)
    " &
  done
wait

#9 . Combine results
printf "Combining results.\n"
Rscript -e "library(sceptre); collect_results('$res_dir');"
