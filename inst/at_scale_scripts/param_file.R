library(ondisc)
param_funct <- function(param) {
        switch(param,
               # modify the parameters below
               storage_dir = "/Users/timbarry/Desktop/example",
               n_processors = 40,
               pod_sizes = c(gene = 3, gRNA = 2, pair = 5),
               expression_matrix = ondisc_matrix(system.file("extdata", "expressions.h5", package = "sceptre")),
               perturbation_matrix = ondisc_matrix(system.file("extdata", "perturbations.h5", package = "sceptre")),
               covariate_matrix = readRDS(system.file("extdata", "covariate_matrix.rds", package = "sceptre")),
               gene_gRNA_pairs = readRDS(system.file("extdata", "gene_gRNA_pairs.rds", package = "sceptre")),
               regularization_amount = 1,
               side = "left",
               seed = 4,
               B = 500
               )}
