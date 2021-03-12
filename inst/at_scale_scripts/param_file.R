library(ondisc)
param_funct <- function(param) {
        switch(param,
               # modify the parameters below
               storage_dir = "~/my_sceptre_dir/",
               expression_matrix = ondisc_matrix(system.file("extdata", "expressions.h5", package = "sceptre")),
               perturbation_matrix = ondisc_matrix(system.file("extdata", "perturbations.h5", package = "sceptre")),
               covariate_matrix = readRDS(system.file("extdata", "covariate_matrix.rds", package = "sceptre")),
               gene_gRNA_pairs = readRDS(system.file("extdata", "gene_gRNA_pairs.rds", package = "sceptre")),
               side = "left",
               pod_sizes = c(gene = 3, gRNA = 2, pair = 5),
               regularization_amount = 1,
               seed = 4,
               B = 500)
}
