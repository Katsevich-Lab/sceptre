library(data.table)

# CellRanger provides a human reference genome, which can be downloaded via the following command:
# curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
# The version of the reference is GRCh38. This script extracts the start position, end position,
# strand, ID, and name of each gene in the reference. The resulting data frame contains 36,572 genes
# and is saved as an internal dataset in the package.

dt <- fread(file = "/Users/timbarry/research_offsite/external/ref/genes.gtf", skip = 5,
            col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
# extract the genes
dt_gene <- dt |> dplyr::filter(feature == "gene")
# keep only those genes on a chromosome
dt_gene_chr <- dt_gene[grep(pattern = "chr", x = dt_gene$chr),]
# extract the gene id and name
attr_split <- strsplit(x = dt_gene_chr$attribute, split = ";", fixed = TRUE)
gene_ids_and_names <- sapply(attr_split, function(elem) {
  gene_id <- strsplit(x = elem[1], split = "\"", fixed = TRUE)[[1]][2]
  gene_name <- strsplit(x = elem[4], split = "\"", fixed = TRUE)[[1]][2]
  c(gene_id = gene_id, gene_name = gene_name)
}) |> t() |> as.data.frame()
# append these columns to dt_gene
gene_table <- cbind(dt_gene_chr[,c("chr", "start", "end", "strand")], gene_ids_and_names) |>
  dplyr::mutate(chr = factor(chr)) |> dplyr::mutate(tss_position = ifelse(strand == "+", start, end)) |>
  dplyr::select(-start, -end, -strand)
data.table::setorderv(gene_table, c("chr", "tss_position"))
gene_table <- gene_table

# save internally
usethis::use_data(gene_table, internal = TRUE, overwrite = TRUE)
