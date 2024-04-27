library(data.table)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)

#############
# hg 38 table
#############
# CellRanger provides a human reference genome, which can be downloaded via the following command:
# curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
# The version of the reference is GRCh38. This script extracts the start position, end position,
# strand, ID, and name of each gene in the reference. The resulting data frame contains 36,572 genes
# and is saved as an internal dataset in the package called gene_position_data_frame_grch38. The
# file `genes.gtf` that we load is contained within the download.

dt <- fread(file = "~/research_offsite/external/ref/genes.gtf", skip = 5,
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
  c(response_id = gene_id, response_name = gene_name)
}) |> t() |> as.data.frame()
# append these columns to dt_gene
gene_table <- cbind(dt_gene_chr[,c("chr", "start", "end", "strand")], gene_ids_and_names) |>
  dplyr::mutate(chr = factor(chr)) |> dplyr::mutate(position = ifelse(strand == "+", start, end)) |>
  dplyr::select(-start, -end, -strand)
data.table::setorderv(gene_table, c("chr", "position"))
gene_position_data_frame_grch38 <- gene_table
usethis::use_data(gene_position_data_frame_grch38, internal = FALSE, overwrite = TRUE)

#############
# hg 19 table
#############
rm(list = ls())
# We obtained the hg37 reference genome from cellranger
# wget ftp://ftp.ensembl.org/pub/grch37/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh37.82.gtf.gz
library(rtracklayer)
dt <- readGFF("~/research_offsite/external/ref/Homo_sapiens.GRCh37.82.gtf.gz")
# retain only genes
dt <- dt |> dplyr::filter(type == "gene")
# keep only those genes on a chromosome
dt <- dt[!grepl(pattern = "^G", x = dt$seqid),]
dt$chr <- paste0("chr", as.character(dt$seqid))
dt <- dt |>
  dplyr::mutate(position = ifelse(strand == "+", start, end)) |>
  dplyr::select(response_id = gene_id,
                response_name = gene_name,
                chr, position)
gene_position_data_frame_grch19 <- dt
gene_position_data_frame_grch19$chr <- factor(gene_position_data_frame_grch19$chr)
usethis::use_data(gene_position_data_frame_grch19, internal = FALSE, overwrite = TRUE)
