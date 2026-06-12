# Builds the gene-position reference tables shipped in data/. Each section
# downloads a public reference genome, extracts the per-gene coordinates,
# and saves the result with usethis::use_data().

library(data.table)

###############
# grch 38 table
###############
# The genes/genes.gtf file is contained within the 10x Cell Ranger human
# reference (GRCh38, 2020-A), which 10x Cell Ranger has used since 2020.
# This script downloads that reference (~11 GB), extracts genes.gtf, and
# pulls out the start position, end position, strand, ID, and name of each
# gene. The resulting data frame contains 36,572 genes and is saved as the
# package dataset gene_position_data_frame_grch38.
url <- paste0("https://cf.10xgenomics.com/supp/cell-exp/",
              "refdata-gex-GRCh38-2020-A.tar.gz")
tarball <- file.path(tempdir(), "refdata-gex-GRCh38-2020-A.tar.gz")
download.file(url, tarball, mode = "wb")
member <- "refdata-gex-GRCh38-2020-A/genes/genes.gtf"
untar(tarball, files = member, exdir = tempdir())
gtf <- file.path(tempdir(), member)

dt <- fread(file = gtf, skip = 5,
            col.names = c("chr", "source", "feature", "start", "end",
                          "score", "strand", "frame", "attribute"))

# extract the genes
dt_gene <- dt |> dplyr::filter(feature == "gene")
# keep only those genes on a chromosome
dt_gene_chr <- dt_gene[grep(pattern = "chr", x = dt_gene$chr), ]
# extract the gene id and name
attr_split <- strsplit(x = dt_gene_chr$attribute, split = ";", fixed = TRUE)
gene_ids_and_names <- sapply(attr_split, function(elem) {
    gene_id <- strsplit(x = elem[1], split = "\"", fixed = TRUE)[[1]][2]
    gene_name <- strsplit(x = elem[4], split = "\"", fixed = TRUE)[[1]][2]
    c(response_id = gene_id, response_name = gene_name)
}) |> t() |> as.data.frame()
# append these columns to dt_gene
gene_table <- cbind(dt_gene_chr[, c("chr", "start", "end", "strand")],
                    gene_ids_and_names) |>
    dplyr::mutate(chr = factor(chr)) |>
    dplyr::mutate(position = ifelse(strand == "+", start, end)) |>
    dplyr::select(-start, -end, -strand)
data.table::setorderv(gene_table, c("chr", "position"))
gene_position_data_frame_grch38 <- gene_table |>
    dplyr::select(response_id, response_name, chr, position)
usethis::use_data(gene_position_data_frame_grch38, internal = FALSE,
                  overwrite = TRUE)

###############
# grch 37 table
###############
rm(list = ls())
library(data.table)
library(rtracklayer)
# The GRCh37 reference genome is obtained from Ensembl (release 84).
url <- paste0("https://ftp.ensembl.org/pub/grch37/release-84/",
              "gtf/homo_sapiens/Homo_sapiens.GRCh37.82.gtf.gz")
gz <- file.path(tempdir(), "Homo_sapiens.GRCh37.82.gtf.gz")
download.file(url, gz, mode = "wb")

dt <- readGFF(gz) |> as.data.table()
# retain only genes
dt <- dt |> dplyr::filter(type == "gene")
# keep only those genes on a chromosome
dt <- dt[!grepl(pattern = "^G", x = dt$seqid), ]
dt$chr <- paste0("chr", as.character(dt$seqid))
dt <- dt |>
    dplyr::mutate(position = ifelse(strand == "+", start, end)) |>
    dplyr::select(response_id = gene_id,
                  response_name = gene_name,
                  chr, position)
gene_position_data_frame_grch37 <- dt
gene_position_data_frame_grch37$chr <-
    factor(gene_position_data_frame_grch37$chr)
usethis::use_data(gene_position_data_frame_grch37, internal = FALSE,
                  overwrite = TRUE)
