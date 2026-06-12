# Data-generation scripts

This directory documents how the datasets shipped with `sceptre` were
created. Per the Bioconductor guidelines, "Bioconductor requires
documentation of these files in either `inst/scripts/` (preferred) or
`data-raw` directory"
(https://contributions.bioconductor.org/data.html).

Each script records the provenance and parameters of a dataset and
reproduces it; some download large public reference data first.

## Scripts and the objects they produce

### `DATASET_gene_table.R`

Downloads two public reference genomes and produces the gene-position
tables in `data/`:

- `gene_position_data_frame_grch38.rda` -- from the 10x Cell Ranger
    `refdata-gex-GRCh38-2020-A` reference (~11 GB download).
- `gene_position_data_frame_grch37.rda` -- from the Ensembl release-84
    GRCh37 GTF.

### `DATASET_highmoi_example_data.R`

A seeded simulation (`RANDOM_SEED <- 1`) modeled on Gasperini et al.
(2019). Reads `gene_position_data_frame_grch38` and produces:

- `data/highmoi_example_data.rda`
- `data/grna_target_data_frame_highmoi.rda`
- `inst/extdata/highmoi_example/`, in 10x Cell Ranger (v2) format,
    split into `gem_group_1/` and `gem_group_2/`.

### `DATASET_lowmoi_example_data.R`

A seeded simulation (`set.seed(123)`) modeled on Papalexi et al. (2021).
Reads `gene_position_data_frame_grch38` and produces:

- `data/lowmoi_example_data.rda`

### `DATASET_parse_example_data.R`

A seeded subset (`set.seed(4)`) of Parse Biosciences' public "Technical
Performance of CRISPR Detect in Cell Lines" demo (split-pipe v1.0.6p),
written in Parse's native DGE format. Produces the example input for
`import_data_from_parse()` in `inst/extdata/parse_example/`:

- `gene_mat.mtx`, `all_genes.csv` (5000 cells x 266 genes)
- `grna_mat.mtx`, `all_grnas.csv` (5000 cells x 3 guides)

The Parse demo is licensed CC BY 4.0; attribution to Parse Biosciences.
Source: https://www.parsebiosciences.com/datasets/
technical-performance-of-crispr-detect-in-cell-lines/
