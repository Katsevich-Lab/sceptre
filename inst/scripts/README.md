# Data-generation scripts

This directory documents how the datasets shipped with `sceptre` were
created. Per the Bioconductor guidelines, "Bioconductor requires
documentation of these files in either `inst/scripts/` (preferred) or
`data-raw` directory"
(https://contributions.bioconductor.org/data.html).

The scripts are documentary. They depend on external reference files and
on helper functions not exported by `sceptre`, so they are not expected
to run as-is in a user session; they record the provenance and
parameters of each dataset.

## Scripts and the objects they produce

### `DATASET_gene_table.R`

Produces the gene-position reference tables in `data/`:

- `gene_position_data_frame_grch38.rda`
- `gene_position_data_frame_grch37.rda`

External inputs (not shipped with the package):

- GRCh38: `genes.gtf` from 10x Cell Ranger `refdata-gex-GRCh38-2020-A`.
- GRCh37: `Homo_sapiens.GRCh37.82.gtf.gz` from Ensembl (release-84).

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

## Not yet covered by a script

`inst/extdata/parse_example/` is example input for
`import_data_from_parse()`. It was committed directly and has no
generation script here; its provenance is not yet documented.
