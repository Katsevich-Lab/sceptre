# Import data from Cell Ranger

`import_data_from_cellranger()` imports data from the output of one or
more calls to Cell Ranger count. Each directory supplied as an input to
this function should be in feature-barcode format, containing the files
`features.tsv.gz` and `matrix.mtx.gz` (and optionally
`barcodes.tsv.gz`). Users can create either a standard `sceptre` object
or an `ondisc`-backed `sceptre` object; the latter is more appropriate
for large-scale data. See [the introductory
chapter](https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-whole_game_import_data)
or [Chapter
1](https://timothy-barry.github.io/sceptre-book/import-data.html) of the
manual for more information about this function.

## Usage

``` r
import_data_from_cellranger(
  directories,
  moi,
  grna_target_data_frame,
  extra_covariates = data.frame(),
  use_ondisc = FALSE,
  directory_to_write = NULL
)
```

## Arguments

- directories:

  a character vector of file paths to directories containing the output
  of one or more calls to Cell Ranger count. Each directory should
  contain the files `matrix.mtx.gz` and `features.tsv.gz` (and
  optionally `barcodes.tsv.gz`).

- moi:

  a string indicating the MOI of the dataset, either "low" or "high".

- grna_target_data_frame:

  a data frame containing columns `grna_id` and `grna_target` mapping
  each individual gRNA to its target. Non-targeting gRNAs should be
  assigned a label of "non-targeting". Optionally,
  `grna_target_data_frame` can contain columns `chr`, `start`, and
  `end`, giving the chromosome, start coordinate, and end coordiante,
  respectively, of each gRNA. Additionally, `grna_target_data_frame` can
  contain the column `vector_id` specifying the vector to which a given
  gRNA belongs.

- extra_covariates:

  (optional) a data frame containing extra covariates (e.g., batch,
  biological replicate) beyond those that `sceptre` can compute.

- use_ondisc:

  (optional; default `FALSE`) a logical indicating whether to store the
  expression data in a disk-backed `ondisc` matrix (`TRUE`) or an
  in-memory sparse matrix (`FALSE`).

- directory_to_write:

  (optional) a string indicating the directory in which to write the
  backing `.odm` files (must be specified if `use_ondisc` is set to
  `TRUE`).

## Value

an initialized `sceptre_object`

## Examples

``` r
data(grna_target_data_frame_highmoi)
directories <- paste0(
  system.file("extdata", package = "sceptre"),
  "/highmoi_example/gem_group_", c(1, 2)
)

# 1. create a standard sceptre_object from Cell Ranger output
sceptre_object <- import_data_from_cellranger(
  directories = directories,
  moi = "high",
  grna_target_data_frame = grna_target_data_frame_highmoi,
)
#> Processing directory 1.
#>  ✓
#> Processing directory 2.
#>  ✓
#> Combining matrices across directories.
#>  ✓
#> Creating the sceptre object.
#>  ✓

# 2. create an ondisc-backed sceptre_object from Cell Ranger output
sceptre_object <- import_data_from_cellranger(
  directories = directories,
  moi = "high",
  grna_target_data_frame = grna_target_data_frame_highmoi,
  use_ondisc = TRUE,
  directory_to_write = tempdir()
)
#> Round 1/2 processing of the input files.
#>  Processing file 1 of 2.
#>  Processing file 2 of 2.
#> Round 2/2 processing of the input files.
#>  Processing file 1 of 2. Computing cellwise covariates. Writing to disk.
#>  Processing file 2 of 2. Computing cellwise covariates. Writing to disk.
```
