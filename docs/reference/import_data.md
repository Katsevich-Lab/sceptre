# Import data

`import_data()` imports data from a collection of R objects to create a
`sceptre_object`. Users can create either a standard `sceptre` object or
an `ondisc`-backed `sceptre` object; the latter is more appropriate for
large-scale data. See [Chapter 1 of the
manual](https://timothy-barry.github.io/sceptre-book/import-data.html#import-data-from-a-collection-of-r-objects)
for more detailed information about this function.

## Usage

``` r
import_data(
  response_matrix,
  grna_matrix,
  grna_target_data_frame,
  moi,
  extra_covariates = data.frame(),
  response_names = NA_character_,
  use_ondisc = FALSE,
  directory_to_write = NULL
)
```

## Arguments

- response_matrix:

  a matrix of response UMI counts, with responses in rows and cells in
  columns. The matrix should be of type `"matrix"`, `"dgCMatrix"`,
  `"dgRMatrix"`, or `"dgTMatrix"`. The row names of the matrix should
  give the response IDs.

- grna_matrix:

  a matrix of gRNA UMI counts, with gRNAs in rows and cells in columns.
  The matrix should be of type `"matrix"`, `"dgCMatrix"`, `"dgRMatrix"`,
  or `"dgTMatrix"`. The row names of the matrix should give the gRNA
  IDs.

- grna_target_data_frame:

  a data frame containing columns `grna_id` and `grna_target` mapping
  each individual gRNA to its target. Non-targeting gRNAs should be
  assigned a label of "non-targeting". Optionally,
  `grna_target_data_frame` can contain columns `chr`, `start`, and
  `end`, giving the chromosome, start coordinate, and end coordiante,
  respectively, of each gRNA. Additionally, `grna_target_data_frame` can
  contain the column `vector_id` specifying the vector to which a given
  gRNA belongs.

- moi:

  a string indicating the MOI of the dataset, either "low" or "high".

- extra_covariates:

  (optional) a data frame containing extra covariates (e.g., batch,
  biological replicate) beyond those that `sceptre` can compute.

- response_names:

  (optional) a vector of human-readable response names; names with the
  prefix "MT-" are taken as mitochondrial genes and are used to compute
  the covariate `response_p_mito`.

- use_ondisc:

  (default `FALSE`) a logical value (i.e., `TRUE` or `FALSE`) indicating
  whether to create an `ondisc`-backed `sceptre_object` (`TRUE`) or a
  standard `sceptre_object` (`FALSE`).

- directory_to_write:

  (optional) a file path to a directory in which to write the backing
  `.odm` files for the response and gRNA expression matrices. Must be
  supplied if `use_ondisc` is set to `TRUE`.

## Value

an initialized `sceptre_object`

## Examples

``` r
data("lowmoi_example_data")
# 1. initialize a standard sceptre_object from R objects
sceptre_object <- import_data(
  response_matrix = lowmoi_example_data$response_matrix,
  grna_matrix = lowmoi_example_data$grna_matrix,
  grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
  extra_covariates = lowmoi_example_data$extra_covariates,
  moi = "low"
)

# 2. initialize an ondisc-backed sceptre_object from R objects
sceptre_object <- import_data(
  response_matrix = lowmoi_example_data$response_matrix,
  grna_matrix = lowmoi_example_data$grna_matrix,
  grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
  extra_covariates = lowmoi_example_data$extra_covariates,
  moi = "low",
  use_ondisc = TRUE,
  directory_to_write = tempdir()
)
```
