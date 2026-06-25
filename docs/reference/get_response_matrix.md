# Data getter functions

The data getter functions (i.e., `get_response_matrix()`,
`get_grna_matrix()`, `get_cell_covariates()`) return of a specified data
field from a `sceptre_object`.

## Usage

``` r
get_response_matrix(sceptre_object)

get_grna_matrix(sceptre_object)

get_cell_covariates(sceptre_object)
```

## Arguments

- sceptre_object:

  a `sceptre_object`

## Value

`get_response_matrix()` returns the response matrix contained within a
`sceptre_object`.

`get_grna_matrix()` returns the gRNA matrix contained within a
`sceptre_object`.

`get_cell_covariates()` returns the cell covariate data frame contained
within a `sceptre_object`.

## Examples

``` r
# 1. create a sceptre_object

data(highmoi_example_data)
data(grna_target_data_frame_highmoi)
sceptre_object <- import_data(
  response_matrix = highmoi_example_data$response_matrix,
  grna_matrix = highmoi_example_data$grna_matrix,
  grna_target_data_frame = grna_target_data_frame_highmoi,
  moi = "high",
  extra_covariates = highmoi_example_data$extra_covariates,
  response_names = highmoi_example_data$gene_names
)

# 2. extract data fields from the sceptre_object
response_matrix <- get_response_matrix(sceptre_object)
grna_matrix <- get_grna_matrix(sceptre_object)
cell_covariates <- get_cell_covariates(sceptre_object)
```
