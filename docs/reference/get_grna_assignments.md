# Get gRNA assignments

`get_grna_assignments()` returns the gRNA-to-cell assignments contained
within a `sceptre_object`. The output is a sparse logical matrix, with
gRNAs in the rows and cells in the columns. A given entry of the matrix
is set to `TRUE` if the given gRNA is assigned to the given cell (and
`FALSE` otherwise).

## Usage

``` r
get_grna_assignments(sceptre_object, apply_cellwise_qc = FALSE)
```

## Arguments

- sceptre_object:

  a `sceptre_object` that has had
  [`assign_grnas()`](https://katsevich-lab.github.io/sceptre/reference/assign_grnas.md)
  called on it

- apply_cellwise_qc:

  a logical value (i.e., `TRUE` or `FALSE`) indicating whether to return
  the gRNA-to-cell assignment matrix after cellwise QC has been applied
  (default `FALSE`)

## Value

a sparse logical matrix containing the gRNA-to-cell assignments

## Details

When using the "maximum" assignment strategy, exactly one gRNA is
assigned to a given cell. In other words, each column of the
gRNA-to-cell assignment matrix contains exactly one TRUE entry.

## Examples

``` r
data(highmoi_example_data)
data(grna_target_data_frame_highmoi)
# import data
sceptre_object <- import_data(
  response_matrix = highmoi_example_data$response_matrix,
  grna_matrix = highmoi_example_data$grna_matrix,
  grna_target_data_frame = grna_target_data_frame_highmoi,
  moi = "high",
  extra_covariates = highmoi_example_data$extra_covariates,
  response_names = highmoi_example_data$gene_names
)
discovery_pairs <- construct_cis_pairs(sceptre_object)
sceptre_object <- sceptre_object |>
  set_analysis_parameters(
    discovery_pairs = discovery_pairs,
    side = "left"
  ) |>
  assign_grnas(method = "mixture") |>
  run_qc()
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> Performing gRNA-to-cell assignments for gRNA ENSG00000224277_grna1 (1 of 60)
#> Performing gRNA-to-cell assignments for gRNA ENSG00000226772_grna1 (5 of 60)
#> Performing gRNA-to-cell assignments for gRNA ENSG00000286326_grna2 (10 of 60)
#> Performing gRNA-to-cell assignments for gRNA candidate_enh_3_grna1 (15 of 60)
#> Performing gRNA-to-cell assignments for gRNA candidate_enh_5_grna2 (20 of 60)
#> Performing gRNA-to-cell assignments for gRNA candidate_enh_8_grna1 (25 of 60)
#> Performing gRNA-to-cell assignments for gRNA candidate_enh_10_grna2 (30 of 60)
#> Performing gRNA-to-cell assignments for gRNA candidate_enh_13_grna1 (35 of 60)
#> Performing gRNA-to-cell assignments for gRNA candidate_enh_15_grna2 (40 of 60)
#> Performing gRNA-to-cell assignments for gRNA candidate_enh_18_grna1 (45 of 60)
#> Performing gRNA-to-cell assignments for gRNA candidate_enh_20_grna2 (50 of 60)
#> Performing gRNA-to-cell assignments for gRNA non-targeting_grna5 (55 of 60)
#> Performing gRNA-to-cell assignments for gRNA non-targeting_grna10 (60 of 60)

grna_assignment_matrix <- get_grna_assignments(
  sceptre_object = sceptre_object
)
grna_assignment_matrix_with_qc <- get_grna_assignments(
  sceptre_object = sceptre_object,
  apply_cellwise_qc = TRUE
)
```
