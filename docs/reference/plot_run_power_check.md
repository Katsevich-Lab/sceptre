# Plot run power check

`plot_run_power_check()` creates a visualization of the outcome of the
power check analysis. Each point in the plot corresponds to a
target-response pair, with positive control pairs in the left column and
negative control pairs in the right column. The vertical axis indicates
the p-value of a given pair; smaller (i.e., more significant) p-values
are positioned higher along this axis (p-values truncated at `clip_to`
for visualization). The positive control p-values should be small, and
in particular, smaller than the negative control p-values.

## Usage

``` r
plot_run_power_check(
  sceptre_object,
  point_size = 1,
  transparency = 0.8,
  clip_to = 1e-20
)
```

## Arguments

- sceptre_object:

  a `sceptre_object` that has had
  [`run_power_check()`](https://katsevich-lab.github.io/sceptre/reference/run_power_check.md)
  called on it

- point_size:

  (optional; default `1`) the size of the individual points in the plot

- transparency:

  (optional; default `0.8`) the transparency of the individual points in
  the plot

- clip_to:

  (optional; default `1e-20`) p-values smaller than this value are set
  to `clip_to` for better visualization. If `clip_to=0` is used then no
  clipping is done.

## Value

a single `ggplot2` plot

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
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
sceptre_object |>
  set_analysis_parameters(
    side = "left",
    positive_control_pairs = positive_control_pairs,
    resampling_mechanism = "permutations",
  ) |>
  assign_grnas(method = "thresholding") |>
  run_qc() |>
  run_calibration_check(
    n_calibration_pairs = 500,
    calibration_group_size = 2
  ) |>
  run_power_check() |>
  plot_run_power_check()
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> Constructing negative control pairs.
#>  ✓
#> Generating permutation resamples.
#>  ✓
#> Analyzing pairs containing response ENSG00000253631 (1 of 96)
#> Analyzing pairs containing response ENSG00000100053 (5 of 96)
#> Analyzing pairs containing response ENSG00000100325 (10 of 96)
#> Analyzing pairs containing response ENSG00000253963 (15 of 96)
#> Analyzing pairs containing response ENSG00000100314 (20 of 96)
#> Analyzing pairs containing response ENSG00000177993 (25 of 96)
#> Analyzing pairs containing response ENSG00000099956 (30 of 96)
#> Analyzing pairs containing response ENSG00000253920 (35 of 96)
#> Analyzing pairs containing response ENSG00000203280 (40 of 96)
#> Analyzing pairs containing response ENSG00000187792 (45 of 96)
#> Analyzing pairs containing response ENSG00000211666 (50 of 96)
#> Analyzing pairs containing response ENSG00000236611 (55 of 96)
#> Analyzing pairs containing response ENSG00000253546 (60 of 96)
#> Analyzing pairs containing response ENSG00000099917 (65 of 96)
#> Analyzing pairs containing response ENSG00000253889 (70 of 96)
#> Analyzing pairs containing response ENSG00000099889 (75 of 96)
#> Analyzing pairs containing response ENSG00000241973 (80 of 96)
#> Analyzing pairs containing response ENSG00000225783 (85 of 96)
#> Analyzing pairs containing response ENSG00000279548 (90 of 96)
#> Analyzing pairs containing response ENSG00000100068 (95 of 96)
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> Generating permutation resamples.
#>  ✓
#> Analyzing pairs containing response ENSG00000224277 (1 of 5)
#> Analyzing pairs containing response ENSG00000226772 (5 of 5)
```
