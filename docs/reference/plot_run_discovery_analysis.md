# Plot run discovery analysis

`plot_run_discovery_analysis()` creates a visualization of the outcome
of the discovery analysis. The visualization consists of four plots:

- The upper left plot superimposes the discovery p-values (blue) on top
  of the negative control p-values (red) on an untransformed scale.

- The upper right plot is the same as the upper left plot, but the scale
  is negative log-10 transformed. The discovery p-values ideally should
  trend above the diagonal line, indicating the presence of signal in
  the discovery set. The horizontal dashed line indicates the multiple
  testing threshold; discovery pairs whose p-value falls above this line
  are called as significant.

- The bottom left panel is a volcano plot of the p-values and log fold
  changes of the discovery pairs. Each point corresponds to a pair; the
  estimated log-2 fold change of the pair is plotted on the horizontal
  axis, and the (negative log-10 transformed) p-value is plotted on the
  vertical axis. The horizontal dashed line again indicates the multiple
  testing threshold. Points above the dashed line (colored in purple)
  are called as discoveries, while points below (colored in blue) are
  called as insignificant.

- The bottom right panel is a text box displaying the number of
  discovery pairs called as significant.

## Usage

``` r
plot_run_discovery_analysis(
  sceptre_object,
  x_limits = c(-1.5, 1.5),
  point_size = 0.55,
  transparency = 0.8,
  return_indiv_plots = FALSE
)
```

## Arguments

- sceptre_object:

  a `sceptre_object` that has had `run_discovery_analysis` called on it

- x_limits:

  (optional; default `c(-1.5, 1.5)`) a numeric vector of length 2 giving
  the lower and upper limits of the x-axis (corresponding to log-2 fold
  change) for the "Discovery volcano plot" panel

- point_size:

  (optional; default `0.55`) the size of the individual points in the
  plot

- transparency:

  (optional; default `0.8`) the transparency of the individual points in
  the plot

- return_indiv_plots:

  (optional; default `FALSE`) if `FALSE` then a list of `ggplot` is
  returned; if `TRUE` then a single `cowplot` object is returned.

## Value

a single `cowplot` object containing the combined panels (if
`return_indiv_plots` is set to `TRUE`) or a list of the individual
panels (if `return_indiv_plots` is set to `FALSE`)

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
discovery_pairs <- construct_cis_pairs(sceptre_object,
  positive_control_pairs = positive_control_pairs,
  distance_threshold = 5e6
)
sceptre_object |>
  set_analysis_parameters(
    side = "left",
    discovery_pairs = discovery_pairs,
    resampling_mechanism = "permutations",
  ) |>
  assign_grnas(method = "thresholding") |>
  run_qc() |>
  run_calibration_check() |>
  run_discovery_analysis() |>
  plot_run_discovery_analysis()
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> Constructing negative control pairs.
#>  ✓
#> Generating permutation resamples.
#>  ✓
#> Analyzing pairs containing response ENSG00000253631 (1 of 98)
#> Analyzing pairs containing response ENSG00000100053 (5 of 98)
#> Analyzing pairs containing response ENSG00000100325 (10 of 98)
#> Analyzing pairs containing response ENSG00000253963 (15 of 98)
#> Analyzing pairs containing response ENSG00000100314 (20 of 98)
#> Analyzing pairs containing response ENSG00000177993 (25 of 98)
#> Analyzing pairs containing response ENSG00000099956 (30 of 98)
#> Analyzing pairs containing response ENSG00000253920 (35 of 98)
#> Analyzing pairs containing response ENSG00000203280 (40 of 98)
#> Analyzing pairs containing response ENSG00000187792 (45 of 98)
#> Analyzing pairs containing response ENSG00000211666 (50 of 98)
#> Analyzing pairs containing response ENSG00000236611 (55 of 98)
#> Analyzing pairs containing response ENSG00000253546 (60 of 98)
#> Analyzing pairs containing response ENSG00000099917 (65 of 98)
#> Analyzing pairs containing response ENSG00000253889 (70 of 98)
#> Analyzing pairs containing response ENSG00000099889 (75 of 98)
#> Analyzing pairs containing response ENSG00000241973 (80 of 98)
#> Analyzing pairs containing response ENSG00000225783 (85 of 98)
#> Analyzing pairs containing response ENSG00000279548 (90 of 98)
#> Analyzing pairs containing response ENSG00000100068 (95 of 98)
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> Generating permutation resamples.
#>  ✓
#> Analyzing pairs containing response ENSG00000100218 (1 of 97)
#> Analyzing pairs containing response ENSG00000099958 (5 of 97)
#> Analyzing pairs containing response ENSG00000211638 (10 of 97)
#> Analyzing pairs containing response ENSG00000211685 (15 of 97)
#> Analyzing pairs containing response ENSG00000220891 (20 of 97)
#> Analyzing pairs containing response ENSG00000187905 (25 of 97)
#> Analyzing pairs containing response ENSG00000235954 (30 of 97)
#> Analyzing pairs containing response ENSG00000244486 (35 of 97)
#> Analyzing pairs containing response ENSG00000225783 (40 of 97)
#> Analyzing pairs containing response ENSG00000224277 (45 of 97)
#> Analyzing pairs containing response ENSG00000286326 (50 of 97)
#> Analyzing pairs containing response ENSG00000211672 (55 of 97)
#> Analyzing pairs containing response ENSG00000203280 (60 of 97)
#> Analyzing pairs containing response ENSG00000233521 (65 of 97)
#> Analyzing pairs containing response ENSG00000100319 (70 of 97)
#> Analyzing pairs containing response ENSG00000211674 (75 of 97)
#> Analyzing pairs containing response ENSG00000099917 (80 of 97)
#> Analyzing pairs containing response ENSG00000100053 (85 of 97)
#> Analyzing pairs containing response ENSG00000229770 (90 of 97)
#> Analyzing pairs containing response ENSG00000253920 (95 of 97)
```
