# Plot

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) creates a plot
depicting the current state of a `sceptre_object`.

## Usage

``` r
# S4 method for class 'sceptre_object'
plot(x, y, ...)
```

## Arguments

- x:

  a `sceptre_object`

- y:

  ignored argument

- ...:

  arguments passed to the plotting function dispatched by
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

a single `cowplot` object containing the combined panels (if
`return_indiv_plots` is set to `TRUE`) or a list of the individual
panels (if `return_indiv_plots` is set to `FALSE`)

## Details

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) is "generic" in
the sense that it dispatches a specific plotting function based on the
pipeline function that was most recently called on the `sceptre_object`.
For example, if `run_assign_grnas()` is the most recently called
pipeline function, then
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) dispatches
`plot_run_assign_grnas()`. Similarly, if
[`run_power_check()`](https://katsevich-lab.github.io/sceptre/reference/run_power_check.md)
is the most recently called pipeline function, then
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) dispatches
[`plot_run_power_check()`](https://katsevich-lab.github.io/sceptre/reference/plot_run_power_check.md),
and so on. Users can pass arguments to the function dispatched by
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) as named
arguments to [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

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
# set analysis parameters, assign grnas
sceptre_object <- sceptre_object |>
  set_analysis_parameters() |>
  assign_grnas(method = "thresholding") |>
  plot()
```
