# Assign gRNAs to cells

`assign_grnas()` performs the gRNA-to-cell assignments. `sceptre`
provides three gRNA-to-cell assignment strategies: the mixture method,
the thresholding method, and the maximum method. The mixture method
involves assigning gRNAs to cells using a principled mixture model.
Next, the thresholding method assigns a gRNA to a cell if the UMI count
of the gRNA in the cell is greater than or equal to some integer
threshold. Finally, the maximum method assigns the gRNA that accounts
for the greatest number of UMIs in a given cell to that cell. The
maximum method is available only in low MOI. See [Chapter 3 of the
manual](https://timothy-barry.github.io/sceptre-book/assign-grnas.html)
for more detailed information about `assign_grnas()`.

## Usage

``` r
assign_grnas(
  sceptre_object,
  method = "default",
  print_progress = TRUE,
  parallel = FALSE,
  n_processors = "auto",
  log_dir = tempdir(),
  ...
)
```

## Arguments

- sceptre_object:

  a `sceptre_object`

- method:

  (optional) a string indicating the method to use to assign the gRNAs
  to cells, one of `"mixture"`, `"thresholding"`, or `"maximum"`. The
  default is `"maximum"` in low MOI and `"mixture"` in high MOI.

- print_progress:

  (optional; default `TRUE`) a logical indicating whether to print
  progress updates

- parallel:

  (optional; default `FALSE`) a logical indicating whether to run the
  function in parallel. `parallel = TRUE` is recommended only on Mac; it
  is not supported on Windows and may behave unreliably on Linux
  clusters.

- n_processors:

  (optional; default `"auto"`) an integer specifying the number of
  processors to use if `parallel` is set to `TRUE`. The default,
  `"auto"`, uses half the physical cores. The fraction may be tuned via
  the `parallelly.availableCores.fraction` R option.

- log_dir:

  (optional; default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html)) a string
  indicating the directory in which to write the log files (ignored if
  `parallel = FALSE`)

- ...:

  optional method-specific additional arguments

## Value

an updated `sceptre_object` in which the gRNA assignments have been
carried out

## Note

See the manual for information about the method-specific additional
arguments.

## Examples

``` r
data("lowmoi_example_data")
# 1. import data, set default analysis parameters
sceptre_object <- import_data(
  response_matrix = lowmoi_example_data$response_matrix,
  grna_matrix = lowmoi_example_data$grna_matrix,
  extra_covariates = lowmoi_example_data$extra_covariates,
  grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
  moi = "low"
) |> set_analysis_parameters()

# 2. assign gRNAs (three different methods)
sceptre_object <- sceptre_object |> assign_grnas(method = "thresholding")
sceptre_object <- sceptre_object |> assign_grnas(method = "maximum")
sceptre_object <- sceptre_object |> assign_grnas(method = "mixture")
#> Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.
#> Performing gRNA-to-cell assignments for gRNA ENSG00000182704_grna1 (1 of 50)
#> Performing gRNA-to-cell assignments for gRNA ENSG00000287679_grna1 (5 of 50)
#> Performing gRNA-to-cell assignments for gRNA ENSG00000257275_grna2 (10 of 50)
#> Performing gRNA-to-cell assignments for gRNA ENSG00000242110_grna1 (15 of 50)
#> Performing gRNA-to-cell assignments for gRNA ENSG00000224311_grna2 (20 of 50)
#> Performing gRNA-to-cell assignments for gRNA ENSG00000260303_grna1 (25 of 50)
#> Performing gRNA-to-cell assignments for gRNA ENSG00000181374_grna2 (30 of 50)
#> Performing gRNA-to-cell assignments for gRNA ENSG00000169836_grna1 (35 of 50)
#> Performing gRNA-to-cell assignments for gRNA ENSG00000109606_grna2 (40 of 50)
#> Performing gRNA-to-cell assignments for gRNA nt_grna5 (45 of 50)
#> Performing gRNA-to-cell assignments for gRNA nt_grna10 (50 of 50)
```
