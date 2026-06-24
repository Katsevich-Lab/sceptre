# Print

[`print()`](https://rdrr.io/r/base/print.html) prints information about
the dataset and the status of the analysis to the console. The output
contains several fields: Attributes of the data (summarizing key
features of the data), Analysis status (indicating the analysis
functions that have been called), Analysis parameters (summarizing the
analysis parameters set in
[`set_analysis_parameters()`](https://katsevich-lab.github.io/sceptre/reference/set_analysis_parameters.md)),
gRNA-to-cell assignment information (summarizing the outcome of the
gRNA-to-cell assignment step), and Summary of results (summarizing the
key analysis results). A subset of these fields may be printed,
depending on the status of the analysis.

## Usage

``` r
# S4 method for class 'sceptre_object'
print(x)
```

## Arguments

- x:

  a `sceptre_object`

## Value

the value NULL

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
print(sceptre_object)
#> An object of class sceptre_object.
#> 
#> Attributes of the data:
#>  • 500 cells
#>  • 100 responses
#>  • High multiplicity-of-infection 
#>  • 50 targeting gRNAs (distributed across 25 targets) 
#>  • 10 non-targeting gRNAs 
#>  • 5 covariates (batch, grna_n_nonzero, grna_n_umis, response_n_nonzero, response_n_umis)
#> 
#> Analysis status:
#>  ✓ import_data()
#>  ✗ set_analysis_parameters()
#>  ✗ assign_grnas()
#>  ✗ run_qc()
#>  ✗ run_calibration_check()
#>  ✗ run_power_check()
#>  ✗ run_discovery_analysis()
#> 
#> Analysis parameters: 
#>  • Discovery pairs: not specified
#>  • Positive control pairs: not specified
#>  • Sidedness of test: not specified
#>  • Resampling mechanism: not specified
#>  • gRNA integration strategy: not specified
#>  • Resampling approximation: not specified
#>  • Multiple testing adjustment: none
#>  • N nonzero treatment cells threshold: not specified
#>  • N nonzero control cells threshold: not specified
#>  • Formula object: not specified
```
