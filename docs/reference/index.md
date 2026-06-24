# Package index

## Package documentation

- [`sceptre`](https://katsevich-lab.github.io/sceptre/reference/sceptre-package.md)
  [`sceptre-package`](https://katsevich-lab.github.io/sceptre/reference/sceptre-package.md)
  : sceptre

## Pipeline functions

The functions that form the basis of the `sceptre` pipeline

- [`import_data()`](https://katsevich-lab.github.io/sceptre/reference/import_data.md)
  : Import data
- [`import_data_from_cellranger()`](https://katsevich-lab.github.io/sceptre/reference/import_data_from_cellranger.md)
  : Import data from Cell Ranger
- [`import_data_from_parse()`](https://katsevich-lab.github.io/sceptre/reference/import_data_from_parse.md)
  : Import data from Parse (experimental)
- [`set_analysis_parameters()`](https://katsevich-lab.github.io/sceptre/reference/set_analysis_parameters.md)
  : Set analysis parameters
- [`assign_grnas()`](https://katsevich-lab.github.io/sceptre/reference/assign_grnas.md)
  : Assign gRNAs to cells
- [`run_qc()`](https://katsevich-lab.github.io/sceptre/reference/run_qc.md)
  : Run QC
- [`run_calibration_check()`](https://katsevich-lab.github.io/sceptre/reference/run_calibration_check.md)
  : Run calibration check
- [`run_power_check()`](https://katsevich-lab.github.io/sceptre/reference/run_power_check.md)
  : Run power check
- [`run_discovery_analysis()`](https://katsevich-lab.github.io/sceptre/reference/run_discovery_analysis.md)
  : Run discovery analysis

## Plotting functions

Functions to visualize the analysis

- [`plot(`*`<sceptre_object>`*`)`](https://katsevich-lab.github.io/sceptre/reference/plot-sceptre_object-method.md)
  : Plot
- [`plot_grna_count_distributions()`](https://katsevich-lab.github.io/sceptre/reference/plot_grna_count_distributions.md)
  : Plot gRNA count distributions
- [`plot_assign_grnas()`](https://katsevich-lab.github.io/sceptre/reference/plot_assign_grnas.md)
  : Plot assign gRNAs
- [`plot_run_qc()`](https://katsevich-lab.github.io/sceptre/reference/plot_run_qc.md)
  : Plot run QC
- [`plot_run_calibration_check()`](https://katsevich-lab.github.io/sceptre/reference/plot_run_calibration_check.md)
  : Plot run calibration check
- [`plot_run_discovery_analysis()`](https://katsevich-lab.github.io/sceptre/reference/plot_run_discovery_analysis.md)
  : Plot run discovery analysis
- [`plot_run_power_check()`](https://katsevich-lab.github.io/sceptre/reference/plot_run_power_check.md)
  : Plot run power check
- [`plot_covariates()`](https://katsevich-lab.github.io/sceptre/reference/plot_covariates.md)
  : Plot covariates
- [`plot_response_grna_target_pair()`](https://katsevich-lab.github.io/sceptre/reference/plot_response_grna_target_pair.md)
  : Plot response-gRNA-target pair

## Pair constructor functions

Functions to facilitate the construction of discovery pairs and positive
control pairs

- [`construct_positive_control_pairs()`](https://katsevich-lab.github.io/sceptre/reference/construct_positive_control_pairs.md)
  : Construct positive control pairs

- [`construct_cis_pairs()`](https://katsevich-lab.github.io/sceptre/reference/construct_cis_pairs.md)
  :

  Construct *cis* pairs

- [`construct_trans_pairs()`](https://katsevich-lab.github.io/sceptre/reference/construct_trans_pairs.md)
  :

  Construct *trans* pairs

- [`gene_position_data_frame_grch38`](https://katsevich-lab.github.io/sceptre/reference/gene_position_data_frame_grch38.md)
  [`gene_position_data_frame_grch37`](https://katsevich-lab.github.io/sceptre/reference/gene_position_data_frame_grch38.md)
  : Gene position data frames

## Information extraction functions

Functions to extract information and results from a `sceptre_object`

- [`print(`*`<sceptre_object>`*`)`](https://katsevich-lab.github.io/sceptre/reference/print-sceptre_object-method.md)
  : Print
- [`write_outputs_to_directory()`](https://katsevich-lab.github.io/sceptre/reference/write_outputs_to_directory.md)
  : Write outputs to directory
- [`get_result()`](https://katsevich-lab.github.io/sceptre/reference/get_result.md)
  : Get result
- [`get_response_matrix()`](https://katsevich-lab.github.io/sceptre/reference/get_response_matrix.md)
  [`get_grna_matrix()`](https://katsevich-lab.github.io/sceptre/reference/get_response_matrix.md)
  [`get_cell_covariates()`](https://katsevich-lab.github.io/sceptre/reference/get_response_matrix.md)
  : Data getter functions
- [`get_grna_assignments()`](https://katsevich-lab.github.io/sceptre/reference/get_grna_assignments.md)
  : Get gRNA assignments

## `ondisc` intput/output functions

- [`read_ondisc_backed_sceptre_object()`](https://katsevich-lab.github.io/sceptre/reference/read_ondisc_backed_sceptre_object.md)
  [`write_ondisc_backed_sceptre_object()`](https://katsevich-lab.github.io/sceptre/reference/read_ondisc_backed_sceptre_object.md)
  :

  Write or read an `ondisc`-backed `sceptre_object`

## Example data

- [`highmoi_example_data`](https://katsevich-lab.github.io/sceptre/reference/highmoi_example_data.md)
  : Example high-MOI data
- [`grna_target_data_frame_highmoi`](https://katsevich-lab.github.io/sceptre/reference/grna_target_data_frame_highmoi.md)
  : gRNA target data frame
- [`lowmoi_example_data`](https://katsevich-lab.github.io/sceptre/reference/lowmoi_example_data.md)
  : Example low-MOI data
