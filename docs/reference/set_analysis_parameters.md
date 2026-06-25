# Set analysis parameters

`set_analysis_parameters()` sets the analysis parameters that control
how the statistical analysis is to be conducted. See [Chapter 2 of the
manual](https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html)
for more detailed information about this function.

## Usage

``` r
set_analysis_parameters(
  sceptre_object,
  discovery_pairs = data.frame(grna_target = character(0), response_id = character(0)),
  positive_control_pairs = data.frame(grna_target = character(0), response_id =
    character(0)),
  side = "both",
  grna_integration_strategy = "union",
  formula_object = "default",
  resampling_approximation = "skew_normal",
  control_group = "default",
  resampling_mechanism = "default",
  multiple_testing_method = "BH",
  multiple_testing_alpha = 0.1
)
```

## Arguments

- sceptre_object:

  a `sceptre_object`

- discovery_pairs:

  (optional) a data frame with columns `grna_target` and `response_id`
  specifying the discovery pairs to analyze

- positive_control_pairs:

  (optional) a data frame with columns `grna_target` and `response_id`
  specifying the positive control pairs to analyze

- side:

  (optional; default `"both"`) the sidedness of the test, one of
  `"left"`, `"right"`, or `"both"`

- grna_integration_strategy:

  (optional; default `"union"`) a string specifying the gRNA integration
  strategy, either `"singleton"`, `"union"`, or `"bonferroni"`

- formula_object:

  (optional) a formula object specifying how to adjust for the
  covariates in the model

- resampling_approximation:

  (optional; default `"skew_normal"`) a string indicating the resampling
  approximation to make to the null distribution of test statistics,
  either `"skew_normal"` or `"no_approximation"`

- control_group:

  (optional) a string specifying the control group to use in the
  differential expression analysis, either `"complement"` or
  `"nt_cells"`

- resampling_mechanism:

  (optional) a string specifying the resampling mechanism to use, either
  `"permutations"` or `"crt"`

- multiple_testing_method:

  (optional; default `"BH"`) a string specifying the multiple testing
  correction method to use; see `p.adjust.methods` for options

- multiple_testing_alpha:

  (optional; default `0.1`) a numeric specifying the nominal level of
  the multiple testing correction method

## Value

an updated `sceptre_object` in which the analysis parameters have been
set

## Note

Every argument to this function is optional, but typically, users want
to specify `discovery_pairs` at minimum.

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

# set analysis parameters
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
discovery_pairs <- construct_cis_pairs(sceptre_object,
  positive_control_pairs = positive_control_pairs,
  distance_threshold = 5e6
)
sceptre_object <- sceptre_object |>
  set_analysis_parameters(
    discovery_pairs = discovery_pairs,
    positive_control_pairs = positive_control_pairs,
    side = "left"
  )
```
