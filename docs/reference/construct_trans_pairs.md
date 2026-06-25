# Construct *trans* pairs

`construct_trans_pairs()` is a helper function to facilitate
construction the set of *trans* pairs. `construct_trans_pairs()` returns
the entire set of possible target-response pairs.
`construct_trans_pairs()` is a useful pair constructor function for
analyses in which we seek to conduct a *trans* analysis, testing each
target against each response. `construct_trans_pairs()` takes as
arguments `sceptre_object` (required), `positive_control_pairs`
(optional), and `pairs_to_exclude` (optional). By default
`construct_trans_pairs()` returns a data frame with columns
`grna_target` and `response_id`, where each gRNA target is mapped to
each response ID.

## Usage

``` r
construct_trans_pairs(
  sceptre_object,
  positive_control_pairs = data.frame(),
  pairs_to_exclude = "none"
)
```

## Arguments

- sceptre_object:

  a `sceptre_object`

- positive_control_pairs:

  (optional) the set of positive control pairs

- pairs_to_exclude:

  (optional; default `"none"`) a string specifying pairs to exclude from
  the *trans* pairs, one of `"none"`, `"pc_pairs"`, or
  `"pairs_containing_pc_targets"`

## Value

a data frame with columns `grna_target` and `response_id` containing the
*trans* discovery set

## Details

The optional argument `pairs_to_exclude` enables the user to remove
specific pairs from the *trans* set and takes values `"none"`,
`"pc_pairs"`, or `"pairs_containing_pc_targets"`. If `pairs_to_exclude`
is set to `"none"` (the default), then no pairs are removed from the
*trans* set. Next, if `pairs_to_exclude` is set to `"pc_pairs"` (and the
`positive_control_pairs` data frame is passed), then then the positive
control target-response pairs are excluded from the *trans* set.
Finally, if `pairs_to_exclude` is set to `"pairs_containing_pc_targets"`
(and `positive_control_pairs` is passed), then *all* pairs containing a
positive control gRNA target are excluded from the *trans* pairs. (In
this sense setting `pairs_to_exclude` to `"pairs_containing_pc_targets"`
is stronger than setting `pairs_to_exclude` to `"pc_pairs"`.) Typically,
in gene-targeting (resp., noncoding-regulatory-element-targeting)
screens, we set `pairs_to_exclude` to `"pc_pairs"` (resp.,
`"pairs_containing_pc_targets"`). See [Section 2.2.2 of the
manual](https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#sec-set-analysis-parameters_construct_trans_pairs)
for more detailed information about this function.

## Examples

``` r
# 1. low-moi, gene-targeting screen
data("lowmoi_example_data")
sceptre_object <- import_data(
  response_matrix = lowmoi_example_data$response_matrix,
  grna_matrix = lowmoi_example_data$grna_matrix,
  extra_covariates = lowmoi_example_data$extra_covariates,
  grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
  moi = "low"
)
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
discovery_pairs <- construct_trans_pairs(
  sceptre_object = sceptre_object,
  positive_control_pairs = positive_control_pairs,
  pairs_to_exclude = "pc_pairs"
)

# 2. high-moi, enhancer-targeting screen
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
discovery_pairs <- construct_trans_pairs(
  sceptre_object = sceptre_object,
  positive_control_pairs = positive_control_pairs,
  pairs_to_exclude = "pairs_containing_pc_targets"
)
```
