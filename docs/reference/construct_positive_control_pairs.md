# Construct positive control pairs

`construct_positive_control_pairs()` is a helper function to facilitate
construction of the positive control pairs. Positive control pairs are
target-response pairs for which we know (or have strong reason to
believe) that there is a regulatory relationship between the target and
the response. We can use positive control pairs to verify that `sceptre`
is sensitive (i.e., capable of detecting true associations) on the
dataset under analysis. `construct_positive_control_pairs()` takes as an
argument a `sceptre_object` and returns a data frame with columns
`grna_target` and `response_id`, where gRNA targets and response IDs
with matching names are paired. Typically, the positive control set
consists of transcription start sites paired to the gene regulated by
those transcription start sites. See [Section 2.2 in the
manual](https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#positive-control-pairs)
for more detailed information about this function.

## Usage

``` r
construct_positive_control_pairs(sceptre_object)
```

## Arguments

- sceptre_object:

  a `sceptre_object`

## Value

a data frame with columns `grna_target` and `response_id` containing the
positive control pairs

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
```
