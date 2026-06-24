# Construct *cis* pairs

`construct_cis_pairs()` is a helper function to facilitate construction
the *cis* pairs. `construct_cis_pairs()` returns the set of
target-response pairs for which the target and response are located on
the same chromosome and in close physical proximity to one another.
`construct_cis_pairs()` is a useful pair constructor function for
screens that aim to map noncoding regulatory elements (e.g., enhancers
or noncoding GWAS variants) to target genes in *cis*.
`construct_cis_pairs()` assumes that the columns `chr`, `start`, and
`stop` are present in the `grna_target_data_frame`, giving the
chromosome, start position, and end position, respectively, of the
region that each gRNA targets. `construct_cis_pairs()` takes several
arguments: `sceptre_object` (required), `distance_threshold` (optional),
`positive_control_pairs` (optional), and `response_position_data_frame`
(optional). By default, `construct_cis_pairs()` pairs each gRNA target
to the set of responses on the same chromosome as that target and within
`distance_threshold` bases of that target. (The default value of
`distance_threshold` is 500,000 bases, or half a megabase.) The
`positive_control_pairs` data frame optionally can be passed to
`construct_cis_pairs()`, in which case the positive control targets
(i.e., the entries within the `grna_target` column of
`positive_control_pairs`) are excluded from the *cis* pairs. One may
want to exclude these from the discovery analysis if these targets are
intended for positive control purposes only. See [Section 2.2.2 of the
manual](https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#sec-set-analysis-parameters_construct_cis_pairs)
for more detailed information about this function.

## Usage

``` r
construct_cis_pairs(
  sceptre_object,
  positive_control_pairs = data.frame(),
  distance_threshold = 500000L,
  response_position_data_frame = NULL
)
```

## Arguments

- sceptre_object:

  a `sceptre_object`

- positive_control_pairs:

  (optional) a data frame with columns `grna_target` and `response_id`
  containing the positive control pairs; if supplied, the positive
  control targets are excluded from the *cis* pairs.

- distance_threshold:

  (optional) target-response pairs located within `distance_threshold`
  bases of one another and on the same chromosome are included in the
  *cis* discovery set.

- response_position_data_frame:

  (optional) a data frame with columns `response_id`, `chr`, and
  `position` giving the genomic coordinate of each response; by default,
  the packaged GRCh38 gene position data frame is used.

## Value

a data frame with columns `grna_target` and `response_id` containing the
*cis* pairs

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
```
