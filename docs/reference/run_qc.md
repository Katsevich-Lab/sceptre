# Run QC

`run_qc()` runs cellwise and pairwise QC on the data. Cellwise QC
involves filtering cells on the covariates `response_n_nonzero`,
`response_n_umis`, and `response_p_mito`. In low-MOI we additionally
remove cells that contain zero or multiple gRNAs. Next, pairwise QC
involves filtering out target-response pairs whose data are too sparse
to be analyzed reliably. In this context we define the “number of
nonzero treatment cells” (resp., the “number of nonzero control cells”)
as the number of cells in the treatment group (resp., control group)
that contain nonzero expression of the response. (We sometimes use the
shorthand `n_nonzero_trt` and `n_nonzero_cntrl` to refer to the number
of nonzero treatment cells and control cells, respectively.) Pairwise QC
involves filtering target-response pairs on `n_nonzero_trt` and
`n_nonzero_cntrl`. See [Chapter 4 of the
manual](https://timothy-barry.github.io/sceptre-book/run-qc.html) for
more detailed information about this function.

## Usage

``` r
run_qc(
  sceptre_object,
  n_nonzero_trt_thresh = 7L,
  n_nonzero_cntrl_thresh = 7L,
  response_n_umis_range = c(0.01, 0.99),
  response_n_nonzero_range = c(0.01, 0.99),
  p_mito_threshold = 0.2,
  additional_cells_to_remove = integer()
)
```

## Arguments

- sceptre_object:

  a `sceptre_object`

- n_nonzero_trt_thresh:

  (optional; default `7L`) an integer specifying the number of nonzero
  *treatment* cells a pair must contain for it to be retained

- n_nonzero_cntrl_thresh:

  (optional; default `7L`) an integer specifying the number of nonzero
  *control* cells a pair must contain for it to be retained

- response_n_umis_range:

  (optional; default `c(0.01, 0.99)`) a length-two vector of percentiles
  specifying the location at which to clip the left and right tails of
  the `response_n_umis` distribution

- response_n_nonzero_range:

  (optional; default `c(0.01, 0.99)`) a length-two vector of percentiles
  specifying the location at which to clip the left and right tails of
  the `response_n_nonzero` distribution

- p_mito_threshold:

  (optional; default `0.2`) a numeric value specifying the location at
  which to clip the right tail of the `response_p_mito` distribution

- additional_cells_to_remove:

  (optional) a vector of integer indices specifying additional cells to
  remove

## Value

an updated `sceptre_object` in which cellwise and pairwise QC have been
applied

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

# set analysis parameters, assign grnas, run qc
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
  ) |>
  assign_grnas(method = "thresholding") |>
  run_qc()
```
