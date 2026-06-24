# Example low-MOI data

The example low-MOI data are a small simulated dataset modeled on that
in the paper "Characterizing the molecular regulation of inhibitory
immune checkpoints with multimodal single-cell screens" by Papalexi et
al., 2021. This dataset includes gene-targeting CRISPRko perturbations
as well as non- targeting perturbations.

## Usage

``` r
data(lowmoi_example_data)
```

## Format

An object of class `list` of length 4.

## Details

The example data are stored in a list containing four components:

- `response_matrix`: the gene-by-cell expression matrix

- `grna_matrix`: the gRNA-by-cell expression matrix

- `grna_target_data_frame`: a data frame containing the columns
  `grna_id` (ID of an individual gRNA) and `grna_target` (genomic target
  of the gRNA)

- `extra_covariates`: a data frame with a single column, `batch`,
  specifying the batch in which each cell was sequenced (`batch_1` or
  `batch_2`)
