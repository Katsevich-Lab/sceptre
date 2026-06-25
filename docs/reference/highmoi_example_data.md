# Example high-MOI data

The example high-MOI CRISPRi data are a small simulated dataset modeled
on that of "A genome-wide framework for mapping gene regulation via
cellular genetic screens" by Gasperini et al., 2019. These data contain
perturbations both of gene transcription start sites (TSSs) and
enhancers, as well as non- targeting perturbations.

## Usage

``` r
data(highmoi_example_data)
```

## Format

An object of class `list` of length 4.

## Details

The example data are stored in a list containing the following
components:

- `response_matrix`: the gene-by-cell expression matrix

- `grna_matrix`: the gRNA-by-cell expression matrix

- `extra_covariates`: a data frame containing a single column, `batch`,
  specifying the batch in which each cell was sequenced (`batch_1` or
  `batch_2`)

- `gene_names`: the human-readable name of each gene
