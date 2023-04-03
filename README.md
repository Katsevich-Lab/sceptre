
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sceptre`: robust single-cell CRISPR screen analysis <img src="man/figures/hex.jpg" align="right" width="180"/>

<!-- badges: start -->

![](https://github.com/scarlettcanny0629/sceptre/actions/workflows/R-CMD-check.yaml/badge.svg)

<!-- badges: end -->

Single-cell CRISPR screens provide unprecedented insights into gene
regulation and other facets of human genome biology. However, the
analysis of these screens poses significant statistical and
computational challenges. `sceptre` (pronounced “scepter”) is a
statistically principled, fast, memory-light, and user-friendly software
for single-cell CRISPR screen analysis. `sceptre` achieves
state-of-the-art calibration and power on single-cell CRISPR screen data
by leveraging several methodological and algorithmic advances in
assumption-lean and computationally efficient differential expression
analysis.

## Installation

You can install `sceptre` from Github with the following command:

    install.packages("devtools")
    devtools::install_github("katsevich-lab/sceptre")

A Bioconductor upload is forthcoming. `sceptre` has been tested in R
versions \>= 4.1 on macOS and Linux systems.

## Using the software

`sceptre` includes modules for low multiplicity-of-infection (MOI) and
high MOI single-cell CRISPR screen analysis. These modules are based on
distinct statistical methods.

- **Low MOI**: If you are working with low MOI data (\< 2 gRNAs per
  cell), see the [low MOI tutorial](docs/articles/lowmoi_tutorial.html).

- **High MOI**: If you are working with high MOI data (\> 5 gRNAs per
  cell), see the [high MOI
  tutorial](docs/articles/highmoi_tutorial.html).

- **Large-scale analysis at low or high MOI**: if you are running a
  large-scale analysis (i.e., the data do not easily fit into memory or
  you are using a high-performance cluster or cloud), see the `sceptre`
  [Nextflow pipeline](). This pipeline currently is under active
  development and is not yet stable. Open a Github issue or send an
  email if you are interested in using this pipeline, and we will do our
  best to provide support.

You can browse the [source code]() on Github.

## Associated papers

The following series of papers introduces the `sceptre` methodology and
— more broadly — explores statistical and computational challenges at
play in single-cell CRISPR screen analysis.

- Low MOI p-values and calibration: T Barry, K Mason, K Roeder, E
  Katsevich. “Robust differential expression analysis for single-cell
  CRISPR screens.” BioRxiv preprint.

- High MOI p-values and calibration: T Barry, X Wang, J Morris, K
  Roeder, E Katsevich. “SCEPTRE improves calibration and sensitivity in
  single-cell CRISPR screen analysis.” [Genome
  Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02545-2).

- Log fold change estimation: T Barry, E Katsevich, K Roeder.
  “Exponential family measurement error models for single-cell CRISPR
  screens.” [arXiv preprint](https://doi.org/10.48550/arXiv.2201.01879).

Please consider citing the most relevant among these papers if you find
`sceptre` helpful in your research. Please also consider starring this
repository to increase its visibility.

## Funding

We are grateful to [Analytics at
Wharton](https://analytics.wharton.upenn.edu/) for supporting the
development of this software.

<img src="man/figures/wharton_analytics.png" align="center" width="200"/>
