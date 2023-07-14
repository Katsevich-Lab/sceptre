
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sceptre`: robust analysis of <br>single-cell CRISPR screen data

<!-- badges: start -->

[![R-CMD-check](https://github.com/Katsevich-Lab/sceptre/workflows/R-CMD-check/badge.svg)](https://github.com/Katsevich-Lab/sceptre/actions)

<!-- badges: end -->

<img src="man/figures/hex.jpg" align="right" width="200"/>

`sceptre` is a method for single-cell CRISPR screen data analysis.
`sceptre` tests for association between CRISPR perturbations and the
expressions of genes and estimates the effect size of perturbations on
changes in gene expression. `sceptre` achieves state-of-the-art
calibration and power by leveraging several theoretical advances in
assumption-lean and computationally efficient differential expression
analysis.

We pronounce `sceptre` as “scepter,” as in a large, decorated staff. But
any pronunciation will do.

## Installation

You can install `sceptre` from Github with the following command:

    install.packages("devtools")
    devtools::install_github("katsevich-lab/sceptre")

`sceptre` has been tested in R versions \>= 4.1 on macOS and Linux
systems.

## Using the software

`sceptre` includes separate modules for low multiplicity-of-infection
(MOI) and high MOI single-cell CRISPR screen analysis.

- Low MOI: If you are working with low MOI data (1 gRNA per cell), see
  the [low MOI
  tutorial](https://katsevich-lab.github.io/sceptre/articles/lowmoi_tutorial.html).

- High MOI: If you are working with high MOI data ($\geq$ 2 gRNAs per
  cell), see the [high MOI
  tutorial](https://katsevich-lab.github.io/sceptre/articles/highmoi_tutorial.html).

- Large-scale analysis: If you are running a large-scale analysis (i.e.,
  the data do not easily fit into memory or you are using a
  high-performance cluster or cloud), open a Github issue to discuss use
  of the beta `sceptre` [Nextflow
  pipeline](https://github.com/timothy-barry/sceptre-pipeline).

## Funding

We are grateful to [Analytics at
Wharton](https://analytics.wharton.upenn.edu/) and the [National Science
Foundation](https://www.nsf.gov/) (DMS-2113072) for supporting the
development of this software.

<img src="man/figures/wharton_analytics.png" align="center" width="400"/>
    <img src="man/figures/nsf.jpeg" align="center" width="109"/>
