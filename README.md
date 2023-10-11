
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/hex.jpg" align="right" width="180"/>

## `sceptre`: Single-Cell CRISPR Screen<br>Data Analysis

[![R-CMD-check](https://github.com/Katsevich-Lab/sceptre/workflows/R-CMD-check/badge.svg)](https://github.com/Katsevich-Lab/sceptre/actions)

Single-cell CRISPR screens (e.g., Perturb-seq, TAP-seq) pose enormous
potential for advancing understanding of human disease biology. However,
the analysis of these screens presents considerable statistical and
computational challenges. `sceptre` is an R package for single-cell
CRISPR screen data analysis, emphasizing statistical rigor,
computational efficiency, and ease of use.

## Key features

`sceptre` offers the following key features:

- Import data from 10X CellRanger or a set of R matrices
- Assign gRNAs to cells using one of three principled methods
- Perform extensive quality control
- Run gRNA-to-gene differential expression analyses
- Verify false discovery rate control using negative control gRNAs
- Verify adequate power using positive control gRNAs
- Visualize each step in the pipeline by rendering an informative plot

`sceptre` is compatible with a broad range of single-cell CRISPR screen
experimental designs:

- low multiplicity-of-infection and high multiplicity-of-infection
- Gene-targeting and noncoding-regulatory-element-targeting
- CRISPRko, CRISPRi, CRISPRa, CRISPR base editing, and CRISPR prime
  editing
- Gene and protein expression readout
- Direct gRNA capture and barcoded gRNA capture

## Get started

The fastest way to get started with `sceptre` is to work through the Get
Started vignette (10 minutes): `vignette("sceptre")`.

## Installation instructions

Users can install `sceptre` with the following command:

    # R version >= 4.1
    install.packages("devtools")
    devtools::install_github("katsevich-lab/sceptre")

`sceptre` will be uploaded to Bioconductor in the near future.

## Featured publications

`sceptre` has been featured in the following publications:

- [Morris et al.,
  2023](https://www.science.org/doi/10.1126/science.adh7699). “Discovery
  of target genes and pathways…”. *Science*.
- [Barry et al.,
  2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02545-2).
  “SCEPTRE improves calibration and sensitivity…”. *Genome Biology*.

See publication list for a full list of publications.

## Bugs and feature requests

We ask that users submit bugs and feature requests as issues to the
[Github repository](https://github.com/Katsevich-Lab/sceptre/issues).

## Funding

[Analytics at Wharton](https://analytics.wharton.upenn.edu/) and the
[National Science Foundation](https://www.nsf.gov/) (DMS-2113072)
support the development of this software.

<img src="man/figures/wharton_analytics.png" align="center" width="200"/>
    <img src="man/figures/nsf.jpeg" align="center" width="109"/>
