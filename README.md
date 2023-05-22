
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sceptre`: robust analysis of <br>single-cell CRISPR screen data

<!-- badges: start -->

[![R-CMD-check](https://github.com/Katsevich-Lab/sceptre/workflows/R-CMD-check/badge.svg)](https://github.com/Katsevich-Lab/sceptre/actions)

<!-- badges: end -->

<img src="man/figures/hex.jpg" align="right" width="200"/>

Single-cell CRISPR screens provide unprecedented insights into gene
regulation and other facets of human genome biology. However, the
analysis of these screens poses significant statistical and
computational challenges. `sceptre` (pronounced “scepter”) is a
statistically principled, fast, memory-light, and user-friendly software
for single-cell CRISPR screen analysis. `sceptre` achieves
state-of-the-art calibration and power on single-cell CRISPR screen data
by leveraging several methodological and algorithmic advances in
assumption-lean and computationally efficient differential expression
testing.

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

- Large-scale analysis at low or high MOI: If you are running a
  large-scale analysis (i.e., the data do not easily fit into memory or
  you are using a high-performance cluster or cloud), see the `sceptre`
  [Nextflow
  pipeline](https://github.com/timothy-barry/sceptre-pipeline). Open a
  Github issue if you are interested in using this pipeline.

The high-MOI function currently is slower and more memory-intensive than
the low-MOI function. Updates to improve the speed and memory efficiency
of the high-MOI function are coming June 2023.

## Associated papers

The following series of papers introduces the `sceptre` methodology,
applies `sceptre` to discover biology, and — more broadly — interrogates
and aims to resolve statistical and computational challenges at play in
single-cell CRISPR screen analysis.

- (Low MOI analyis) T Barry, K Mason, K Roeder, E Katsevich. “Robust
  differential expression testing for single-cell CRISPR screens at low
  multiplicity of infection.” [BioRxiv
  preprint](https://www.biorxiv.org/content/10.1101/2023.05.15.540875v1), 2023.
  <span style="color:red">(new!)</span>

- (Application of `sceptre` to GWAS loci) J Morris, C Caragine, Z
  Daniloski, J Domingo, T Barry, L Lu, K Davis, M Ziosi, D Glinos, S
  Hao, E Mimitou, P Smibert, K Roeder, E Katsevich, T Lappalainen, N
  Sanjana. “Discovery of target genes and pathways at GWAS loci by
  pooled single-cell CRISPR screens.”
  [Science](https://www.science.org/doi/10.1126/science.adh7699), 2023.
  <span style="color:red">(new!)</span>

- (High MOI analysis) T Barry, X Wang, J Morris, K Roeder, E Katsevich.
  “SCEPTRE improves calibration and sensitivity in single-cell CRISPR
  screen analysis.” [Genome
  Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02545-2), 2021.

- (Effect size estimation) T Barry, K Roeder, E Katsevich. “Exponential
  family measurement error models for single-cell CRISPR screens.”
  [arXiv preprint](https://doi.org/10.48550/arXiv.2201.01879), 2022.

Please consider citing the most relevant among these papers if you find
`sceptre` helpful in your research. Please also consider starring this
repository to increase its visibility. Finally, please contact the
authors if you seek to copy and paste portions of the source code into
your own package.

## Funding

We are grateful to [Analytics at
Wharton](https://analytics.wharton.upenn.edu/) and the [National Science
Foundation](https://www.nsf.gov/) (DMS-2113072) for supporting the
development of this software.

<img src="man/figures/wharton_analytics.png" align="center" width="400"/>
    <img src="man/figures/nsf.jpeg" align="center" width="109"/>
