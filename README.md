
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sceptre`: robust single-cell CRISPR screen analysis <img src="man/figures/hex.jpg" align="right" alt="" width="180" />

<!-- badges: start -->

![R-CMD-check](https://github.com/scarlettcanny0629/sceptre/actions/workflows/R-CMD-check.yaml/badge.svg)
<!-- badges: end -->

Single-cell CRISPR screens provide unprecedented insights into gene
regulation and other facets of human genome biology. However, the
analysis of these screens poses significant statistical and
computational challenges. `sceptre` (pronounced “scepter”) is a
methodology and associated R package for rigorously identifying
regulatory relationships in single-cell CRISPR screen experiments.
`sceptre` tests whether a given perturbation is associated with the
change in expression of a given gene using the robust, powerful, and
intuitive conditional randomization test.

**Update March 2022**: We are excited to release `sceptre` version
0.1.0, a major update that significantly improves the speed and
ease-of-use of the software. Please download the latest version of
`sceptre` (see below) and check the updated tutorial and news page for
further details.

# Installation

You can install the development version of the package from Github with
the following command:

    install.packages("devtools")
    devtools::install_github("katsevich-lab/sceptre")

You can browse the source code on Github
[here](https://github.com/katsevich-lab/sceptre). `sceptre` has been
tested in R versions \>=3.5 on macOS and Linux systems.

# Using the software

`sceptre` has several interfaces, which you can choose between based on
the size of your analysis.

**Small or moderately-sized analysis**: If you are running an analysis
of small or moderate size (i.e., the data fit into memory and you are
using a single computer), see the standard `sceptre` tutorial
[here](https://katsevich-lab.github.io/sceptre/articles/using_sceptre_v2.html).

**Large-scale analysis**: If you are running a large-scale analysis
(i.e., the data do not fit into memory or you are using a
high-performance cluster or cloud), see the at-scale tutorial here.
<span style="color:blue">(Currently under construction; will be
available soon.)</span>

**Note**: `sceptre` currently applies to high multiplicity-of-infection
(MOI; \>5 gRNAs/cell) single-cell CRISPR screen data. `sceptre` has not
yet been carefully vetted in low-MOI settings. We are working on
developing such an extension, which we expect to be available in 2022.

# References

Please consider starring this repository and citing the following if you
find `sceptre` helpful in your research.

**Methods papers**

T Barry, X Wang, J Morris, K Roeder, E Katsevich. “SCEPTRE improves
calibration and sensitivity in single-cell CRISPR screen analysis.”
[Genome
Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02545-2).

T Barry, E Katsevich, K Roeder. “Exponential family measurement error
models for single-cell CRISPR screens.” [arXiv
preprint](https://doi.org/10.48550/arXiv.2201.01879).

**Application paper**

J Morris, Z Daniloski, J Domingo, T Barry, M Ziosi, D Glinos, S Hao, E
Mimitou, P Smibert, K Roeder, E Katsevich, T Lappalainen, N Sanjana.
“Discovery of target genes and pathways of blood trait loci using pooled
CRISPR screens and single cell RNA sequencing.” Preprint available on
[bioRxiv](https://www.biorxiv.org/content/10.1101/2021.04.07.438882v1).

# Funding

We are grateful to [Analytics at
Wharton](https://analytics.wharton.upenn.edu/) for supporting the
development of this software.

<img src="man/figures/wharton_analytics.png" align="center" alt="" width="350" />
