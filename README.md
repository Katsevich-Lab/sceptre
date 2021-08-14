
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sceptre`: robust single-cell CRISPR screen analysis <img src="man/figures/hex.jpg" align="right" alt="" width="180" />
![R-CMD-check](https://github.com/scarlettcanny0629/sceptre/actions/workflows/R-CMD-check.yaml/badge.svg)

<!-- badges: start -->
<<<<<<< HEAD
[![R-CMD-check](https://github.com/scarlettcanny0629/sceptre/workflows/R-CMD-check/badge.svg)](https://github.com/scarlettcanny0629/sceptre/actions)
=======
>>>>>>> d60da2c8f0bd49e89ff231867da570e1120aa74b
<!-- [![R build status](https://travis-ci.com/timothy-barry/sceptre.svg?branch=main)](https://travis-ci.com/timothy-barry/sceptre)
<!-- badges: end -->


Single-cell CRISPR screens provide unprecedented insights into gene
regulation and other facets of human genome biology. However, the
analysis of these screens poses significant statistical and
computational challenges. `sceptre` (pronounced “scepter”) is a
methodology and associated R package for rigorously identifying
regulatory relationships in single-cell CRISPR screen experiments.
`sceptre` tests whether a given perturbation is associated with the
change in expression of a given gene using the intuitive and powerful
conditional randomization test.

# Using the package

You can interact with the package in several ways.

**Demo**: A small demo illustrating the core features of `sceptre` is
available
[here](https://katsevich-lab.github.io/sceptre/articles/sceptre-small-example.html).
We recommend working through the demo to get a quick (5 minutes) feel
for how the method works.

**Moderately-sized analysis**: If you are running an analysis of
moderate size (i.e., the data fit into memory and you are using a single
computer), see the tutorial [running a moderately-sized
analysis](https://katsevich-lab.github.io/sceptre/articles/sceptre-on-moderately-sized-data.html).

**Large-scale analysis**: If you are running a large-scale analysis
(i.e., the data do *not* fit into memory or you are using multiple nodes
on a computer cluster), see the tutorial [running sceptre at
scale](https://katsevich-lab.github.io/sceptre/articles/sceptre-at-scale.html).

**Note**: `sceptre` currently applies to high multiplicity-of-infection
(MOI; \>5 gRNAs/cell) single-cell CRISPR screen data. `sceptre` has not
yet been carefully vetted in low-MOI settings. We are working on
implementing such an extension, but for the time being, please be
cautious in attempting to apply `sceptre` to low-MOI data.

# Installation

You can install the development version of the package from Github with
the following command:

    install.packages("devtools")
    devtools::install_github("katsevich-lab/sceptre")

You can browse the source code on Github
[here](https://github.com/katsevich-lab/sceptre).

Note that sceptre has been tested only in R versions >=3.5 in both macOS systems and Linux systems.

# References

**Methods paper**: T Barry, X Wang, J Morris, K Roeder, E Katsevich.
“Conditional resampling improves calibration and sensitivity in
single-cell CRISPR screen analysis.” Preprint available on
[bioRxiv](https://www.biorxiv.org/content/10.1101/2020.08.13.250092v6).

**Application paper**: J Morris, Z Daniloski, J Domingo, T Barry, M
Ziosi, D Glinos, S Hao, E Mimitou, P Smibert, K Roeder, E Katsevich, T
Lappalainen, N Sanjana. “Discovery of target genes and pathways of blood
trait loci using pooled CRISPR screens and single cell RNA sequencing.”
Preprint available on
[bioRxiv](https://www.biorxiv.org/content/10.1101/2021.04.07.438882v1).

# Funding

We are grateful to [Analytics at
Wharton](https://analytics.wharton.upenn.edu/) for supporting the
development of this software.

<img src="man/figures/wharton_analytics.png" align="center" alt="" width="350" />
