
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sceptre`: robust single-cell CRISPR screen analysis <img src="man/figures/hex.jpg" align="right" alt="" width="180" />

<!-- badges: start -->

<!-- badges: end -->

Single-cell CRISPR screens provide unprecedented insight into gene
regulation and other facets of human genome biology. However, the
analysis of these screens poses significant statistical and
computational challenges. `sceptre` (pronounced “scepter”) is a
methodology and associated R package for rigorously identifying
regulatory relationships in high multiplicity-of-infection single-cell
CRISPR screens. For a given gRNA and gene, `sceptre` produces a
well-calibrated *p*-value for the null hypothesis that the gRNA has no
impact on the expression of the gene. `sceptre` is robust and reliable:
it seamlessly accounts for technical factors (e.g., sequencing batch and
sequencing depth) and works well even when the expression model is
misspecified.

Tutorials available on the [package
website](https://timothy-barry.github.io/sceptre/) demonstrate how to
apply `sceptre` to small and
[large](https://timothy-barry.github.io/sceptre/articles/sceptre-at-scale.html)
datasets.

# Installation

You can install the development version of the package from Github with
the following command:

    install.packages("devtools")
    devtools::install_github("timothy-barry/sceptre")
