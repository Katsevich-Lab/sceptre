
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `sceptre`: robust single-cell CRISPR screen analysis <img src="man/figures/hex.jpg" align="right" alt="" width="180" />

<!-- badges: start -->

![R-CMD-check](https://github.com/scarlettcanny0629/sceptre/actions/workflows/R-CMD-check.yaml/badge.svg)
<!-- badges: end -->

Single-cell CRISPR screens provide unprecedented insights into gene
regulation and other facets of human genome biology. However, the
analysis of these screens poses significant statistical and
computational challenges. `sceptre` (pronounced “scepter”) is a
statistically principled, fast, memory-light, and user-friendly software
for single-cell CRISPR screen analysis. `sceptre` is powered by several
statistical and algorithmic innovations in assumption-lean and
computationally efficient differential expression analysis.

# Installation

You can install the `sceptre` package from Github with the following
command:

    install.packages("devtools")
    devtools::install_github("katsevich-lab/sceptre")

`sceptre` has been tested in R versions \>= 4.1 on macOS and Linux
systems.

# Using the software

`sceptre` includes modules for low multiplicity of infection (MOI) and
high MOI single-cell CRISPR screen analysis. These modules are based
distinct statistical methods.

**Low multiplicity of infection**:

**High multiplicity of infection**:

**Large-scale analysis**:

**Large-scale analysis**: If you are running a large-scale analysis
(i.e., the data do not easily fit into memory or you are using a
high-performance cluster or cloud), see the `sceptre` Nextflow pipeline
[here](https://github.com/timothy-barry/sceptre-pipeline). The
documentation for the `sceptre` Nextflow pipeline currently is sparse;
please open a Github issue if you are interested in using this pipeline,
and we will provide support.

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

### Funding

We are grateful to [Analytics at
Wharton](https://analytics.wharton.upenn.edu/) for supporting the
development of this software.

<img src="man/figures/wharton_analytics.png" align="center" width="200" />
