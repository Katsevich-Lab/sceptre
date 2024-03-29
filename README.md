
<!-- README.md is generated from README.Rmd. Please edit that file -->

<div style="margin-top: 5px;">

<img src="man/figures/hex.jpg" align="right" width="150"/>

</div>

## <span style="font-size:60px;">`sceptre`</span>

An R package for single-cell CRISPR screen data analysis, emphasizing
statistical rigor, computational efficiency, and ease of use.

<!-- badges: start -->

[![R-CMD-check](https://github.com/Katsevich-Lab/sceptre/workflows/R-CMD-check/badge.svg)](https://github.com/Katsevich-Lab/sceptre/actions)

<!-- badges: end -->

## Release of `sceptre` v0.9.2

We are excited to announce the release of `sceptre` v0.9.2, a
substantial update to the package. This milestone includes the following
developments:

- A reimagined user experience based on a modular, object-oriented
  workflow.
- Further improvements in speed and memory efficiency.
- A unified interface for low- and high-MOI analyses.
- A suite of plotting functions facilitating visualization of each step
  in the pipeline.
- Expanded support for gRNA assignment and quality control.
- Interoperability with output from 10X Cell Ranger and Parse
  Biosciences CRISPR Detect.
- An [e-book](https://timothy-barry.github.io/sceptre-book/) guiding
  users through the entire process of analyzing their data using
  `sceptre`.

`sceptre` v0.9.2 facilitates an entire analysis pipeline for single-cell
CRISPR screens, starting from UMI count data obtained from tools like
10X Cell Ranger.

<!--
&#10;## About `sceptre`
&#10;Single-cell CRISPR screens (e.g., Perturb-seq, TAP-seq) combine CRISPR and single-cell sequencing to survey the effects of genetic perturbations on individual cells. Despite their promise, single-cell CRISPR screens present considerable statistical and computational challenges. `sceptre` is an R package for single-cell CRISPR screen data analysis, emphasizing statistical rigor, computational efficiency, and ease of use.
&#10;### Key features
&#10;-   Import data from 10X CellRanger or a set of R matrices
-   Assign gRNAs to cells using one of three principled methods
-   Perform extensive quality control
-   Run gRNA-to-gene differential expression analyses
-   Verify false discovery rate control using negative control gRNAs
-   Verify adequate power using positive control gRNAs
-   Visualize each step in the pipeline by rendering an informative plot
&#10;### Compatible experimental designs
&#10;-   low multiplicity-of-infection and high multiplicity-of-infection
-   Gene-targeting and noncoding-regulatory-element-targeting
-   CRISPRko, CRISPRi, CRISPRa, CRISPR base editing, and CRISPR prime editing
-   Gene and protein expression readout
&#10;## Get started
&#10;The fastest way to get started with `sceptre` is to work through the [Get Started vignette](articles/sceptre.html) (30 minutes). Users also can read the `sceptre` [manual](https://timothy-barry.github.io/sceptre-book/).
-->

## Featured publications

- [Morris et al.,
  2023](https://www.science.org/doi/10.1126/science.adh7699). “Discovery
  of target genes and pathways…”. *Science*.
- [Barry et al.,
  2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02545-2).
  “SCEPTRE improves calibration and sensitivity…”. *Genome Biology*.

## Bug reports, feature requests, and software questions

For bug reports, please open a [GitHub
issue](https://github.com/Katsevich-Lab/sceptre/issues). For feature
requests, please start a discussion under [feature
requests](https://github.com/Katsevich-Lab/sceptre/discussions/categories/feature-requests).
For questions about `sceptre` functionality, documentation, or how to
apply it to your data, please start a discussion under
[Q&A](https://github.com/Katsevich-Lab/sceptre/discussions/categories/q-a).
