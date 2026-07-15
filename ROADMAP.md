# sceptre roadmap

Last reviewed: July 2026

## Purpose

`sceptre` is an open-source R package for statistically rigorous, scalable analysis of single-cell CRISPR screens. This roadmap describes the maintainers' development priorities for the next two years. It is a living document: the maintainers review it at least annually and update it as priorities and capacity evolve.

## Overview of development priorities

Since we initially developed `sceptre`, its usage has increased, Perturb-seq datasets have grown in size, and predictive modeling has become an increasingly common Perturb-seq application. Our main development priorities are to mature and extend the package to support these evolving needs. Among others, these priorities include improving the package's distribution, scalability, and extensibility.

## Specific development priorities

### 1. Bioconductor submission

**Goal.** We plan to submit `sceptre` to Bioconductor.

**Rationale.**  `sceptre` has been available for download only via GitHub since its inception, rather than being included in any software repository. We expect the review process to mature our software, and the inclusion in the repository to increase its visibility and enrich its community.

**Rough timing.** We plan to submit `sceptre` to Bioconductor in time for the Fall 2026 release.

**Status.** In progress.

### 2. Improved gRNA assignment

**Goal.** We plan to overhaul the mixture-based gRNA assignment module to improve its statistical and computational performance.

**Rationale.** We discovered a [bug](https://github.com/Katsevich-Lab/sceptre/issues/220) in our gRNA assignment mixture model (in the function [`run_reduced_em_algo_cpp()`](https://github.com/Katsevich-Lab/sceptre/blob/main/src/mixture_functs.cpp#L45)). This bugged version still performs well in our benchmarks (perhaps the reason why it took so long to discover the issue), better in fact than the correct mixture implementation. Furthermore, it was [demonstrated recently](https://www.biorxiv.org/content/10.64898/2026.01.22.701179v2) that statistically accurate gRNA assignment can be attained at much lower computational cost than traditional mixture models. Both of these factors motivated us to overhaul our mixture-based gRNA assignment module.

**Rough timing.** Summer-Fall 2026

**Status.** In progress.

### 3. Hardening of the Nextflow pipeline

**Goal.** We plan to harden the companion `sceptre` [Nextflow pipeline](https://github.com/timothy-barry/sceptre-pipeline) into a more mature execution path for large Perturb-seq studies on distributed computing infrastructure. We plan to harden it by adding features including containerization, testing and continuous integration, array job support, and feature parity with the `sceptre` R package.

**Rationale.** Modern Perturb-seq datasets can exceed the memory and compute capacity of a single worker. The existing Nextflow pipeline establishes the feasibility of distributed `sceptre` execution, but is still a prototype.

**Rough timing.** 2026-2027, as capacity and funding permit.

**Status.** Planned.

### 4. Augmented output with per-pair power

**Goal.** We plan to augment `sceptre` association results with an estimate of power and/or minimum detectable effect for each perturbation-gene pair.

**Rationale.** A non-significant result can reflect either a well-powered lack of association or an underpowered analysis. Distinguishing these cases is important for conventional Perturb-seq interpretation and for downstream workflows that use association results as training or evaluation labels. This capability will build on related power-analysis methods developed in [PerturbPlan](https://katsevich-lab.github.io/perturbplan/).

**Rough timing.** 2027, as capacity and funding permit.

**Status.** Planned.

### 5. Maintainable and extensible software architecture

**Goal.** We plan to refactor the `sceptre` code base around clearer, modular interfaces.

**Rationale.** As `sceptre` has grown, its internals have become more complex and more difficult to test, maintain, and extend. Clear interfaces for data access, statistical test execution, and result generation will make it easier to add capabilities such as power-aware outputs and distributed execution without coupling them tightly to the current implementation. The refactor will also make the code base easier for community contributors to understand and support stable reference outputs for downstream implementations in other languages or on other hardware.

**Rough timing.** 2027-2028, as capacity and funding permit.

**Status.** Planned.
