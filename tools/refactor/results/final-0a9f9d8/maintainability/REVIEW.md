# Final maintainability refresh: 0a9f9d8

Frozen integration source **0a9f9d879db6a52e3daef237daa1711b9d542976**, Git tree de329c06874db5353e27676140b2ebe60f7e53e8. This independent engineering scorecard refresh uses the actual post-Air R/source bytes. No main/docs/production edits or package tests/builds were performed. Source/caller review and prior exact component qualification are separate from automated scores.

## Outcome and provenance

The fixed source-tracing exercises and testability denominator still support the refactor's ownership/locality improvement. Raw complexity distributions are **unchanged from32b2076**; the final result is not a claim of uniformly reduced complexity, measured human comprehension, perfect source quality or a resource pass. Native validation/ownership surface is an explicit retained cost, and formatting alone does not establish maintainability.

Used the established cyclocomp **1.1.2** library `/private/tmp/sceptre-rf0-metric-library` and Lizard **1.24.1** with bundled nd under `/private/tmp/sceptre-rf0-metric-python`; analyzer file hashes match the prior established environment. R4.6.0/Python3.9.6 are retained in raw records. No analyzer was installed or substituted. The exact frozen `measure.py`, `measure-r.R` and `measure-sites.py` bytes match the prior run. `measure.py` ran its AST-only R script **once**; the site scanner ran once. Neither evaluates the measured package.

Verified every archived R/src Git blob against0a **before** measurement and rechecked source hashes after. Per-function/per-file hashes, exact commands, environment, scripts, source.tar and records are retained. Current handwritten aggregate SHA256: `c67f2f24f45ba19289718d84a1b97d356ed12021c25b829131f917524a7423c8` (32R+16native source/header files). Generated bindings/tests are excluded from analyzer denominators.

Frozen baseline denominators and detailed measurements are reused unchanged from the prior verified21f9ea0 run: both detailed R/C++ files had reproduced byte-for-byte against RF0. The same35 baseline handwritten files were copied into `baseline-source/` and **all35 frozen hashes verified now**; aggregate baseline identity remains `d405a883a1e01c3d1c69554c5121d0b1742d6d4429d38e966b85c81d0596d2ef`. Baseline analyzer measurements were not needlessly rerun. `prior-32-proof/`, `frozen-baseline/`, `baseline-rerun/` and setup.json preserve these distinctions. One report-assembly attempt initially referenced a missing local copy of the prior sites file; it was corrected by copying that unchanged record, without any source/analyzer failure or analyzer rerun (reporting-attempts.json).

## Within-language raw measurements

| Measurement | Frozen original21f | Final0a |
|---|---:|---:|
| R exclusive function literals |184|221|
| R exclusive max / p90 / p95 |49 /9.70 /13.85|44 /9.00 /15.00|
| R named top-level owners |122|160|
| R inclusive named-owner p95 |22.75|19.25|
| R AST control-depth max |7|7|
| C++ functions, including header methods |43|64|
| C++ cyclomatic max / p90 / p95 |17 /8.80 /13.60|19 /10.70 /15.85|
| Lizard nd maximum |5|11 (heuristic)|

No sum across split functions or languages is used. R nested functions are separate exclusive units; inclusive owner values are not added to them. The permutation driver remains3/7 exclusive/inclusive versus7/15 original; CRT remains4/11 versus9/23. Each of four actual execution workers is3/3. Native entry9, staged tail8, null-score segment9 and skew-normal15 are distinct units, not a combined score. Calibration pair sampling19 and sparse pair counts18 show the added explicit safety/ownership surface; R assignment validation still has maximum44 and native-count validation29. Namespace-aware hotspot rows are retained separately because the fixed display selector misses some namespace-qualified names; raw measurements contain them.

Lizard nd11 is **not11 literal nested blocks**. Its pinned token heuristic is formatting/logical-operator sensitive; the retained prior static probe gives6 for three sequential unbraced logical guards and3 for the braced equivalent, while a true three-nested-if example gives3. The same raw heuristic is preserved, not replaced by a favorable adjusted score. Air formatted R only, so this C++ result does not change. R's documented AST depth rule is different. Neither measures cognitive complexity.

LOC is supporting only: R physical/token lines change from baseline8662/6086 to final10193/6276; prior32 was9383/5732. Air line layout and comment wrapping therefore raise line counts without a complexity change. Native physical/nloc is1211/849 baseline,1923/1259 at32,1925/1260 final: the small final increase is the validated D-pointer cache/comment. These are maintenance-surface diagnostics, not readability scores.

## Four fixed source-tracing exercises

`fixed-source-tracing.json` preserves the original four questions, current declaration/caller locations, observed paths, constraints and no resulting implementation diff. `actual-caller-search.txt` records the current search. Refreshed relevant owning bodies/callers and used the reviewed Air AST/source carry for unchanged surrounding paths. Reviewer has substantial prior package/refactor familiarity; there is no blind human experiment, elapsed comprehension time, inferred reading count or hypothetical changed-file count.

1. **Eligibility:** pair_eligibility.h owns inclusive valid-count and R missing-value compatibility. Discovery's plain vector adapter and actual calibration sampler call this owner. Measurements remain in storage/pair_qc. No second policy is hidden in a wrapper; new scientific inputs/policy would need an explicit future contract.
2. **Response fitting:** fitting.R owns Poisson/coefficient/dispersion and assignment logistic fits; score_preparation.R owns arrays/eigensystem. Actual shared-permutation and shared-CRT misses plus both NT guide-specific paths call the fitter. Context/views define aligned rows; NT discovery remains uncached by guide. Response-first shared permutation banks and guide-first CRT banks are necessary distinct traversal/reuse semantics.
3. **Storage:** assignment thresholds, pair QC, calibration group counts and CRT reads use storage operations. Native counts consume canonical analysis groups; external ondisc alone receives its localized/owned compatibility inputs. Stored zeros, repetitions/order, one-row materialization and existing all-response summary scan are explicit. Separate sparse/backed kernels remain legitimate implementations.
4. **Cache:** supported parameter/design/contrast/control-view, low-MOI assignment and retained-QC-cell changes resolve three explicit events into one invalidator. Lookup/append/commit have one owner and three analysis callers. Cache response-ID keys do not independently guarantee view validity; arbitrary manual mutation remains outside supported setters.

These are observed source ownership paths rather than an automatically computed semantic coupling score. They do not simulate a future batched fitter or backend.

## Site denominator/classification after Air

| Identical scanner | Original lines/files | Prior32 | Final0a |
|---|---:|---:|---:|
| Private ondisc |6/4|6/1|6/1|
| Raw slot candidates |51/9|26/2|35/2|
| Narrow view-remap candidates |17/4|1/1|1/1|
| Cache-reset implementation |3/1|1/1|1/1|
| Cache-writeback implementation |3/1|1/1|1/1|

Every final matching line is classified with its actual function in `site-classification.json`, and the original per-line classifications remain retained. Final raw slots comprise33 adapter lines and2 fresh assignment-matrix construction lines. **The26→35 change is line splitting:** the same narrow regex has40 raw-slot occurrences at both32 and0a, in the same two responsibility locations. Neither count is a new semantic operation or architectural regression. The sole remapping match is a legitimate plot/QC display subset. Three mutation-event and three successful commit callers remain, despite one implementation each.

Patterns are deliberately narrow matching-line diagnostics. They omit other constructor fields/map names/public routes and are not complete coordinate/coupling inventories. Private ondisc remains an explicitly selected compatibility dependency, not a removed debt or newly public storage API.

## Seven fixed testability operations and actual qualification

**7/7 facade-independent implementations,6/7 direct ordinary-input production R entries,1 private C++ component.** Baseline remains5/7 facade-free and4/7 ordinary-R. Ordinary inputs include supported Matrix/backing handles and compact fit/settings lists; not zero preparation. The fixed operations are response fit, assignment fit, score pieces, eligibility, staged inference, family result finalization and response-row read.

Staged inference is `sceptre::inference::evaluate_staged_tail`, taking Rcpp numeric arrays and a const C++ DrawBank. The installed production Rcpp `run_association_test` **still takes XPtr**. A developer-only bridge accepts ordinary R draw lists and constructs/copies the bank; it is neither an ordinary-R production interface nor a newly public API. Bridge compilation, Rcpp runtime, bank preparation and nonmutation copy remain real burdens.

This reviewer executed the fresh -O2 actual-source component qualification at workerde275681 (integration2ab):112 expectations/3tests across all17 frozen domains, including selected stages/distributions, compact outputs, RNG and input/bank nonmutation. Every relevant source/header/bridge/runner byte is verified identical from that executed source to0a; copied records and per-file hashes are in `component-qualification/` and `independent-testability.json`. No new component test ran during this metric refresh. Current direct R test definitions for fits/score/eligibility/results/row access remain independent ordinary-input tests; their source is retained. This is genuine interface/component evidence, not exhaustive numerical-state coverage or a percentage coverage claim.

## Qualification limits and follow-ups

The coordinator is running final0a installed/native/package checks. No completion or resource result from those active jobs is inferred here. The declared resource evidence plan uses the full160-process matrix at **2ab**, independently reviewed source carry, and48 actual **0a** targeted processes for crossed/Papalexi settings. Do not relabel the full intermediate matrix as measured at0a. Prior studies and observed crossings/exceptions remain separately pinned.

Nonblocking follow-ups: improve the fixed metric display's namespace matching; use a trustworthy pinned native AST analyzer if literal nesting becomes a criterion; optionally reflow prose paragraphs to avoid short orphan comment continuations introduced by conservative per-line wrapping. No production reopening or new implementation is requested. Cognitive complexity, timed human comprehension, clone coverage, complete dynamic coupling and an aggregate 'pristine' score remain unmeasured. Compatibility/static metrics alone do not prove broad source quality.

All analyzer R work is complete. After packaging/hash readback, the reviewer will send an explicit quiet acknowledgment and perform no further runtime/compiler/profile/archive/broad hashing until released.
