# Final SCEPTRE refactor acceptance — October 16, 2026

**Complete, with two explicitly accepted bounded memory costs.** Production
source is frozen at `0a9f9d879db6a52e3daef237daa1711b9d542976`; the original
authority is `21f9ea098b69c444a884a49bc254de711968e3b5`. Finishing changes
contain developer documentation and evidence only. Nothing is merged or
released; [draft PR 234](https://github.com/Katsevich-Lab/sceptre/pull/234)
remains the review surface.

## Result

Public exports, methods, formals/defaults and the complete 56-slot saved-object
schema are unchanged. The frozen suite preserves supported numerical results,
cell/result order, RNG, cache validity, saved objects and sparse/backed workflows.
Separately reviewed corrections address named safety, ownership, import and
plotting defects in the [behavior ledger](cleanup/behavior-dispositions.md).
This is not universal equivalence on invalid or previously broken inputs.

The implementation now has concrete owners for source/analysis/fit coordinates,
storage, eligibility, fitting, score preparation, cache invalidation/commit,
separate permutation/CRT plans, native inference and result finalization.
Actual callers use those owners; obsolete routes are removed. The
[architecture guide](../../ARCHITECTURE.md) explains current change paths; the
[maintainability record](MAINTAINABILITY.md) reports benefits and retained costs.
No extension registry, new method, backend capability or saved-schema replacement
was introduced.

## Qualification and source pins

The [evidence directory](results/final-0a9f9d8/README.md) supplies reports,
content-addressed raw records, exact manifests and verified restoration tools.
Earlier trials and failures remain tied to their actual pins. Resource evidence
is a **complete 160-process study at 2ab3f5a plus a 48-process targeted study at
final 0a9f9d8**, with independently reviewed Air source carry. It is not a
complete 160-process study measured at 0a9f9d8.

| Gate | Result and scope |
| --- | --- |
| Build | Fresh immutable Git archive, preclean O2 compilation; source blobs, compiler/SDK/flags, dependency environment and DLL hashes verified. R 4.6.0, Apple arm64; private source ondisc 1.3.5, global installation unchanged. |
| Original authority | Verified unchanged original-source build and independent same-environment oracle reused. Eight original output hashes remain exact after explicit check-environment alignment; no oracle or scientific expectation replaced to clear a mismatch. |
| Installed suite at 0a | **4,740/4,740 expectations in 154 blocks**, zero failures/errors/warnings/skips; backed dependency required. |
| Native at 0a | **112 component expectations** across 17 frozen domains; 9 sampler, 214 count-ownership, 107 exact dispersion, 14 paired dispersion-failure and 959 result-ownership controlled cases. Declared negative controls retained. |
| Full package check | Final candidate and verified unchanged original: **Status OK, zero errors/warnings/notes** under `R CMD check --no-manual`, including examples and vignette rebuilding. Candidate tarball SHA256 `e442622f397248ba40788ef088b55c393d616a07602044c8cf0c82cbc40e5413`. |
| Actual pipeline scripts | Final candidate **30/30 stages**, no captured warnings; ordered payload exactly matches verified unchanged original outputs. |
| Actual Nextflow | Final standard and nuclear workflows complete, 34 R tasks each; exact paired payload/pod plan under the recorded local serialized `fair=true` shared-arrival/seed configuration. |
| Selected book cells | **14/15**, same original empty-barcode example failure and exit 1 retained. Completed payloads exactly equal. |
| Resources | All meaningful timing settings pass under source carry and targeted repeat. Papalexi cold serial RSS remains a crossing with Eugene's separate bounded acceptance; details below. |
| Maintainability | Four fixed traces, every narrow site match and seven fixed operations reviewed. **7/7 facade-independent implementations, 6/7 ordinary-input production R entries**; mixed complexity and code-surface costs disclosed. |
| Independent source/guide review | Whole 32b2076 source and architecture review plus separately reviewed native D-access and Air deltas; no unresolved P0/P1/material P2 source finding. Exact source/AST/consumer carry retained, rather than another whole-source review at each pin. |
| Actual BiocCheck at 0a | Original **0 errors/1 warning/5 notes**, final **0/1/6**; matched BiocCheck 1.48.1/R 4.6.0/Bioconductor 3.23 (BiocVersion 3.23.1), no disabled checks or introduced errors/warnings. |
| Hosted CI at 0a | **All six jobs pass** in [push run 38954047191](https://github.com/Katsevich-Lab/sceptre/actions/runs/38954047191): Windows/macOS/Linux release, Linux oldrel-1, no-ondisc and Air. Raw metadata/logs retained. |

Final consumer execution compares new 0a outputs to unchanged, source/DLL/input/
dependency-verified original outputs from the 32b2076 qualification. The
original was not rerun during this last candidate-only consumer refresh.
Pipeline source is `38dc760a5b4105d1c3f339cd45f2c98d971087f4`; book source is
`63704b0cca4b334f5b957da928ce718d38184f20`.

## Resources and explicit dispositions

Thresholds remain at most 10% elapsed increase for reference medians of at least
one second, and at most the greater of 10% or 64 MiB sampled RSS increase.
Preparation is outside the discovery timer; cold means an empty fit cache, not
a purged OS filesystem cache. Each setting uses one warm-up and three measured
fresh processes per role, counterbalanced within settings. Scientific/input/B/
pair/configuration/warning checks remain exact. Profiling is separate.

The complete 2ab matrix covers five cases, cold/warm and one/two workers. The
final 0a repeat covers every crossing plus the other two meaningful Papalexi
timing settings. Independent source carry establishes the Air delta's literal/
formal and narrowly normalized executable-AST preservation; 864 native/binding/
test/harness/manual/interface Git entries are unchanged. Air source proof is
not itself a performance measurement.

| Actual 0a targeted setting | Reference/candidate median seconds | Main RSS increase | Outcome |
| --- | --- | --- | --- |
| Papalexi cold, 1 worker | 2.645 / 2.517 | +162.30 MiB | Time passes; memory crossing accepted separately |
| Difficult sparse cold, 1 worker | 2.361 / 2.310 | +10.45 MiB | All gates pass |
| Shared CRT cold, 2 workers | 1.299 / 1.291 | -0.88 MiB | All gates pass; worker/tree RSS also decrease |
| Backed warm, 1 worker | 0.478 / 0.419 | -26.62 MiB | Memory passes; timing diagnostic below one second |
| Papalexi cold, 2 workers | 1.926 / 1.541 | +36.25 MiB | All gates pass; worker/tree RSS decrease |
| Papalexi warm, 1 worker | 1.203 / 1.134 | +10.39 MiB | All gates pass |

The other sparse/CRT/backed crossings in the 2ab matrix did not repeat; their
exact causes remain unknown. Those observations are retained. A required native
timing failure at 32b was isolated to repeated Rcpp matrix access in the guarded
null-score loop. Caching the existing D data pointer at 2ab restores performance
without changing guards, iteration/accumulation order, draw banks or results;
native replay, actual-source tests and independent review are retained.

**Papalexi whole-workflow memory exception:** final cold serial sampled main/tree
maxima are **504.96875 versus 667.265625 MiB**, an increase of **162.296875 MiB**.
The complete 2ab batch previously showed about **83.5 MiB**. Eugene explicitly
selected **“Accept this bounded cost and finish”** after seeing both results
and the unresolved causal limit. The generic gate remains crossed, with a
separate accepted exception in [D2026-10-16-04](cleanup/decisions.md). This does
not extend the earlier conversion acceptance or establish an unavoidable minimum
or larger-input bound. Preparation/heap-history memory investigation is follow-up;
restoring caller mutation is unacceptable.

The final isolated ownership diagnostic has eight fresh unprofiled controls
and two separate preparation profiles. Both controls load/byte-compile the same
storage function; the diagnostic variant deliberately restores unsafe borrowed
numeric-triplet sorting. Science/configuration hashes match. Production-minus-
borrowed paired RSS differences are **-15.30, +2.55, +6.14 MiB**; production and
borrowed measured maxima are **570.52 and 582.80 MiB**. This does **not** attribute
the original-versus-final 162 MiB gap to copying. Separate R allocation traffic
is 123,225,368 versus 73,534,952 bytes, a **47.39 MiB** increase; all 56 retained
slot sizes match and pre-call live vector use is essentially unchanged.
Large allocations appear under anonymous storage frames: the lexical
`copy_stack` category's equal 35,736-byte totals do not measure ownership copies.
Rprofmem is allocation traffic, not retained/peak memory, and excludes native/
HDF5 allocations. The unsafe variant is not a proposed production change.

Separate historical aggregation/conversion studies retain their actual pins.
B6 uses complete selected rows in bounded dense column blocks, preserving the
platform's original dense accumulator, order and repeated metadata multiplicity;
exact-output/RSS and meaningful timing gates pass. Dense output/transpose
workspace remains; one selected column can exceed the nominal 8 MiB block target.
The numeric eight-million-entry owned-triplet cases retain accepted
**+121.938/+123.141 MiB** whole-child costs under D2026-10-16-02. Logical/small
and separate dense-copy gates pass. Air changed R text: carry uses retained
syntax/source proof, not false final R-byte identity against 783a250 or 909aa3a.

The 16 earlier separate fit/bank/allocation diagnostics retain exact expected
reuse: shared permutation response fits and one full-budget bank, CRT response
prepass/guide traversal and whole-guide banks, warm cache behavior and complete
staged budgets. They are pinned historical diagnostics, not extra final timing
trials. Invalidated attempts remain excluded, including the first 0a partial run
started despite a foreign R process in preflight. Only this chat's own benchmark
tree was stopped; the foreign job ended independently. The clean restart receipt
records no competing R/compiler/compression jobs.

## Maintainability, hygiene and limits

The fixed-denominator result improves from 5/7 to 7/7 facade-free implementations
and 4/7 to 6/7 ordinary-R entries. Staged inference is a private C++ component
with Rcpp arrays and a const draw bank; the installed Rcpp entry still uses an
external pointer. The ordinary-list bridge is test-only.

R exclusive maximum decreases 49 to 44 and inclusive named-owner p95 decreases
22.75 to 19.25; exclusive p95 rises 13.85 to 15. C++ maximum rises 17 to 19 and
p95 13.60 to 15.85 as explicit safety/ownership surface grows. Lizard nd 11 is
a formatting/operator-sensitive heuristic, not eleven literal nested blocks.
Final R physical/token-bearing **lines** are 10,193/6,276 versus original
8,662/6,086; Air wrapping changes line counts. Raw slot matching lines are 35
in two files (33 adapter, two owned getter lines), versus 51 in nine files
originally. The pre-Air 26-to-35 increase is line splitting: both have 40 regex
occurrences. No comprehension timing, cognitive complexity or aggregate quality
score is claimed.

Air 0.9.0 is the single formatter, with four spaces and preferred width 80.
Independent strict review accepts only 20 added single-expression body braces
across seven files; other executable ASTs, literals/formals and help meaning are
preserved. Comment wrapping is separate; generated bindings and immutable
fixtures/tooling are excluded. BiocCheck's indentation note clears; 146 long
lines remain versus 145 original, including intentional code/literals, URLs and
generated bindings. One extra nested filename-pattern `paste0` note in
`R/execution.R` is a nonblocking hygiene follow-up. Existing no-Bioc-dependency
warning, administrative checks and inherited notes remain. This is not
Bioconductor submission approval.

Native probes certify declared buffer/payload boundaries, not global Rcpp
precious-list cleanup, every R error path or real OS OOM. Sampled RSS can miss
peaks; aggregate fork RSS double-counts shared pages. Cases include actual
Papalexi counts and synthetic stress inputs, not full Replogle/out-of-RAM scale.
The backed dense counterfactual exceeds a declared 256 MiB materialization
budget, not an OS-enforced memory cap.

Nextflow qualifies the disclosed local arrival-controlled configuration;
uncontrolled arrival, historical containers, arrays, nf-test/hosted schedulers
and whole-book rendering remain outside the evidence. Repeated file coordinates
are rejected by SCEPTRE's owned Cell Ranger/Parse constructors under Eugene's
choice; delegated ondisc parsing is unchanged and unqualified for that policy.
Valid triplet duplicates and repeated guide metadata are separate contracts.

Backed CRT/aggregation, SPA, Wilcoxon, accelerated fitting, replacement pairwise
QC, Bioconductor data classes, new outcomes, PR233's guide cache and new Python
execution formats remain future work. Broader malformed-input precedence,
numerical fallback changes, Bonferroni plot interpretation and global plot
sampling remain deferred. Nonblocking source/tooling residue retains explicit
dispositions; completion does not mean absence of every possible defect.
