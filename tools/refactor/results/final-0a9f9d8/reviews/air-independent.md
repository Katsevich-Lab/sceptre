# Independent bounded Air formatting review

Signed source judgment: **no P0, P1 or material P2 finding identified**. The exact unit is acceptable within its authorized formatting scope. This is not combined refactor acceptance, runtime equivalence proof, performance qualification or a BiocCheck result.

Reviewer: `/root/whole_package_review`, independent of the implementer. Reviewed immutable worker commit `c05f4b90e7c51c4320903030274f80d7ebffbaec` against base `2ab3f5a2cf3d9f19242cbe8c419ada913aa721a8`. Production remained read-only. The companion `coverage.json` binds every audited R byte and all four non-R changed paths to Git blobs and SHA256 hashes.

## Exact scope and configuration

Git confirms 31 changed files: 27 handwritten files under `R/`, plus `air.toml`, `.github/workflows/check-standard.yaml`, `CONTRIBUTING.md` and the narrow authorization append to `tools/refactor/cleanup/decisions.md`. No native implementation/header, generated R/C++ binding, NAMESPACE, DESCRIPTION, test, frozen fixture, generated manual or other tool/evidence source changed. All 33 R files (32 handwritten plus the unchanged generated binding) were included in the parse/source comparison. The existing 56-slot S4 source and exported interface declarations remain unchanged.

`air.toml:4-10` retains four-space indentation and preferred width 80, explicitly enables default exclusions, and adds rooted exclusions for frozen fixture and tooling directories. The new Ubuntu job at `.github/workflows/check-standard.yaml:12-21` installs Air **0.9.0** and checks only `R` with `air format --check R`. It adds no package build or mutation step and leaves all prior jobs intact. The contributor guidance accurately distinguishes Air from BiocCheck and explains the scoped commands and retained long literals/URLs/examples. Configuration behavior agrees with the official [Air configuration documentation](https://posit-dev.github.io/air/configuration.html), including the generated `RcppExports.R` default exclusion; setup/version syntax agrees with [setup-air's documented specific-version input](https://github.com/posit-dev/setup-air#install-a-specific-version).

I independently ran the existing `/Users/ekatsevi/.local/bin/air` read-only check, confirmed `air 0.9.0`, and obtained exit 0 with no formatter differences. The implementer's disposable exclusion smoke test proves invalid generated/tooling/fixture sentinels are ignored while an unformatted handwritten positive control fails. Its byte-idempotence and YAML readback artifacts were read, not rerun by me; the YAML artifact records existing-job equality and only the new pinned formatting job.

## Expression, literal, default and brace audit

I independently wrote and ran `/private/tmp/sceptre-air-style-independent/audit.R` using R parse only. Baseline copies were first independently verified byte-for-byte against immutable base Git source. Final committed R bytes were verified against the exact after-bytes audited.

Results:

- All 33 files parse; 26 have raw-identical expressions with source retention disabled.
- The only remaining tree differences are **20 newly added single-expression braces in seven files**: `analysis_context.R`, `association_execution.R`, `execution.R`, `fit_cache.R`, `import_functs.R`, `preprocessing_functions.R` and `storage.R`.
- The pairwise comparator accepts insertion only at function-body or if/for/while/repeat body positions. It accepts no brace removal, generic block normalization, literal replacement, changed call/argument, formals/default change or quoted-language normalization. `quote`, `substitute`, `expression`, `bquote`, `alist` and formula subtrees must remain raw-identical. Function formals/default expressions must remain raw-identical.
- Seven negative controls verify rejection of changed numeric defaults, changed strings, quoted if/function bodies, braces outside permitted bodies, brace removal and quoted language within defaults.

I inspected all 20 recorded positions and their retained module contexts. They are ordinary executable branch/function bodies: coordinate/view selection; progress and maintenance; cache invalidation; sparse-column return; calibration validation; and the ondisc NT localization callback. No added block occurs in a quoted/language-data/formula/default expression or an expression-capturing analysis argument. The anonymous `lapply` callback and result-binding helper remain ordinary evaluated closures. These changes preserve ordinary evaluation, values, visibility, branch selection and order; no source/body introspection equivalence is asserted, because added braces and formatting necessarily change source/body representation.

The independent proof is stronger than erasing all single-expression blocks on both sides. It permits only the actual added blocks and keeps every other parsed component identical. The implementer's separate `strict-audit.R` independently agrees on the same 20 positions and additionally guards common quosure/expr heads. Its literal/formals audit agrees for all 33 files.

## Comments, documentation and generated help

I independently compared full-line comment blocks in all 33 files: block boundaries, word sequence and order remain identical after whitespace normalization. The retained R-parser comment-token audit additionally covers inline comments; all streams match. Roxygen text separately matches in all 33 files. The 267 retained individual prose-wrap records preserve text and introduce no list/heading/quote/fence/tag-start continuation. Code examples, executable diagnostics, URLs and mathematical expressions were not rewritten to satisfy the width preference. Existing generated manuals are byte-unchanged, and inspected public documentation wraps retain their meaning.

The residual long-line inventory is honest: 100 executable/literal lines, 17 URLs and one indivisible parameter-name list remain. This does not claim a hard 80-character limit or that BiocCheck style notes disappear. I found no unrelated style or inline `paste0` cleanup.

## Limits and acceptance boundary

I ran only the independent parse audit and pinned Air check; no package tests, oracle reconstruction, consumer/Nextflow work, runtime benchmark, native compiler, full suite or resource matrix. No new native source review is necessary for this formatting-only delta, and the earlier whole-package/native reviews remain preserved rather than replaced.

The base `2ab3f5a` resource matrix reportedly completed exact scientific comparisons but retained material CPU/RSS crossings (Papalexi cold serial memory, sparse cold serial time, CRT cold fork time/tree RSS, backed warm serial memory). Those results remain attached to their measured pin. Air formatting supplies no performance remedy or resource acceptance. The final combined immutable pin still requires the owner's required installed/package/native/consumer/resource/platform/CI/BiocCheck qualification and disposition of material failures.

### Bounded resource evidence carry judgment

I find **no material source-based inference gap** in the owner's proposed bounded carry, provided the final combined pin contains this exact formatting unit and no additional unreviewed implementation changes. This is a source/scope judgment, not new empirical resource qualification. The companion `resource-carry-source-proof.json` verifies exact Git-entry/blob equality across native sources/headers, generated bindings, tests, tool/harness/reference sources, interfaces and generated manuals; the sole permitted tool-tree delta is the non-executable authorization append. Formatted handwritten R is explicitly **not byte-identical**: its separately audited equivalence is limited to unchanged parsed expressions/literals/defaults apart from the 20 approved executable body blocks, plus unchanged comment/help meaning.

The proposed qualification is: retain the complete **160-process matrix measured at 2ab3f5a**, retaining every crossing; at the final combined pin, run **48 fresh counterbalanced processes** covering all four crossed settings and the two additional meaningful Papalexi timing settings (cold fork and warm serial), with the unchanged original baseline. Carry only the previously passing settings untouched by an implementation change. The owner additionally plans full installed/native/package/consumer/actual BiocCheck and six CI checks at the combined final pin. These are plans at this review's signing, not reviewer-executed or completed results.

That scope exercises every unresolved measured failure at final source and the actual Papalexi cost remedy while carrying only measurements across the bounded formatting-only delta. It is adequate as explicitly qualified mixed-pin evidence; it must **not** be described as a 160-process matrix at the Air/final pin, or as byte identity for formatted R. Fresh-process variability, source/introspection metadata and unmeasured final-pin settings remain limits. The source proof cannot predict timings/RSS or turn a retained crossing into a pass; actual targeted results must meet the owner's original gates or receive an explicit disposition. If the final pin adds another implementation change, this carry judgment must be reassessed for that change and affected settings.

No production change is requested by this review. No new nonblocking cleanup is proposed in this bounded unit; prior deferred P3 items remain follow-ups.

Signed on 2026-10-16 by `/root/whole_package_review`. Source judgment is limited to `c05f4b90e7c51c4320903030274f80d7ebffbaec` relative to the stated base, with report/coverage hashes supplied at handoff.
