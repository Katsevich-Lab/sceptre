# Native result lifetime repair (B8)

Production checkpoint: **11502f871ebeb78d7cf084c45ffdb895d74ac937**, on the separately preserved CL4 checkpoint dc093d27acefd2b3897eaf8b5a0bb4e9cbe7282f.
Final reusable probe checkpoint: **6dbd1d0621fd88716423f8d77cacf176dce66fd3**; production files are identical to 11502f8.
The [readback](native-result-lifetime-evidence.json) records exact SHA256/MD5 pins, before/after checks, transformations, dependencies/threads, observed boundaries and all workload totals.
Independent signed review and the final combined-package gate remain separate. This does not claim universal Rcpp/OOM safety or completed whole-refactor acceptance.

## Known defect and bounded repair

Raw R allocation errors can bypass C++ destructors. Native calibration sets/results and TailResult::null_statistics previously remained live during R output allocation/wrapping.

- Groups allocate the exact integer matrix and materialize its pointer before owning the set/sample. The same seed, rejection loop, sorted membership and unordered-set row visitation execute within a native-only scope.
- Pairs materialize the input pointers and protect a continuation before native scratch. The original candidate permutation, draws, eligibility and append order populate only accepted-result vectors. The complete native scope ends before a trivial raw pointer transfers into protected conversion.
- The association entry checks shapes before native ownership and finishes owned observed-index/XPtr temporaries before allocating its tail payload. Staging, arithmetic, full-bank borrowing, segment order and fallback remain unchanged. The duplicate shape check in null_score_segment was removed: D.ncol() itself uses Rcpp protection. Every real owning caller now checks shapes before native ownership.
- Pair/tail conversion uses public R_UnwindProtect with a previously protected non-NULL token, raw-C-API-only leaves and allocation-free noexcept cleanup. Cleanup deletes the payload once on success or before the original R jump resumes. No C++ exception crosses an R protection frame. The tiny continuation guard unwinds normally on a C++ compute error, and its entire scope ends before protected conversion.
- Output retains actual types/fields/names and exact accepted/selected length. There are no request-sized rejected-pair output buffers or unreached B2/B3 allocations.
- Three raw matrix row strides widen to R_xlen_t before multiplication, preserving former Rcpp matrix offsets for valid long matrices.
- The queued D contract is corrected: rows are eigenbasis coordinates in eigen() order, dimension equals design coefficient count, and columns follow fit cells.

R public exports/formals, S4 schema, R wrappers, frozen fixtures, statistical decisions, budgets, RNG order and file effects are unchanged. Only the private generated C++ calibration-pair return declaration becomes SEXP; its R return stays the same named list. Both bindings were generated together; the generated R wrapper is byte-unchanged.

## Source/API support and review

The active senior independently approved this design after inspecting supported minimum R 4.1: a supplied token skips allocation before begincontext/SETJMP; cleanup precedes R_ContinueUnwind. Leaves have no Rcpp precious/native C++ owners, native allocation or C++ throws. Guard transfer changes only the protection obligation; no broad exception catch is introduced.

Primary references: [R extensions manual](https://stat.ethz.ch/R-manual/R-devel/doc/manual/R-exts.html#Condition-handling-and-cleanup-code), [R 4.1 context implementation](https://github.com/wch/r-source/blob/R-4-1-branch/src/main/context.c), [public declarations](https://github.com/wch/r-source/blob/R-4-1-branch/src/include/Rinternals.h).
Downloaded source hashes remain at /private/tmp/sceptre-native-output-design-sources/source-manifest.json; context.c SHA256 is f76ddb376b6c99896d97bfa9901c2d3fc739ce5e4c4eb64ac3caefcd45a39015.
No R/Rcpp minimum is raised. Local execution is R 4.6/Rcpp 1.1.1.1.1; source support verification is not execution on minimum versions.

| Source/functions | Inspected responsibility |
| --- | --- |
| result_lifetime.h: ResultContinuation constructor/destructor/release | Pre-owner protection, C++ failure cleanup, nonallocating transfer, final explicit unprotect |
| negative_control_functions.cpp: sample_calibration_groups | Exact output first, unchanged subset stream/visitation, wide column-major offsets |
| same: sample_calibration_pairs, calibration_pair_record, delete_calibration_pair_result | Cached R inputs; accepted vectors; unchanged candidate/count policy; exact output and normal/error cleanup |
| low_level_full_test.cpp: run_association_test | Pre-owner validation, observed conversion scope, borrowed bank, ended native ownership scope |
| inference_stages.cpp: association_record_leaf, delete_association_tail, association_result_record | Exact record/selected segment; trivial context; raw conversion and exactly-once cleanup |
| inference_scores.cpp: null_score_segment | Aligned inputs required before native computation; no duplicate R shape protection |
| inference.h and entry/score comments | Ownership/protection obligations and precise score/eigenbasis coordinates |
| component bridge: component_segment, component_tail_record | Direct callers check before ownership and use the consuming record interface |

No unrelated native cleanup is part of this repair.

## Bounded execution evidence

The final immutable archive is /private/tmp/sceptre-native-lifetime-ci-source. All R/compiler processes exited before the shared slot was released. Root owns the combined package install/suite; those were not duplicated for this unit.

The reusable command has one explicit original-source interface, without fallback:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    Rscript tools/refactor/test-result-ownership.R \
      IMMUTABLE_CANDIDATE_SOURCE VERIFIED_COMPILED_ORIGINAL_SOURCE FRESH_EVIDENCE_DIR

It verifies the original's complete pinned 21f9ea0 SHA256 manifest/inventory, loads that source via pkgload::load_all(compile=FALSE, helpers=FALSE, export_all=FALSE), and checks loaded DLL bytes equal its prebuilt DLL. It does not install/rebuild the original or select candidate expectations. Local original: /private/tmp/sceptre-b6-precision-9631ef2-attempt2/reference/source.

The probe compiles actual candidate source after mechanical test-only allocator/lifetime instrumentation, preserving owner scopes/algorithms. It passes **959 cases: 12 successes, 858 controlled raw R errors, 89 native bad_alloc points**, plus **one separately counted safe negative control**. Every main case ends with zero live tracked native blocks, payloads and contexts, and no unsafe boundary. An error with an owned leaf live requires exactly one result cleanup. The negative control detects deliberately unsafe ordering through a C++ error before any raw jump and leaves zero live counters. No real OOM occurs.

Groups/pairs are identical() to independent original results, including types/order. Tail fields/types/order match unchanged frozen expectations at their existing 1e-8 tolerance. Input nonmutation and R RNG state are checked.

**Nine boundary kinds were observed:** vector allocation, INTEGER/REAL materialization, names, attributes, precious preservation, callable lookup, continuation creation and PROTECT. Four available interpositions were unobserved: allocMatrix, R_alloc, coercion and direct R_PreserveObject. The readback distinguishes availability from observation and lists exact baseline/all-visit counts.

The 12 workloads include accepted/all-rejected pairs, both accepted skew-normal tails, failed-fit/B3-zero and empirical fallback, variable CRT cardinality, compact output and stage 1 with enormous unused B3. That stage-1 case tracks **3,992 native container bytes**, identical to ordinary stage 1, without allocating unused budgets. The all-rejected 50,000-request pair case tracks **200,000 bytes**, exactly its existing permutation and no accepted-vector capacity. These allocator counts exclude payload object/header/R memory; they are not RSS/general resource acceptance.

The uninstrumented component bridge also links and passes **3 tests / 112 expectations** across all 17 frozen native domains, including segment ordering, effects/scores, compact records and RNG. Log: /private/tmp/sceptre-native-lifetime-components.log.

Final artifacts and exact hashes are in the committed readback. Main RDS: /private/tmp/sceptre-native-lifetime-ci-probe/result-ownership-results.rds, SHA256 db90738b1378558c174b10b32ac67a0a2b10c37c7ba8450eac6da26ab6dfb597. Full log: /private/tmp/sceptre-native-lifetime-ci-probe.log.

## Preserved attempts and evidence correction

1. At 30e7aef compilation succeeded, then probe accounting incorrectly counted Rcpp's separate post-unwind error-translation cleanup as result cleanup. Failed/diagnostic logs remain under /private/tmp/sceptre-native-lifetime-30e7aef-*.
2. At e65a6ed compilation and group/pair cases progressed; the safe probe detected duplicate D.ncol() protection with live tail ownership. Logs/state remain under /private/tmp/sceptre-native-lifetime-{final,debug2}*. Senior-approved pre-owner correction: 11502f8.
3. The first two passing runs (qualification and qualified) exposed evidence retention failure: final copying of original basenames overwrote transformed source after compilation. Their pre-copy transform hashes/results remained. All six exact transforms per run were reconstructed from immutable sources/original include paths; **all 12 MD5s matched their pre-copy records**. Recovered bytes are separate. [Recovery readback](native-result-transform-recovery.json).
4. Harness 6dbd1d0 separates transforms/ and originals/, refuses nonfresh evidence directories and rechecks retained transform hashes after copying. Its explicit verified-source CLI reran all 959 cases successfully. Production remains 11502f8.

Final combined build/oracle/suites, cross-platform CI, consumers and general resource gates belong to the integration checkpoint. This repair does not imply another unchanged triplet resource run.

The bounded evidence archive and per-file manifest are described in [the bundle receipt](native-result-lifetime-bundle.json). The final raw source transforms and original copies are retained separately; no compiled probe DLL persistence is claimed.
