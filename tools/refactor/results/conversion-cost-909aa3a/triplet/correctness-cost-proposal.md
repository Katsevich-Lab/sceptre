# Optimized triplet ownership cost disposition proposal

**Status: proposed measured exception for independent review; the declared gates remain unchanged and their crossings remain recorded. No exception is accepted by the measurement report.**

Reference 21f9ea098b69c444a884a49bc254de711968e3b5; candidate 909aa3ad76ee21e7f9ab83a63c645bc0fb6a43c7 (sparse ownership production from reviewed da6a079, dense repair from reviewed bc53bc1). Fresh separate release-O2 installs have matching compiler prefixes. Source/build/DLL/harness/fixture hashes, all warmups, four alternating pairs, separate profiles and required five-conversion batches are bound in the complete summary. All original 8841efa trials remain preserved.

Measured report: `/private/tmp/sceptre-triplet-qualified-909aa3a/summary.json` (SHA256 `9dfd92d6cf625b8889da8931f5f811488fef32b7a8a86cda98c9564db7db2868`).

## Measured gates

|nnz/type|view|OS original MiB|OS current MiB|OS delta MiB|sampled phase delta MiB|gate|
|---|---|---:|---:|---:|---:|---|
|1000000-numeric|row|258.562|274.000|15.438|-3.781|pass|
|1000000-numeric|column|268.578|283.656|15.078|0.125|pass|
|1000000-logical|row|270.312|272.734|2.422|-1.562|pass|
|1000000-logical|column|265.688|269.953|4.266|10.672|pass|
|8000000-numeric|row|495.594|617.531|121.938|122.203|CROSSING|
|8000000-numeric|column|503.484|626.625|123.141|123.125|CROSSING|
|8000000-logical|row|525.875|587.766|61.891|61.062|pass|
|8000000-logical|column|528.781|586.969|58.188|58.203|pass|

All recorded sampled phase values are lower bounds; OS whole-child peak is authoritative. Three of 112 warmup/timing children lack an in-conversion sample, all in the short 1M cases; all 8M children have conversion observations. Missing observations are counted per case/build in the separate final qualification summary. They do not become fabricated phase passes. No real or sampled process exceeded the 2 GiB safety stop.

## Allocation evidence and improvement

|nnz/type|view|original traffic MiB|current traffic MiB|delta bytes/nnz|8841efa current traffic MiB|
|---|---|---:|---:|---:|---:|
|1000000-numeric|row|13.005|28.264|16.000144|28.413|
|1000000-numeric|column|13.120|28.379|16.000144|28.528|
|1000000-logical|row|20.635|28.264|8.000096|36.043|
|1000000-logical|column|20.749|28.379|8.000096|36.157|
|8000000-numeric|row|99.790|221.860|16.000018|222.009|
|8000000-numeric|column|99.904|221.974|16.000018|222.123|
|8000000-logical|row|160.825|221.860|8.000012|283.044|
|8000000-logical|column|160.939|221.974|8.000012|283.159|

These are NET original-to-current allocation traffic differences, not live memory peaks or isolated sorting workspace. The profiler records one newly owned numeric i,j,x vector at 1M/8M positions (four/four/eight bytes per position plus R headers); for logical inputs x is already freshly replaced by one-values, so only the borrowed i,j columns require duplication. The removed logical x copy is directly visible as one fewer 64,000,048-byte allocation at 8M. Numeric-valued input bits, including signed zero and custom NaN payloads, are copied without arithmetic.

## Why the copies exist and alternatives

The original by-reference sort changes supplied sparse slots, as reproduced on all eight triplet correctness cases. The current conversion must own every borrowed column it sorts to preserve caller aliases. It reuses only constructor-proven fresh logical values, expanded compressed coordinates and dense-coercion columns; no NAMED/refcount heuristic is used. It retains stable coordinate order, repeated stored positions without coalescing, explicit zeros, fractional values, class/dimnames and logical stored-FALSE presence semantics.

- Retaining the original borrowed in-place sort would reintroduce demonstrated caller mutation; rejected.
- Skipping copies based on refcounts or using arithmetic to manufacture copies would weaken alias/bit contracts; rejected.
- Recopying constructor-proven fresh columns was avoidable and has been removed, measured and independently reviewed.
- Public Matrix conversion can coalesce repeated triplet coordinates and therefore cannot silently replace this existing stored-entry contract.
- A new native stable sparse-construction/counting-sort implementation could trade coordinates, permutation storage and passes differently. The measured 16 bytes/position is NOT a mathematical lower bound across all possible algorithms. Such a replacement needs its own duplicate/stability/bit/nonmutation proof and matched resource experiment; no unmeasured claim that it would be faster or smaller is made here.

## Proposed bounded disposition

Accept only the recorded 8M NUMERIC triplet row/column OS and sampled-phase increase crossings for this exact reviewed source and declared fixture family/platform. Proposed empirical ceilings: additional whole-child OS peak <= 128 MiB above the matched original and current whole-child OS peak <= 650 MiB for each declared 8M numeric setting. The same 128 MiB delta ceiling applies to the observed phase lower-bound comparisons, without upgrading them into exact peak measurements. All ordinary 1M memory gates and all logical memory gates retain their original thresholds; all meaningful batch timing gates retain <= 10% increase.

These ceilings are specific to this completed experiment and allow the observed ownership cost while preserving a large margin below the independently enforced2GiB process safety stop. They are not a 16 bytes/nnz automatic allowance, a waiver for larger or different matrices, or acceptance of arbitrary future code changes. The two-size profile supports approximately linear copy traffic only over the measured 1M-to-8M range. This does not establish peak memory at larger sizes.

Final integrated workflow/backed-memory/worker/caching/consumer gates remain required. Scientific fixtures and numeric tolerances remain unchanged. Dense conversion is separately characterized, its near-symmetric correctness repair is independently tested, and no sparse/backed path densifies.

## Implementation completion limitation

Source correctness, nonmutation, native exceptional ownership, exact frozen workflow/API/index/native contracts and the full bounded experiment are complete on the qualified source. Resource acceptance is still open until the independent reviewer/root explicitly adopts or rejects this measured disposition. All crossings, outliers and earlier trials remain intact.
