# Aggregation precision and bounded-workspace qualification

Original 21f9ea098b69c444a884a49bc254de711968e3b5; candidate
783a2501295acc0b3fddab7bd2a615e545e5ae37. Fresh preclean release-O2 source/build/DLL
pins are in `builds/`. Installed tests pass 613 checks in 19 blocks with no failure,
error, warning or skip; separate installed original/candidate replays match all
228 complete outputs exactly and preserve both inputs. The unchanged original
function/input fixtures and reference provenance are independently reviewed.

The complete 80-child resource experiment comprises 8 exact-output children,
32 fresh single-call children,32 separate 20-call timing children and 8 allocation
profiles. All 4 complete outputs/RSS gates and all 3 meaningful timing gates pass.

| Setting | Original/candidate median 20-call seconds | Time gate | Original/candidate max measured single RSS MiB |
| --- | --- | --- | --- |
| Numeric CSC 3200 x 24000 | 21.796 /15.284 | pass | 1354.469 /622.375 |
| Numeric triplet + repeated metadata | 2.090 /2.011 | pass | 326.891 /300.172 |
| Logical CSC + stored FALSE | 1.952 /1.917 | pass | 296.219 /286.594 |
| Actual Papalexi counts; synthetic vector groups | .665 /.540 | inconclusive under 1s | 296.000 /288.062 |

Use the exact full JSON as numerical authority; the table rounds medians.
All inputs/dependency/thread/RNG settings and build hashes match; safety stops,
return failures and warning lists are empty. Both DLL hashes remain unchanged.
OS whole-child RSS is primary; active samples are lower bounds. Rprofmem shows
no full-guide-assay dense allocation. Existing dense output/transpose remain.
The measured groups contain few guides. The larger-guide column-block branch is
exact-tested across its 2^20-element boundary, not a separately claimed resource
speed result. One whole selected column is the minimum workspace and can exceed
the nominal 8MiB double-element target. No general memory-scaling claim follows.

`passing-experiment/summary.json.gz` preserves every warmup, primary/batch trial,
phase/sample/environment and profile record. Its raw-file manifest also binds
large inputs/outputs left unchanged at the recorded `/private/tmp` paths.
All 8 raw profile traces are retained. `transport-manifest.json` binds every
committed file and each uncompressed original; compression is verified lossless.

`failed-9631ef2/` preserves the previous repaired-source resource attempt: the
large setting passes, but the medium triplet timing crosses 15%, so the declared
harness stops before logical/actual settings. No exception is accepted. The
subsequent small-group source fast path removes redundant output/subsetting
work without changing the dense accumulator. The original 88ed37c Mac-only
study remains separately preserved and its cross-platform acceptance superseded.

Preparation/runner evidence retains the initial source-only test using an
unmatched user-installed support namespace (12 public-name failures; source
helper comparisons passed), corrected declared-support-namespace runs, and the
initial archive-preparation script's unsupported Python tar `filter` argument
(corrected before any R compilation). None establishes candidate acceptance.
The passing installed 613/228 evidence uses the actual 783 candidate namespace.

Independent source/resource reports and reconciliation JSON are in `reviews/`.
Hosted 3d9fafc numerical tests subsequently pass in all five lanes, while standard
R CMD check finds the SHA256 tests' undeclared `digest` dependency. The one-line
Suggests correction is separate; failed hosted runs remain failures. Final
combined installed/hosted/consumer/workflow-resource and whole-source gates remain.
