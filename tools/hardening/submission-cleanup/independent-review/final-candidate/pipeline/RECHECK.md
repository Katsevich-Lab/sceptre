# Pipeline compatibility-refresh review

**No new material pipeline issue was found for the exact01bb39c artifact under the bounded local scope. New image/scientific/scheduler identities agree, with both representative Docker workflows and the distributed failure guard passing.**

Reviewed pipeline: `/Users/ekatsevi/.codex/worktrees/bioc-compatibility-refresh/sceptre-pipeline`, branch`codex/bioc-compatibility-refresh`, commit`818a156a8a53cf3e34e58b05bd37d115eb4ced49`. The production baseline remains818a156. Evidence-only commits `6d35036ffab74f9e1201199f11824e90beb06044` and final `bfce9a812c45a7a9e389f1650ff2b42745c7ea89` now preserve the refresh receipts. The final evidence checkout is clean. No production/config/test/containerrecipe paths changed. A transient two-file runtime metadata omission was corrected inbfce9a. No production implementation or Git mutation was performed by this reviewer. The dirty72cc arrays/HPC checkout was not reviewed or qualified. All23 production pipeline files are unchanged from the previously reviewed implementation/guard boundary.

## Exact artifact and independent runtime readback

- New package source: `01bb39ccf8373ae8965a2679ee4813a63ea2d85d`.
- Default-build archive: `/private/tmp/sceptre-bioc-cleaned-01bb39c/sceptre_0.99.0.tar.gz`, independently hashed SHA256 `606363fce32238354fc0d126d79f515973da7c8287ca43b14ce53ce2d00df4b7`.
- New local image: `sha256:1510c16a1d5f0da6a266f8685c4b06fa1865e0116b159956af9625dd1cc249b2`.

An independent temporary container with networking disabled and no host mounts verified the embedded source manifest and actual embedded archive bytes against this exact source/SHA. Runtime: LinuxARM64 R4.5.1, SCEPTRE0.99.0, ondisc1.3.1, Rhdf5lib1.30.0 and OpenBLAS. Installed SCEPTRE SHA256 is `8f59580aaf43aee9218bb379acda39d7ce1004e8fe65e4a2012a0a7ce98dec9b`; installed ondisc SHA256 is `b5a4a7c23c8b6469b14c787b4a46e482646cb756d012f67bb1def568d8eae552`.

Actual loaded native SHA256 values are SCEPTRE`b7ea2540d26514efdb0442ac2d418554d3f139ae0f7f550bb107c22f8fb5b6cf` and ondisc`e3a65d52160de9d271853780fac71b80478c4d1e55864a3d49da4a28a6c2f53a`. Their equality to earlier native binaries was verified, not presumed. The new installed R package identity differs. Full original runtime readback is in `image-runtime-readback.json`.

## New terminal receipts independently reviewed

These suites were executed by the worker; this reviewer inspected saved logs, RDS reports, run signatures/manifests and traces without repeating broad science or builds.

| New-artifact receipt | Verified result |
| --- | --- |
| Direct R suite |9memory cases plus2required backed cases; no skip |
| Distributed R suite |9memory cases plus2required backed cases; no skip |
| Docker distributed high-MOI CRT |same-image direct full tables/plans equal; protocol1,2discovery units;18resume tasks=17cached+1fresh validation |
| Docker reference low-NT empirical permutation |same-image direct full tables/plans equal;8resume tasks=7cached+1fresh validation |
| Injected distributed worker failure |marker present; no collector submitted; no completed discovery family published |
| Guards in both Docker harnesses |nested literal-child change, dynamic include, missing/cyclic config, legacy flag, output collision, changed seed and altered source all reject as intended |
| Installed executable example |prepare/execute/collect public APIs passed |

Both new Docker signatures name image1510c16, and the result manifests and both scientific reports carry the independently observed installed identity8f59580a. This binds the receipts to the new image; earlier169c16e/9d150d receipts were not substituted. The literal-local external include policy and host-side hashing rule are unchanged; dynamic/interpolated includes remain unsupported. The new guard harnesses passed for this image. No new external-body marker probe or nf-test run is claimed in this refresh.

New source/archive equivalence receipt reports176scoped files equal with source01bb39c/archive606363, and unchanged native/tests/DESCRIPTION/NAMESPACE relative to169. That is a worker/primary source-equivalence receipt, distinct from the independent actual-image archive readback.

## Historical saved-payload comparison

I inspected `/private/tmp/extract-pipeline-compatibility.R` and independently compared `/private/tmp/pipeline-bioc01bb-prior-comparison/{prior,current}.rds`: all72payloads are identical. These cover nine in-memory cases at eight saved stages, including standalone power/discovery afterQC, and compare result tables, realized plans/budgets, assignments and R RNGkind/state. This reuses saved historical169computation and compares it to new01bbcomputation; it does **not** rerun the old source, prove newly executed historical checks or substitute old results for the new image's same-runtime suites.

## Evidence and limits

Compact proof files: `verification.json`, `identities.json`, `pipeline-source-manifest.json`, `image-runtime-readback.json`, `scientific-receipt-summary.json`, `receipt-index.json` and `artifact-manifest.json`. Receipt index hashes the actual new logs, RDS reports, signatures, traces, source-equivalence receipt and paired payload artifacts.

The final evidence readback verified all19committed hash-manifest entries, with no missing files or hash mismatch. Thirteen original build/scientific/Nextflow/runtime/source/extraction receipt files are byte-identical to the temporary originals already reviewed. Installed-package and session-info originals are now committed and match the recorded R4.5.1/SCEPTRE0.99.0/ondisc1.3.1 versions. The committed execution record binds package evidence HEAD`9a4694450b5d456beedc11f3876ecb1cd84cb05c` separately from source01bb39c. See `evidence-commit-readback.json`.

This refresh qualifies bounded pipeline compatibility against the named archive/image on LinuxARM64. R full-check/Bioconductor qualification is handled by the separate primary review; this pipeline review does not rerun or independently certify those checkers or submission approval. No new hosted CI, nf-test, large-screen, cluster/cloud, Apptainer, AMD64 or cross-platform qualification is claimed. No image was published. Durable evidence is bound to cleanbfce9a while the qualified production baseline remains818a156. Any subsequent evidence-only commit should preserve these original bytes and identities.
