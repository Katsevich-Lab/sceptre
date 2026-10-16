# Same-environment workflow reference

The committed RF0 outputs record the original Mac numerical environment. The
unchanged original package already gives different NB dispersion estimates and
some downstream calibration decisions on hosted Linux/R/BLAS combinations.
Those differences precede the refactor. The stored files and their tolerance
remain unchanged.

CI therefore runs the immutable fixture inputs through original sceptre commit
`21f9ea098b69c444a884a49bc254de711968e3b5` in a fresh R process before testing
the candidate. It verifies every original R/native/data source SHA256 against
the frozen source manifest, rejects additional source files, forces a clean
native rebuild, and rechecks source and input hashes after execution. Candidate
code is never loaded in that reference process. No new random inputs or cases
are generated.

The temporary oracle includes workflows, cache sequences, continuation of the
frozen saved object, edge observations, and the expanded permutation workflow.
When ondisc is already installed, it also runs the original backed workflow
and records its separate import warnings/RNG state. The generator never installs
ondisc; the dedicated in-memory lane remains dependency free.

Tests require pinned and complete provenance, matching numerical dependencies,
source/input/harness/output integrity, and the current platform/BLAS setting.
The fingerprint also records effective R compiler/preprocessor/linker flags,
compiler executable identity/version, build environment flags, MD5s of effective
user/site Makevars and R configuration inputs, and `RNGkind()`. Cached readback
rechecks the live fingerprint and explicit debug mode before reuse. Batched
compiler probes are reused only while flags, configuration hashes and compiler
binary identity remain unchanged; RNG changes are checked on every readback.
User/site Makevars with includes or variable expansion are queried afresh so
changes in untracked included files or externally supplied values cannot hide
behind the compiler-query cache. Environment-field ordering uses fixed byte
order, independent of the test runner's reproducible-output locale.
The developer generator records SHA256 plus MD5; package-test readback uses
base R's MD5 support without additional package dependencies. A missing,
incomplete, stale, or unpinned oracle fails CI. A temporary oracle supplies only
scientific expectations; all inputs remain the committed fixtures.

Candidate scientific outputs use the original package's same-environment
outputs with the existing `1e-8` tolerance and exact IDs, order, budgets,
QC/significance decisions and stages. API, slot/formal inventories, index maps,
edge errors, and the independent frozen native sampler oracle keep their
portable fixed checks. The stored Mac scientific outputs are additionally
checked whenever their recorded numerical environment and debug build apply.
That qualification includes the compiler command/version recorded in the
native fixture and the workflow manifest's RNG kind.

The full package-check lanes compile the original with `debug = FALSE`, matching
normal package installation. The `test_local()` lane uses `debug = TRUE` on both
sides. Native optimization settings are part of the declared reference mode;
they are not used to loosen any comparison.

For a local paired run, archive the pinned original into a fresh directory and
run the following from the candidate checkout:

```sh
Rscript --vanilla tools/refactor/generate-runtime-oracle.R /path/to/original /path/to/fresh-oracle
SCEPTRE_REFACTOR_ORACLE=/path/to/fresh-oracle Rscript --vanilla -e 'testthat::test_local(".", stop_on_failure=TRUE)'
```

Set `SCEPTRE_REFACTOR_ORACLE_DEBUG=false` in both commands when validating an
installed package rather than a default debug `load_all()` build. Keep the
verified source archive available for the test readback. CI's
`prepare-runtime-oracle.R` fetches and archives the exact original commit into
the runner's persistent temporary directory and exports the oracle path to
later steps.
The archive command disables checkout line-ending conversion locally, and
source verification still checks the actual pinned bytes without normalization.
