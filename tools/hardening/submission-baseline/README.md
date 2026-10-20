# Actual Bioconductor submission baseline

The submission baseline is **21f9ea098b69c444a884a49bc254de711968e3b5**,
submitted August 31, 2026. It is not the later internal-refactor baseline.
The original official run used **BiocCheck 1.49.30**, Bioconductor 3.24 and
R 4.6.1, and reported **0 errors, 1 warning, 3 notes**.

The authoritative chain is:

1. [Submission issue 176](https://github.com/Bioconductor/BiocContributions/issues/176)
   records the initial clone, first build and authors' diagnostic explanation.
2. [First source build](https://github.com/r-universe/biocstaging/actions/runs/33438717492/job/99641478053)
   explicitly resets to the submitted commit.
3. [Original Bioconductor check](https://github.com/r-universe/biocstaging/actions/runs/33438717492/job/99642866650)
   supplies the preserved `official/00BiocCheck.log` and `official/output.log`.
4. `official/sceptre_0.99.0.tar.gz` is the original hosted source artifact,
   SHA256 `5300b3e3adfbb73a9a2173384ea3ec7e57b8cad8612f3fdc8116d5db5156eaf4`.
   Its `RemoteSha` identifies the submitted commit. The 75 checked packaged
   R, native, man, vignette and NAMESPACE files match Git exactly. DESCRIPTION
   differs only through build-time wrapping and R-universe metadata.

The live Compass preflight searched `sceptre Bioconductor` and
`sceptre submission`, then read the exact `projects/Bioconductor submission.md`
note. That note identifies issue 176 and the local planning tracker. Its
historical tracker says five notes; the actual submitted official build and
authors' September 1 issue comment say three. The official build is primary.
No Compass content was changed.

## Matched checker runs

| Source | Checker | Local ecosystem | Errors / warnings / notes | Long functions | Long lines |
| --- | --- | --- | --- | --- | --- |
| Original official submission | 1.49.30 | Hosted R 4.6.1 / Bioc 3.24 | 0 / 1 / 3 | 49 | 145 |
| Exact original source artifact | 1.47.26 | macOS ARM64 R 4.6.0 / Bioc 3.23 | 0 / 1 / 3 | 48 | 145 |
| Exact original source artifact | 1.49.30 | macOS ARM64 R 4.6.0 / Bioc 3.24 | 0 / 1 / 3 | 49 | 145 |

The two local runs preserve separate process logs, complete checker logs,
serialized/JSON result objects, session details and commands in
`matched-1.47.26/` and `matched-1.49.30/`. They use the exact official wrapper
options `no-check-R-ver` and `no-check-bioc-help`; no coding, formatting,
documentation or dependency check is disabled. The original wrapper is pinned
to r-universe action commit `e859fe2adc2553a5c8b0b2ec6044215d60d3d9c7`.

BiocCheck 1.49.30 was installed into an isolated temporary library from official
source commit `ca5c23e019dd695b7a5e29286f593a0e9b13639a`; its source archive,
DESCRIPTION and install log are retained in `checker-1.49.30/`. No user library
was modified. The original hosted R/Linux environment was not recreated; the
local checker's version and flags match, while platform/dependency differences
remain explicit. Candidate qualification must use the same isolated library.

An independent audit also ran the unchanged submitted archive under 1.47.26
with **all default checks enabled**. That run reported **1 error, 1 warning,
5 notes**, including a support-site HTTP timeout and the unchanged R-minimum
and mailing-list notes. Its exact invocation, terminal log, complete diagnostic
JSON and result are retained in `external-audit/unsuppressed-1.47.26/`. This is
not a clean unsuppressed check. The official-wrapper comparison above retains
the original submission's two documented options; it must not be presented as
all default checks passing. The independent authority/fix-scope reports are
preserved verbatim under `external-audit/` with their original manifest.

The 1.47.26 function-length checker overwrites duplicate function names. The
1.49.30 checker retains both `compute_group` methods and flags the inherited
52-line `StatQqBand` method, explaining the one-function difference. Inventory
scripts preserve duplicate entries by index.

## Diagnostic instances, not just categories

The existing warning is absence of a direct Bioconductor dependency. The
authors justified it in the September 1 issue comment; it is not a new warning
from this refactor and no meaningless dependency is added to silence it.

The original display note covers eight `cat()` calls in the S4
`print(sceptre_object)` method and one `print()` used to serialize the summary
under `sink()`. `inventories/display-call-lineage.csv` maps each original call
to the pre-cleanup qualified source. The summary call moved from
`R/s4_analysis_functs_2.R:723` into `R/output_resources.R:18`; cleanup added
connection/sink restoration without creating a new logging use. This is
inherited functional lineage, not an additional diagnostic instance.

The original long-function note has complete inventories for both checkers.
`inventories/qualified-long-function-lineage.csv` compares the **pre-cleanup**
qualified archive from source `169c16e` with the submitted source. It records
explicit rename/move evidence, so renamed permutation/CRT workers and the
S4 print method are not counted as new merely because their names/lines moved.
Under 1.47.26 the 69 current sites comprise 44 inherited sites, 22 new sites
and 3 new threshold crossings. The latter are the formula builder (47 to 52),
backed object reader (39 to 63), and show method (50 to 54). Extracting a new
long helper from an old long body is still classified as a new site. Semantic
lineage does not assert AST identity or excuse new behavior within a function.

The pre-cleanup archive has 473 long lines: 110 exact unchanged submitted line
instances and 363 new or modified instances. The complete source-text inventory
and conservative exact-text delta are in `inventories/`. A lower total or an
old formatting-note category alone would not establish a clean delta.

The following pre-cleanup diagnostics were absent from both matched baseline
runs and must be resolved:

- `1:3` in `R/distributed_execution.R` at lines 290 and 343.
- Paste inside the progress-message signal in `R/execution.R:42`.
- Three superassignments: `R/perturbation_fits.R:51` and
  `R/schema_migration.R:389,393`.
- Missing runnable examples for `prepare_execution`, `execute_work_unit`, and
  `collect_execution`.
- New or newly long functions and new/modified overlong lines described above.

Line references in this section intentionally identify the frozen pre-cleanup
source. The final candidate comparison below records the actual cleanup result.

`review-helpers.R` and its log preserve an independent review of the structural
helper extraction against the frozen pre-cleanup source: 40 exact comparisons
passed across low-MOI complement/NT permutation and high-MOI complement CRT,
including stage identities, plans, ordered contexts, unit values, collection,
RNG/caller immutability and representative failure messages. This comparison
loads functions in isolated environments over the same pre-cleanup installed
class/runtime. It does not replace final frozen-source installed validation.

## Frozen cleanup comparison

Implementation commit **01bb39ccf8373ae8965a2679ee4813a63ea2d85d** produced
archive SHA256
`606363fce32238354fc0d126d79f515973da7c8287ca43b14ce53ce2d00df4b7`.
BiocCheck **1.49.30** on that exact archive completed with **0 errors,
1 warning, 3 notes**, using the original official wrapper options and the same
local checker environment as the matched original-source baseline.

The remaining warning is the original dependency warning. All nine flagged
display/serialization calls map to the original submitted sites. The remaining
45 long functions map to already flagged submission functions; the remaining
108 long lines reuse exact submitted source text, with multiplicity preserved.
There are **zero new flagged function sites and zero new or modified long
lines**, as well as no new diagnostic categories. The literal sequence,
condition-paste, superassignment and missing-example notes are absent.

Complete candidate receipts and every remaining instance's original source
location are in
[`../submission-cleanup/qualification/remaining-diagnostic-map.json`](../submission-cleanup/qualification/remaining-diagnostic-map.json),
[`instance-gate/`](../submission-cleanup/qualification/instance-gate/) and
[`bioccheck-1.49.30/`](../submission-cleanup/qualification/bioccheck-1.49.30/).
The source helper comparison was repeated against this exact frozen archive;
all 40 comparisons passed. Final installed tests, executable examples,
historical scientific parity and R CMD check have separate integration-owner
receipts and are not inferred from BiocCheck success.

`baseline.json` records source/runtime scope. `files-sha256.json` pins captured
artifacts without rewriting their bytes. The original submission archive and
checker source are historical evidence and must remain separate from any new
candidate archive.
