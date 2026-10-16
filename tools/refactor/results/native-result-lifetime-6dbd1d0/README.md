# Scoped native result lifetime evidence

Production 11502f871ebeb78d7cf084c45ffdb895d74ac937;
probe 6dbd1d0621fd88716423f8d77cacf176dce66fd3.
`evidence.tar.gz` retains 106 files, including every bounded attempt, actual
original/transformed source, corrected source-only original oracle boundary,
compile logs, source/DLL checks, mechanical transform recovery and final results.
Archive SHA256: 5cab85f53155f3284ab1083ed4e6e3ab23f0ea0a19dbfb8a7bb4142acb44e0d9.
The unchanged per-file manifest is also readable here; its SHA256 is
d5e8ecc8c65437435eacb332ecc4ded9ae255c9c2aa566155b767851d24ddd3f.

The actual-source probe passes 959 cases across 12 workloads: 12 successful exact
outputs, 858 controlled raw R errors and 89 native allocation errors, plus one
separate safe unsafe-order negative control. Every result/scratch payload is
released at the declared boundaries. Component tests pass 112 expectations in
three blocks, covering all 17 frozen native cases. No actual OOM is attempted.
All-rejected requests and stage 1 with huge unused B3 exercise exact/reached-only
allocation. Full scientific/RNG/installed/platform qualification remains separate.

The final driver accepts verified already-compiled original SOURCE, not an
installed-library guess. It checks the complete frozen 21 source inventory,
namespace path and loaded-DLL equality, then obtains only original group/pair
expectations. Original fixtures/expected tail cases are unchanged. CI receives
this source path from the existing original-runtime-oracle bootstrap, avoiding
another full original build.

The first attempt exposed a duplicate protected dimension read under native
ownership; its minimal pre-owner-validation repair is retained. A later report
copy overwrote instrumented source files: exact transforms were recovered and
matched against their recorded pre-overwrite hashes, then the final runner
separated original/transformed directories. Independent review verified recovery
and all final retained hashes; that evidence P2 is closed, not omitted.

[Source/design/probe review](../../cleanup/native-result-lifetime-review.md)
and [signed closure](../../cleanup/implementation-reviews/sceptre-native-output-evidence-closure-6dbd1d0.md)
qualify named native payload/scratch and raw conversion leaves only. Pre-leaf
local/input Rcpp precious preservation and universal Rcpp/OOM cleanup are outside
that scope; no such global guarantee is made. No unused-stage or rejected-request
output capacity, statistical change, new method or new R/Rcpp minimum is added.
