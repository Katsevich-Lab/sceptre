# Independent W08 scientific source review and tail diagnosis

Reviewed W08 source files are identified by SHA256 in the adjacent manifest;
W08's checkout was still uncommitted, so this is not a revision-wide certification.
No W08 checkout file or ordinary installed library was modified by W03.

## Null and scope

The fixed-total multinomial example has exactly1000 UMIs per cell, so CP10K
is exactly a factor10 and log1p is strictly monotone. Its null probability is
.01 in both arms in the null/effect scenarios. Increasing the effect probability
from.02 to.04 takes mass only from the reservoir. Therefore the tested null
marginal remains Binomial(1000,.01), exchangeable across the fixed independent
cell arms even in the planted-effect scenario. This is a valid known null for
the relative-distribution rank test. The two treated guides are unioned; NT
controls are distinct. It does not assess stochastic depth, assignment errors,
clustered/donor dependence, broad integration strategies or causal raw-count
inference. The confounded scenario makes batch prevalence .9 versus .1 and null
probability .04 in batch1; conditional on batch there is no treatment effect,
but the marginal distributions differ. Its rejection is expected and must not
be labelled type-I error under a true rank-distribution null.

Independent screens (not genes or cells within one screen) are the MC units.
The planned200 screens per scenario give theoretical MC SE .0154 at rate.05.
Wilson intervals stratified by scenario/response are appropriate pointwise
intervals; screen seeds are paired across scenarios, so differences should be
analysed as paired if compared and scenarios must not be pooled as independent.
A raw-alpha summary is not a two-pair BH/FDR calibration certificate.

The pilot source has no per-replicate checkpoint/error capture and permits an
existing output directory. It therefore needs a fresh explicitly selected
output directory, terminal verification and preservation of the exact installed
candidate source/binary identity. SessionInfo and script hashes alone identify
neither the candidate commit nor binary (all program versions currently share
package0.99.0). W03 has not executed the W08 candidate package or200-screen pilot.
Statistical acceptance remains conditional on actual terminal results, scoped to
an experimental unadjusted cell-level relative-distribution method.

## Independently reproduced numerical limitation

The adjacent source script regenerates the actual seed821 fixed-depth count
outcomes and computes the rank statistic and tie-corrected variance directly,
without invoking sceptre. On this R4.6.0 host:

- Base `wilcox.test(exact=FALSE, correct=TRUE)` two-sided p=0.
- Independent continuity-adjusted z=17.19461.
- Stable normal-tail arithmetic `2*pnorm(-abs(z))`=2.9143e-66.

Local primary authority: inspected the actual installed stats internal
`.wilcox_test_two_pval_asymp` body. Its two-sided branch uses
`2*min(p <- F(z), 1-p)`. For this positive statistic, F(z) rounds to1, so the
subtraction returns0. This is **cancellation**, not representational underflow
of a2.9e-66 probability. Ties select the normal-tail branch, with the usual
continuity adjustment. The stable value is another evaluation of the same
normal approximation; it is NOT an exact finite-sample p-value, nor independent
validation of that approximation's far tail.

Acknowledge the inherited base-R/version-specific numerical limitation if the
experimental entry point intentionally delegates exactly to base stats. Do not
claim calibrated extreme tails or general FDR from base equality. Any stable-tail
replacement would be an explicit opt-in method behavior decision with its own
R-version/tie/sidedness tests, not an unnoticed patch to the current default.
This lightweight400-cell rank diagnostic was run without package compilation,
large simulation, parallel workers or taking another worker's heavy lease.
