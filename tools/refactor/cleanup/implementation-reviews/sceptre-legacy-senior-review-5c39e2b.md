# Independent legacy-boundary source review

Reviewer: rf0_senior_review. Date: 2026-10-16.

Reviewed `5c39e2be5f1c010c93212427da782ec028b7a140` versus `627e7e4c619b06359980b9e04fb246faa67c398d`: e25e2e2, 7b65cef, 7f0d3cc and 5c39e2b. Worktree HEAD, clean tree and all changed source/test/help blobs were checked. **No P0, P1 or material P2 source finding.** This is bounded source review, not final combined-branch qualification. No production edits, package builds, R test execution, plot execution/rendering or benchmarks were performed during the active count resource window.

## Runtime changes

- The guide histogram adds exponential bins only above the already constructed max(10, rounded_threshold-1) single-count bound. This prevents the characterized negative/NaN expansion request at high explicit thresholds without changing the ordinary branch's bin upper bounds or existing display labels. Existing label quirks are specifically retained, not silently repaired as a new binning method. Independent regression expectations cover both high thresholds and unchanged default labels/counts.
- Volcano rejection color now reads the existing significant flag directly. The horizontal p threshold remains a display line; clipping, transformations and table order remain. This fixes largest-significant-value ties without recomputing family significance or changing scientific results. The test intentionally includes equal p-values with different existing flags and verifies their exact colors plus nonmutation.
- QC forwards downsample_pairs, size and transparency to the actual pairwise plot; defaults are identical to the removed hardcoded defaults. The unused private return_indiv_plots argument is removed with its sole call. List/combined output and the no-discovery single-cellwise-plot exception are retained. Tests assert sampled row identities, layer settings, default/direct equivalence, return class, serialization and RNG.
- QQ setup selects the minimum ggplot group independently within each PANEL, retaining input row order. Equal-size validation is unchanged; beta order-statistic ranks/scales/grid sampling are untouched. The new two-panel regression checks both bands against independent ppoints/qbeta expectations and the old malformed-size error. It does not introduce a global point cap or new confidence-band method.
- Mixture perturbation validation now checks its own upper range against20, correcting the mistaken pi variable under the already stated error contract. Public defaults remain within that range. Other shape/NA/fractional/nonfinite behavior and validation order are retained. The processor helper is the same literal predicate/message at the same three validation positions, not a stronger general scalar validator. The unused lazy calibration_check_run formal had no actual caller argument or executable use; removing it does not force a formerly lazy expression. The unused guide-store read and misspelled histogram local are similarly bounded cleanup.
- Maximum assignment consumes two explicit source-cell summary vectors instead of the whole facade. Its sole dispatcher extracts the same slots; which(), <= fraction, < total, initial membership and union ordering are unchanged. Tests independently specify repeated requested guides, exact one-based memberships/flag order, nonmutation and RNG. Public assignment formals, saved fields, pipeline boundaries and scientific assignment rules remain.

## Duplicate theme removal

In plot_cellwise_qc the first get_my_theme() remains. Between it and the removed duplicate, the code only adds scales/labels/axis guide; no theme element is modified. Both calls return the same complete theme_bw-based object, and the explicit legend and axis-title overrides still follow the retained theme. Thus source inspection identifies no changed theme/visual rule and no reason to demand a revert solely because the duplicate was removed.

A tiny old/new theme/build-data comparison is **pending**, explicitly held at the coordinator's request until the count timing/profile window is released. This review does not claim that comparison or rendered visual parity has run. The lane's default pairwise plot comparison is relevant but does not independently compare old/new cellwise themes. Complete the tiny comparison after the quiet window before recording an executed default-visual-equivalence result.

## Contracts/help/evidence

The new contracts cover actual source/analysis/NT-local positions, repeated assignments/request order, validators' real acceptance/error behavior, storage effects and retained pipeline QC phases. QQ prose accurately says per-occupied-grid sampling; guide/target Bonferroni interpretation and a new global downsampling rule remain deferred. The five whitelisted Rd diffs only correct argument/return prose; names/usages/examples are unchanged, including QC's single-plot exception. Documentation commit7f0d3cc is reported with five-file AST parity; source diffs inspected here are consistent with documentation-only changes in that commit.

The retained five original-source reproduction script/JSON/RDS/log SHA256 values independently match the report manifest. Their reported observations motivate the precise fixes, not a general claim about all plots. The worker reports39 tests/808 expectations on a fresh immutable -O2 install, including existing scientific/saved/backed/fork contracts; those tests were not rerun by this reviewer. They predate the count allocation remedy and cannot qualify the final integrated branch. Final immutable package/CI/consumer/resource checks and the small held theme comparison remain outstanding; no performance acceptance or pristine-source claim is made.

## Reviewed file SHA256

- `R/assign_grna_functions.R`: `d9ce47de1d9ba371ef366deda9c26744d6616d220273b3e19c0f32c688a5f48c`
- `R/check_functions.R`: `3b021705dcf19128b263e871fe65ee75bd45ed0a8630b238b859367560731f34`
- `R/plotting_functions.R`: `66b1547678d785da9d14ed2fd784ac8d5ed42bb88059a73742bb65c8de98bcfa`
- `R/qq_plot_helpers.R`: `348fc17cde0864bea749054fa755324d7f1b24dd746395cd6eb8c6b19ec7d568`
- `R/s4_analysis_functs_1.R`: `52df97f69fe6f5a9b9876ac4bbccb92fab7c03df88b1ad31d8d8ce1c73b8d3d1`
- `man/plot_assign_grnas.Rd`: `238cc947609a012f4e0e22aefaa86391fa7165cfe8579d95c54ad48b8c5ba33b`
- `man/plot_covariates.Rd`: `ae932e9556e87bdf1424461f9ad14b8fbc5aa14827d652d5dade48303177c575`
- `man/plot_run_calibration_check.Rd`: `4deca21e047ff64067143ff70db4ba769cea4979be7b70a8ffe73851856680bc`
- `man/plot_run_discovery_analysis.Rd`: `2295beab96e9b78d9b68aa803bb9d8123d5af9dc98c8813707f052c78ccd2877`
- `man/plot_run_qc.Rd`: `6324c26263ec32c9266f61957e2a6ccda2cc88de83f708ab11fdb2da027983fd`
- `tests/testthat/test-refactor-legacy-boundaries.R`: `2c2f8c99b2ef852ddadf9b5f6f2c891ff931fff39fb050d9f5f0f9b7c5421506`
