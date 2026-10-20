# Contributor site sources and review status

The contributor pages are part of the existing pkgdown package website. `CONTRIBUTING.md` is the entry point and retains the contribution policies. `ARCHITECTURE.md` is the internals guide; its primary visual is the object and workflow map.

## Source authority

The component schema, object fields, and legacy field inventory remain authoritative in `tools/components/schema3/`. The script `tools/architecture/prepare-pkgdown.py` generates website-only R Markdown articles from these and their linked supporting documents, rewriting internal document links to rendered article pages. Generated articles are not separate maintained copies. Run the generator after updating canonical Markdown, then rebuild pkgdown. The existing `.Rbuildignore` excludes `vignettes/articles` and maintainer tools from the R source package.

The main [discovery walkthrough](contributor-walkthrough.md) and [change map](contributor-change-map.md) are prose reading guides. They do not run an analysis or establish an executed worked extension.

## Figure provenance and scope

The object map is adapted from the newer October 19 interactive SCEPTRE visualization. It was checked against the combined stage-one candidate: five persistent owners (`ScreenData`, `AnalysisSpec`, `PreparationRecord`, `FitCache`, `ResultRecord`), transient `RunPlan` and `CellView`, import populating `ScreenData` only, and calibration/power/discovery as sibling analyses. The copied figure uses the current export function name and resampling wording covering permutation and CRT. It is a conceptual illustration, not an execution trace.

The older module map remains a complementary reference for requests between functions; it is not labeled incorrect or used as the central object model. Historical supporting receipts and old walkthroughs retain their named revisions and limitations. Their presence on the site does not qualify a later candidate.

## Remaining review

Approachability, a worked contribution exercised by a newcomer, and fresh-checkout setup still need review before contributor readiness is signed off. Website rendering and package/scientific qualification are separate checks. Nothing here publishes the site or establishes hosted/platform qualification.
