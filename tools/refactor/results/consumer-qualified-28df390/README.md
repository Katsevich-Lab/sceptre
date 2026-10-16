# Final consumer qualification: immutable21f9ea0 versus28df390

Read qualification-comparison.json first. Actual30-script payloads and the
explicit fair/serialized local Nextflow control are exactly equal, with exact
identifiers, decisions, schema, ordering and actual pod/arrival plans. Floating
comparisons also pass1e-8. Source, installed-library and native-build provenance
is pinned and verified; seven thread environment variables were fixed1.

Both original two-task Nextflow workflows execute successfully, but their
ordered gathered output comparison fails because completion order changes
factor levels and tied-row ordering. Serialization alone does not fix emission
order. Both failed comparisons are retained; all five raw worker-pod payloads
are exactly equal. The qualified control adds fair=true before the original
collect, without editing pipeline sources or sorting published outputs.

Book runs each pass14/15 unchanged uncached cells. Both preserve the identical
known CellRanger example-data failure; no whole-book success is claimed.

Reports retain absolute paths of original executions. Portable comparison
payloads, source runners/inventory, build manifests/logs and task trace/logs are
copied here for durable review. Large raw work directories/plots/source archives
remain in the parent staging directory. artifact-manifest.json covers every
file here except itself; its external SHA256 is reported to the parent agent.

This evidence is specific to production28df390, not later source edits.
