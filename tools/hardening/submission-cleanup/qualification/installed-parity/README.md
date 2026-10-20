# Installed archive parity qualification

Candidate: source commit `01bb39ccf8373ae8965a2679ee4813a63ea2d85d`;
archive SHA256 `606363fce32238354fc0d126d79f515973da7c8287ca43b14ce53ce2d00df4b7`.
The archive digest was independently verified before and after these checks.

The unchanged `tools/hardening/validation/run-parity.R` ran with:

- `SCEPTRE_PARITY_LIBRARY=/private/tmp/sceptre-bioc-cleaned-01bb39c/sceptre.Rcheck`
- `SCEPTRE_PARITY_OUTPUT=tools/hardening/submission-cleanup/qualification/installed-parity`
- `SCEPTRE_PARITY_ARCHIVE_SHA256` set to the full digest above.

All seven scientific cases match the existing matched-budget reference at zero
tolerance, with 33 comparisons per case (including warnings, RNG state, results,
preparation payloads and response fits). See `parity-summary.json` and the raw
RDS artifacts. These are bounded regression cases, not a scale or coverage claim.

A separate fresh R process loaded this exact installed package (asserting the
resolved package path), extracted its help with `tools::Rd_db`, converted all
three distributed API examples with `tools::Rd2ex`, and executed each standalone
workflow. All three passed; `examples.log` records package path and session.
No source-loaded package was substituted for this installed-example check.

`receipt.json` records terminal status and identity. Source file hashes still
match the 190-file full-source-suite before manifest. R CMD check and final
BiocCheck/baseline diagnostics are independently owned by the integration team.
