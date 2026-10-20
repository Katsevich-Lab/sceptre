# Owned cleanup checkpoint (pre-final)

This evidence covers the four R files and three generated help pages listed in
`checkpoint.json`. Other agents concurrently changed other files. This is a
source-level checkpoint, **not** final archive installation/BiocCheck qualification.

Changes:

- Precompute the worker log filename pattern outside the condition signal.
- Capture perturbation warnings in a function-local environment; keep their
  order, caller propagation and stored diagnostics.
- Return migration stage transitions explicitly rather than superassigning.
- Use `seq_len(3L)` for output-level membership.
- Give all three distributed APIs a standalone runnable bundled-data workflow.
- Decompose long routines into helpers with scientific phase boundaries:
  physical schema admission, field normalization, component copy, request
  reconstruction, historical results and readiness; implementation certificates,
  manifest/proposal validation, plan selection, numerical unit work, compact fit
  selection, ordered proposal ownership, cache reconciliation and result assembly.
- Keep collector updates local until all admission, fit checks and full-family
  result assembly finish; retain the same result/fits commit sequence.
- Wrap owned comments/help and split long literal `stop()` messages into arguments.
  `condition-text.log` verifies concatenated literal condition text against HEAD.

`function-lengths.log` uses BiocCheck's own parsed-source length implementation.
Every owned R function is <=50 source lines; all seven files are <=80 columns.
No diagnostic suppression, native/statistical algorithms, budgets, seed choices,
row traversal or public messages were deliberately changed. The actual-code
execution fingerprint necessarily changes with implementation source; saved
execution manifests from different builds remain incompatible by design.

Focused regression logs cover distributed/direct exact parity, permutation/CRT,
calibration/power/discovery, reversed and corrupt proposals, ownership, warm cache
conflicts, worker failures, backing/reopen and all saved-object migration layouts.
The generated Rd examples execute preparation, every work unit and collection.

Reproduce from the package root with `Rscript --vanilla
 tools/hardening/submission-cleanup/owned/reproduce.R`. Final package-wide and
submission-baseline comparison evidence belongs to the integration owner.
