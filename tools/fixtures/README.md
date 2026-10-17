# Maintaining test inputs

`generate-inputs.R` recreates four deterministic mathematical/workflow input
files in a fresh directory outside the checkout. It never installs SCEPTRE,
executes a historical implementation or generates expected analysis outputs.
See [fixture provenance](../../tests/testthat/fixtures/README.md) before updating
an input. This is permanent maintainer tooling, excluded from the R source
package; CI reads the committed inputs and does not regenerate them.
