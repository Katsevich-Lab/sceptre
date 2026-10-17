# Permanent and temporary documentation layout

The refactor specification moved from the root into `tools/refactor/`, with
relative links fixed and its historical-plan role stated. The contributor guide
maps permanent package/maintainer material and temporary refactor records.
The lasting architecture guide no longer depends on the temporary directory.
The now-unused root-spec exclusion was removed; the whole refactor directory
remains excluded from R package builds. Frozen historical evidence was retained.

Validation used the immutable Git tree in `validation.json`, after the layout
changes and before this receipt was added. `R CMD build --no-build-vignettes
--no-manual` succeeded. The resulting tarball excludes both the temporary
refactor material and developer-only permanent guides/checks, and retains the
package's permanent regression tests. All packaged R/C++/test files match the
checkout. Across the reorganization, 181 implementation/test/CI/configuration
files remain byte-identical to the preceding commit. Current relative document
links resolve. No full package check or runtime/resource rerun is claimed.

The source tree, exact command flags, package hash, build log and tarball file
list are recorded here. All files in this receipt remain outside the package.
