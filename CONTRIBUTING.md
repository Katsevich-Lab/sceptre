# Contributing to sceptre

Contributions to `sceptre` are welcome, including code, tests, and documentation.

Participation in the project is subject to the [Code of Conduct](CODE_OF_CONDUCT.md).
Please do not include private data, credentials, protected health information,
or other sensitive information in pull requests, examples, or test fixtures.

For bug reports, feature requests, and questions about using `sceptre`, see
[Getting help with sceptre](SUPPORT.md).

## Contributing code and documentation

Small bug fixes and documentation corrections may be submitted as pull requests
without prior discussion.

Before implementing a substantial feature, public API change, or change to
statistical methodology, start a [feature-request discussion](https://github.com/Katsevich-Lab/sceptre/discussions/categories/feature-requests).
Early discussion helps the maintainers assess scope, compatibility, and
validation requirements before substantial work is invested. Maintainers may
promote a sufficiently scoped proposal to an issue for implementation tracking.

### Code, tests, and documentation

For internal ownership, cell-coordinate conventions, numerical boundaries and
component compatibility, see [Internal architecture](ARCHITECTURE.md).
Permanent tests in `tests/testthat/` check current package behavior without an
old installation. CI also runs [current-source native checks](tools/checks/README.md).
Temporary [refactor review tools](tools/refactor/README.md) hold the explicit
old-versus-new comparison and downstream/resource evidence. They are excluded
from package builds and never run in CI.

- Keep changes focused on a single problem or feature.
- Add or update tests for behavior changes and bug fixes.
- Format handwritten R code with Air as described below.
- Update user-facing documentation when public behavior changes.

Changes affecting model fitting, calibration, power, numerical results, or other
statistical behavior generally require more than unit tests. Describe the relevant
simulation, benchmark, or reproducible analysis in the pull request.

### Formatting R code

Use Air 0.9.0 as the single R formatter. The repository's `air.toml` specifies
four-space indentation and a preferred width of 80 characters. Format only the
handwritten package source, then check it with the same command used by CI:

```sh
air format R
air format --check R
```

Keep Air's default exclusions for generated bindings. Immutable fixtures and
`tools/` (including reference code and retained results) are explicitly excluded;
do not recursively format tests or reference/evidence directories.

Air does not guarantee that every line is at most 80 characters. Wrap ordinary
prose comments and documentation separately where useful, preserving code
examples, fenced blocks, URLs, mathematical notation and executable strings.
BiocCheck remains a separate Bioconductor checker; Air does not replace its
checks or guarantee that every style note is resolved.

### Before submitting a pull request

Before submitting a pull request, run a package check:

```r
devtools::check()
```

Open the pull request against `main` and explain:

- the problem being addressed;
- the proposed change;
- the tests and checks performed;
- relevant documentation changes; and
- any statistical, compatibility, or performance implications.

A maintainer will review the contribution and may request revisions, additional
tests, or further validation before merging. Maintainers may decline changes
that fall outside the project's scope or current roadmap.

## License

Contributions accepted into the repository are distributed under the project's
[GPL-3 license](LICENSE.md).
