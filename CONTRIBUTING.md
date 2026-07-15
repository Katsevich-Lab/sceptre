# Contributing to sceptre

Contributions to `sceptre` are welcome, including code, tests, and documentation.

Participation in the project is subject to the [Code of Conduct](CODE_OF_CONDUCT.md).
Please do not include private data, credentials, protected health information,
or other sensitive information in pull requests, examples, or test fixtures.

For bug reports, feature requests, and questions about using `sceptre`, see
[Getting help with sceptre](SUPPORT.md).

## Contributing code and documentation

Small bug fixes and documentation corrections may be submitted as pull requests without prior discussion.

Before implementing a substantial feature, public API change, or change to
statistical methodology, start a [feature-request discussion](https://github.com/Katsevich-Lab/sceptre/discussions/categories/feature-requests).
Early discussion helps the maintainers assess scope, compatibility, and
validation requirements before substantial work is invested. Maintainers may
promote a sufficiently scoped proposal to an issue for implementation tracking.

### Code, tests, and documentation

- Keep changes focused on a single problem or feature.
- Add or update tests for behavior changes and bug fixes.
- Follow the existing R style. The repository's `air.toml` file records the formatting configuration.
- Update user-facing documentation when public behavior changes.

Changes affecting model fitting, calibration, power, numerical results, or other statistical behavior generally require more than unit tests. Describe the relevant simulation, benchmark, or reproducible analysis in the pull request.

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

Contributions accepted into the repository are distributed under the project's [GPL-3 license](LICENSE.md).
