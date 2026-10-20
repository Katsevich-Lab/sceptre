# Frozen foundation platform receipt

Source0b901c059a9ce9e55234d672587cfd5d523e7383, archiveSHA256
f9172ba52aeab123d5473a0b189d877527704d86cc6019ff7b37fb69c2d45a6f.
No real-screen inputs were used; these are existing public/synthetic package tests.

| Lane | Actual runtime | Passed | Failed/errors/warnings | Explicit skips | Scope |
| --- | --- | ---: | --- | ---: | --- |
| Macarm64 | R4.6.0/testthat3.3.2 | 5645 | 0/0/0 | 0 | Build, fresh archive installation, full packaged suite, R CMD check Status:OK |
| LinuxAMD64 emulated | R4.5.1/testthat3.2.3 | 5622 | 0/0/0 | 1 | Same archive/private installation/full suite; rhdf5 corruption fixture skipped |
| Linuxarm64 floor | R4.1.3/testthat3.1.8 | 5509 | 0/0/0 | 10 | Same archive/private installation/full suite; ondisc absent |

Original matched checker1.49.30/Bioc3.24, exact official two flags:0errors,
1warning,3notes; zero new categories or flagged instances. All45 long functions,
108 long lines and9 display sites inherit from the actual submitted source21f9ea0.
Complete function, line and display lineage CSVs are included. The final bound
invocation/summary/archive hashes agree. Four negative probes correctly reject
wrong archive, stale summary, weakened options and a new message inside an
inherited diagnostic category. Initial unbound gate remains only in local history.

The final summary-hash binding was finalized after terminal checker output using
the exact new finalizer statements; this did not rerun inference or installation.
See provenance.json. Local full logs, source archive, dependency inventories,
image manifests and earlier receipts are SHA-bound in source-and-artifacts.json;
raw libraries/RDS/archives/logs stay outside the publication branch.

This is foundation qualification, not the later final integrated candidate or
scientific calibration certification. No current Bioconductor ecosystem, hosted
Windows/devel, nativeAMD64, final SCE or accelerated-fitter support is inferred.
Actual baseline PR239 separately checks source840f334; it is not this revision.
Repeat the parent README's exact-archive commands on the final candidate and
refresh images/dependencies when Imports/LinkingTo change.
