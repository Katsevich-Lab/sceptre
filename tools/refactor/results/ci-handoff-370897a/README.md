# CI and maintainer handoff validation

`qualification.json` records local permanent checks at 708f37a and the explicit
old/new command at 370897a. It binds every permanent source/test/workflow input
to unchanged SHA256 bytes across those commits; only comparison wrapper and Git
attributes changed. Production remains unchanged from 0a9f9d8.

The archive retains logs, detailed test results, original runtime observations,
current native instrumentation/results, and failed authoring attempts. It does
not contain duplicate installed libraries or full source archives. Recreate those
from the immutable Git refs recorded in the receipt. No new performance or
consumer qualification is implied. Hosted checks are recorded separately in
`../../CI_HANDOFF.md` after completion.

To restore and verify:

```sh
mkdir /absolute/fresh/evidence
tar -xzf local-evidence.tar.gz -C /absolute/fresh/evidence
```

Compare the archive SHA256 and each restored file with the receipt. The creation
script read back every archive member and verified all recorded hashes.
