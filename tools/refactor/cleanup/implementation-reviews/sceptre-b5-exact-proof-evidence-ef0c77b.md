# B5 exact-source proof evidence addendum

Reviewer: rf0_senior_review. Date: 2026-10-16.

After the source review (`/private/tmp/sceptre-b5-senior-review-ef0c77b.md`, SHA256 7d5b85cda34cf14717723de5c42edeac45484f12c096f2a471c25765333736fd), I inspected the completed exact-run log and read only the retained result metadata. I did not rerun the test, compile code or edit production.

The result binds the original/candidate/probe/runner/input digests to the same immutable source reviewed at ef0c77b, including the runner with stage-reference crosschecks:

- Original: 4add3e36820f9af6b794b715862c0981f44cb51f11750cf5107d39c2cd4ee776
- Candidate: 27ea8665a984dcef83f35c646ab4b38989ea8cacafa28975aa7ffdedd7f9aeb2
- Probe: f9c8a6e991adce64e09446c876ec4a25cf3bb5999328270de6453583124b29cc
- Runner: 6968a8278bba86da49529bdf1e0c63db4526ed6e29c5319be5d446b540a1fa6d
- Inputs: 0d30f0ba5612e920a136e9f6236e2fef5469b6d4811c4b58de274db7835e6f2f

The retained payload contains 107 numerical and 14 controlled paired cases, and the completed log reports that all pass. Policy facts give binary64 radix2/digits53/effective precision53, long-double evaluation digits53 on this host, domain/pole/overflow throw_on_error and underflow/denorm ignore_error. Thus the earlier pending exact local test-only rerun is now supported by pinned evidence. Installed package/scientific workflows and remote supported lanes remain separate gates; this addendum does not claim them passed.

Evidence hashes:

- `/private/tmp/sceptre-dispersion-fault-proof-ef0c77b/dispersion-exception-results.rds`: `a3f572d2bc2ebd9d3ad1ad2b8aebb784c00d2fe3c790e6b4477f23dad414c23d`
- `/private/tmp/sceptre-dispersion-fault-proof-ef0c77b.log`: `f11d89acfa82423566ef2ab54058c30ecfb8f9481fce3fa524d75f32915df086`
