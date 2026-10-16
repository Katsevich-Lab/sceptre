# Assignment getter source probe

24 independent source-evaluated checks and six exact original/candidate cases pass. Saved inputs remain serialized-identical; single/multiple/empty guide memberships, zero-based slots, repeated calls and existing post-QC one-guide logical-vector dimension drop are covered. The baseline source is pinned to 21f9ea0. This is not an installed-package qualification.

invalid-fixture-attempt.log retains the first harness attempt: its one-guide memberships were unsorted, producing the pre-existing invalid ngRMatrix extraction error, and its expected pointer retained names absent in the production slot. The corrected fixture uses valid increasing memberships and an unnamed expected pointer. No production sorting or new behavior was introduced.
