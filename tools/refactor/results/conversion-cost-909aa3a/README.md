# Revised conversion qualification at 909aa3a

Original 21f9ea0; candidate 909aa3ad76ee21e7f9ab83a63c645bc0fb6a43c7.
The staged source manifest SHA256 is
`a2d085699bbe3df5e150c83a496e44f2e6c5057d75de050cd6d41c628f4ab822`.
It is retained unchanged. `transport-manifest.json` maps all 65 staged files (including the unchanged source manifest) to
committed bytes. JSON larger than 100kB is losslessly gzip-compressed; original
uncompressed SHA256 and size are retained and verified after decompression.
Large experimental inputs/outputs, archives and installed libraries remain at
the original recorded paths; this bundle contains complete trial/sample/phase/
allocation records, plans, source/build provenance and bounded correctness data.

Triplet: 16 warmups, 64 primary trials, 16 separate profiles and 32 five-call
batch children. Sixteen conversions/eight exact output pairs pass and inputs
remain unchanged. All four meaningful timings pass (+8.27%, +6.93%, +7.55%,
+4.39%). Small cases and large logical memory gates pass. Large numeric row/
column whole-child increases are 121.938/123.141MiB; both **cross the unchanged
generic gate**. Eugene subsequently accepted these measured bounded costs in
[decision D2026-10-16-02](../../cleanup/decisions.md). Neither the original
proposal nor pre-acceptance reports have been rewritten to pass. Three short
warmup/timing children lack active-phase samples; all eight-million-entry
children have samples. OS whole-child RSS is the primary memory authority.

Dense: all eight conversions/four exact pairs and four OS memory gates pass;
net large-path increases versus original are 50.562/49.453MiB. Profiles support
an extra full dense payload copy through the public general-shape conversion.
All timings are inconclusive (8–75ms); 19/40 short timing/warmup children lack
an active-phase sample. No sampled peak was fabricated. This is a bounded net
cost comparison, not isolated workspace or a timing waiver. Independent review
is retained under `reviews/`.

Scope correction to the original triplet proposal's completion paragraph:
“native exceptional ownership ... complete” applies only to the specifically
controlled CL2 count/construction and prior sampler ownership boundaries. It
does not qualify every native return path. The fresh whole-package review
subsequently identified separate retained calibration and inference return-
wrapping lifetime issues; these remain open and require their own repairs/tests.
Final combined source, consumers and general resources remain unqualified here.
