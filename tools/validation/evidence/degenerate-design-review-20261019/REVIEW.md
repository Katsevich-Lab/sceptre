# Independently confirmed G1 blocker: non-estimable treatment gives discoveries

W03 independently read W02's actual frozen-source4d6d86d diagnostic receipt,
preserved here with SHA256 in `review-receipt.log`. No new inference fit or
simulation was performed during this review. The reproducible reader also
recomputes design rank and projection directly from the saved covariates.

The source diagnostic uses160cells, four independently generated Poisson(5)
response genes, low-MOI two-guide union versus complement controls, and includes
the exact binary treatment indicator in the covariate design. All genes are
true nulls; all four pass the default nonzero QC thresholds. Under BOTH the
public permutation AND CRT workflow, three genes have NaN observed statistics,
p=.004987531 atstage3, and BH-significant results. The fourth has a tiny finite
score-.00438,p1. Variances are signed floating cancellation aroundzero. This is
one diagnostic dataset with four endpoints, not a precise calibration estimate.

Independent computation gives rank(Z)=rank([Z,X])=2 and maximum absolute residual
of projecting X onZ=1.09e-14. Treatment lies in the nuisance-design column space.
For positive GLM weights, efficient information is
`X'WX - X'WZ (Z'WZ)^(-1) Z'WX = 0`.
The tested coefficient cannot be separated from a nuisance coefficient. A
rounded variance of either sign and an arbitrary tiny residual numerator must
not turn that absence of information into a discovery. The probability-domain
result guard misses this failure because the invalid p-values lie in[0,1].

**Admit a separate explicit identifiability precondition correction**, owned by
W02 and reviewed/integrated by W10. Check the treatment contrast on the exact
post-QC fitting/control view per guide group, including the NT restriction when
applicable, before inference/draws. Error clearly with the group/covariate
aliasing/no residual treatment variation. Make any numerical rank tolerance
explicit: rank decisions also classify near aliases and should not be called
symbolic exact equality. Do not clip variance toepsilon, ridge, modify counts,
or perturb treatment/design to force a test. Show ordinary identified-case
results and post-run RNG unchanged, and cover direct/distributed entry points.

A separate nonfinite **observed**-score rejection is also scientifically
justified and can protect other degenerate observed designs. This does not
justify blanket replacement/dropping of nonfinite **null resample** scores;
CRT empty/full resamples have distinct inherited semantics. No empirical-tail
counting change is proposed under this guard. A rejected non-estimable design
has no valid p-value; do not replace it with a small p or treat it as a completed
analysis. A clear error is appropriate for this unsupported inference request.

The actual corrected installed candidate must rerun this public diagnostic and
the positive identified controls. Feature scientific acceptance is held pending
that terminal evidence. W03 has not modified package scientific code or claimed
the guard implemented. The original issue237 all-ties cap is a separate already
justified numerical correction and does not fix this NaN-score failure.
