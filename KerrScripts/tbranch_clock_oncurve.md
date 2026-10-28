# t-branch clock: ENTIRELY ON-CURVE (the frame-dragging cover was spurious)

## Correction
The paper (old eq:clock-t) gave dt/dr = ρ_t/√R6 + c_β√(2Mr/(r²+a²)) with a
**double-cover** term (frame-dragging cover). **THIS IS WRONG.**
The t-branch clock lies **entirely on the genus-2 curve** R6, like the τ-branch.

## Key identity (exact, SymPy)
    (r-2M)Q2 + DE(J(r-2M)+2Ma)² = E² r³ Δ
with DE=(E²-1)r+2M and E²-f = DE/r.  Hence the Randers term
    n_t α = E² r³/(f √R6)   EXACTLY (no cover).

## Corrected closed form
    dt/dr = ρ_t/√R6,   ρ_t = [E² r³ - (2Ma/r) K]/f = P_3(r) + R_Δ(r)/Δ
- the pole at the ergosphere r=2M is **removable** (the numerator vanishes there)
- third-kind poles ONLY at the horizons (Δ=0)

**Second-kind part (polynomial):**
    P_3 = E²r³ + 2E²M r² + 4E²M² r + 8E²M³ - 2(E²-1)JMa
**Clock vector b** = (8E²M³-2(E²-1)JMa, 4E²M², 2E²M, E², 0)
  (the t analogue of the τ-branch b=(0,0,-2M,1,0))
**Third-kind horizon part** (numerator of degree 1):
    R_Δ = -2M(2E²JMar - E²Ja³ - 8E²M³r + 4E²M²a² + 2E²Ma²r + Ja³ - 2Ma²r)

## Verification
- Identity G=E²r³Δ: exact (SymPy)
- On-curve dt/dr against the **Hamilton flow** (H of the t-branch): agreement 1e-10 at r=4, 6, 9, 11
  (the sign only reflects the ingoing orientation)
- Why the old "verified to 1e-14" with the cover was wrong: c_β√(2Mr/(r²+a²))
  is about 0.85 at r=6 (NOT small); if it were in dt/dr it would shift it by 0.85, but the
  on-curve form agrees with the flow to 1e-10. Different radicals (√R6 against the cover) cannot
  coincide with c_β≠0.

## Consequence (simplification)
The t-branch is now **fully parallel** to the τ-branch: both are on-curve, and δφ = ψ is
assembled with the W_kj of R6 plus (t-branch) a third-kind log-σ at the horizons. No
double cover. The τ-branch has a weight-2 dipole at ∞; the t-branch has a third-kind letter
at the horizons.

Scripts: `reproduce_reductions.py` (identity + P_3 + b + R_Δ),
`kerr_tbranch_separatrix_weierstrass.py` (Q2, Hamilton flow).
Paper: corrected eq:clock-t (main + PRD), tab:branch-adiab, ψ assembly.

## Closing the loop: t-branch ψ assembly (end to end)
Clock t = eta2 + eta3:  eta2 = sum_j b_j U_j (second kind), eta3 = int R_Delta/(Delta sqrt R6) (third kind, horizon).
- Assembly eq:psi-split (second-kind part): I_direct = I_assembled = 10.31674453,
  diff **8.85e-13** (W_kj reorganisation of R6 with the vector b^t).
- Full clock: int dEF*t = I_assembled + int dEF*eta3 = 10.41916368, diff **8.85e-13**.
- horizon piece int dEF*eta3 = 0.1024 (weight-2 log-sigma at r_pm).
The orbit at J=2.5 is scattering (no turning point), integrated down to rp+0.6 (near the horizon).
Script: `kerr_tbranch_psi_assembly.py`. Paper: ψ assembly + tab:script-map.
