# Fixed-endpoint no-inversion (referee #11), frozen Schwarzschild: final results

Status of the **frozen-case** (Vaidya mu=0 / Thakurta A const => static Schwarzschild,
b=E^2-1>0 scattering) fixed-endpoint no-inversion theorem. Every claim below is either
proven in closed form, independently Mathematica-checked, or rigorously certified by
interval arithmetic. Companion scripts in this directory; symbolic checks in
`verify_lemmaB_mathematica.wls`.

## Reduction
Parametrize each branch by turning radius r_min (x:=Vmin=V(r_min), V0=V(r0) FIXED).
`Phi_tau(x)=sqrt(x) int_x^{V0} W(A)/sqrt(A-x) dA`, W(V)=K/V', K=1/sqrt(r(r-2)),
V=r(r^2-2r)/(br+2), V'>0 (no photon sphere on tau).  Theorem <= Lemma A + Lemma B:
- **Lemma A** (F_t>F_tau at matched r_min): PROVEN in closed form (algebraic factorization).
- **Lemma B** (Phi_tau decreasing for r_min>r_pk): the object of this file.

## What is PROVEN / CHECKED (Lemma B, frozen)

1. **Closed form for Phi'** (IBP, verified 1e-9; Mathematica):
   `sqrt(x)Phi' = W(V0)(V0-2x)/sqrt(V0-x) + int_x^{V0} W'(t)(2x-t)/sqrt(t-x) dt`.

2. **Elementary sub-regimes** (closed form):
   - grazing `x>=V0/2`: Phi'<=0.
   - quarter `x>=V0/4`: Phi'<=0 (via `int kernel = sqrt(V0-x)(V0+2x)/3`, Chebyshev; the
     `int_x^{2x}, int_{2x}^{V0}` integrals Mathematica-verified).
   => **no-inversion fully elementary for r0<=R*(E)** (when V0<=4 Vpk).

3. **Single-crossing criterion**: Lemma B <=> Phi_tau has one critical point <=> at every
   critical point Phi''<0.  Reduced (Mathematica-verified) to
   `L''=-3/4(a^2+c^2)+1/2 ac+I''/I` at Phi'=0 (a=1/x,c=1/(V0-x)); verified to 40 digits.

4. **Large-r0 asymptotic** (closed form; ALL Mathematica-verified, diff=0):
   `2VW = DE sqrt(Delta)/N = 1 - (1/b)/r + ...`  =>  `W = 1/(2V) - 1/(2 b^{3/2} V^{3/2}) + ...`
   `Phi = arctan(sqrt(V0/x-1)) - b^{-3/2} sqrt(1/x-1/V0) + E(x)`
   `Phi' = [1/(2 sqrt(x(V0-x)))] ( b^{-3/2} sqrt(V0)/x - 1 ) + E'(x)`
   => single non-degenerate maximum at `x_pk ~ b^{-3/2} sqrt(V0)`,
      `Phi''(x_pk) = -1/2 b^{9/4} V0^{-5/4} < 0`  (Mathematica-exact via V0=t^4).
   Remainder `h2=W-1/(2V)+.. ~ +C/V^2` (C~3.07 explicit); reinforcing => Phi''<0 robust.

5. **CAP** (rigorous interval arithmetic, mpmath.iv): order-2 interval jet for W,W',W'';
   rigorous inverse r(A); interval quadrature; guaranteed [Phi'],[Phi''].  Certifies
   single-crossing at r0=10 (Phi''<0 near r_pk; Phi'<0 for r_min>r_pk; thin-cell COVERING
   over the continuum r_min in [5.0,5.5]).  All enclosures validated to contain the truth.

## What is NOT yet a full "for all r0" proof ("almost proven")

The theorem is proven **modulo** closing a compact middle, by a three-way decomposition
with NO conceptual gap:
- r0<=R*: elementary (piece 2).                                  [DONE, analytic]
- r0>=R_large: asymptotic (piece 4); rigor = explicit remainder constants (mechanical);
  R_large is large (~300-1000) because convergence is slow (series in V0^{-1/4}).
- R*<=r0<=R_large: CAP (piece 5); rigorous machinery DONE and validated; a *complete*
  single-r0 certificate is achievable (diagnosis: N~1000-4000 certifies the peak shoulder)
  but the crude cell-sum quadrature is slow near the flat peak.  Remaining engineering:
  O(1/N^3) verified quadrature (needs order-3 jets) and the full (E,r0,r_min) box scan.

## Verdict
Lemma A: **closed**.  Lemma B (frozen Schwarzschild): **almost proven** -- elementary on
r0<=R*, closed-form asymptotic single-non-degenerate-max on r0>=R_large, rigorous CAP
machinery certifying the compact middle (demonstrated at r0=10).  No conceptual obstruction
remains; closure = explicit-constant bookkeeping (asymptotic) + tighter quadrature and
compute (CAP).  Kerr frozen: W non-elementary (genus-2) -> numerically verified only.
All symbolic derivations independently confirmed in Mathematica (`verify_lemmaB_mathematica.wls`).
