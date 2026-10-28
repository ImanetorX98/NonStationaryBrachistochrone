# Adiabatic correction in Vaidya — set-up (frozen = genus-2 Schwarzschild, correction in dM/dv)

Programme: closed forms of the t/τ brachistochrones in **Vaidya** (dynamical
Schwarzschild, with a radiating mass M(v)), by analogy with Thakurta-Kerr.

## 1. Metric and brachistochrone
Ingoing Vaidya:  ds² = −(1−2M(v)/r)dv² + 2 dv dr + r²dΩ².
Spherical, NO spin, NO ergosphere; DYNAMICAL horizon. Null coordinate v (advanced
time). Branches: **v-branch** (minimum advanced time v, the analogue of t) and **τ-branch**
(proper time).

## 2. FROZEN family = Schwarzschild (the a=0 limit of Thakurta-Kerr)
Setting a=0 in the genus-2 forms of Thakurta-Kerr:
- τ: F_τ = J√Emu/√S,  Emu=(E²−1)r+2M,
  S = r(r−2M)Emu(r²(r−2M)−J²Emu)  [SEXTIC, genus 2].
- **NOT elliptic**: 6 distinct roots (discriminant ≠ 0). The CONSTRAINED
  brachistochrone has genus 2 even in Schwarzschild (the elliptic curves in the literature belong to the
  GEODESICS). It is elliptic only on the separatrices |J|=Jc.
- The structure simplifies compared with Kerr: horizons r_+=2M=r_e, r_-=0; no frame dragging.

## 3. Adiabatic expansion (slow M(v))
Slow parameter: Ṁ = dM/dv. Frozen-coefficient WKB:
  φ(r) = φ_0(r; M(v(r))) + Ṁ · δφ + O(Ṁ²),
  δφ = ∫ ∂_M F · v(r) dr · (factor),  v(r) = advanced time along the orbit (clock).
STRUCTURE IDENTICAL to Thakurta-Kerr, with ∂_M F ↔ ∂_E F, v(r) ↔ η(r), Ṁ ↔ A'/A.

## 4. Reduction of ∂_M F — VERIFIED ✓ (`vaidya_dMF_reduction.py`)
∂_M F = N_M/S^(3/2), N_M = S ∂_M K − K ∂_M S/2 is POLYNOMIAL (degree 6),
N_M = J r Emu·(polynomial of degree 4). Reduction:
  ∂_M F = d(𝒜/√S) + Σ_{k=0}^4 c_k^M r^k/√S,
verified dM F(direct) = d(𝒜/√S) + M_poly/√S to **1e-15**.
c_k^M at (M=1,E=7/5,J=5/2) = [−0.8715,−0.3597,1.0354,0.2129,−0.1894].
(Check: N_E=EJr⁴(r−2M)²Emu = the a=0 limit of Thakurta-Kerr ✓.)

## 5. What follows (the same pipeline as for Thakurta-Kerr, verified there)
- clock v(r): advanced time along the frozen orbit (the analogue of η=t/τ). To be derived
  explicitly (dv/dr for the Vaidya brachistochrone).
- decomposition: δφ = ½ · Ṁ · Σ_{k<j} Q_kj^M W_kj + (clock × abelian products) +
  (elementary logs), Q_kj^M = c_k^M b_j − c_j^M b_k, W_kj = genus-2 polylogs (THE SAME
  as for the curve S; they depend only on the curve). ALGEBRAIC coefficients.
- The W_kj are already those of Thakurta-Kerr at a=0 (the same genus-2 curve S).

## 6. Differences from Thakurta-Kerr
- NO conformal transfer theorem (M(v) is not a conformal factor): the frozen metric
  is Schwarzschild(M), not Kerr(E_eff). But the adiabatic structure is the same.
- Null coordinate v → different causality; Vaidya remains a BLACK HOLE (constraint set by the author).
- The frozen metric is structurally simpler (no spin) but has the SAME genus (2).

Script: `vaidya_dMF_reduction.py`. References: `../ThakurtaMetric/
phi_adiabatic_closed_form.md`, `../KerrMetric/kerr_psi_explicit_verified.py`.

## 7. Clock v(r) — VERIFIED ✓ (`vaidya_clock.py`) + a NEW result
Advanced time (M=M(v)): v=t+r_*, hence
  dv/dr = (E/f) dτ/dr + 1/f = E r³/√S + r/(r−2M)   (symbolic check: diff=0).
  dτ/dr = r²(r−2M)/√S (a=0 limit, τ-branch).
v(r) = E R_3 + r + 2M ln(r−2M),  R_3=∫r³/√S (second kind).
The tortoise term r/(r−2M) → **horizon logarithm** 2M ln(r−2M).

**δφ_V = Ṁ ∫∂_M F·v dr** splits into THREE pieces (verified to 1.8e-15):
1. genus-2 polylog (from E R_3) — as in Thakurta-Kerr.
2. elementary (from r) — second kind + polynomial.
3. **∫∂_M F·ln(r−2M) dr — HORIZON DILOGARITHM, NEW.**

RESULT: the advanced time of Vaidya injects a horizon log (v~ln(r−2M) at
r_+) → a hyperelliptic dilogarithm with a dlog "letter" AT R_+ (the horizon),
distinct from the Thakurta-Kerr dilog (which came from the third-kind dipole at ∞).
Physics: the mass variation is weighted by the log divergence of the advanced time
at the horizon. Algebraic coefficients (c_k^M); the W_kj are the same (curve S at a=0).

## 8. EXPLICIT CLOSED FORM of δφ_V — VERIFIED ✓ (reusable)
`vaidya_explicit.py`. Verified to 9.4e-14 (coefficients GIVEN, no fit):

  δφ_V/Ṁ = A_M(r) v(r) − E ∫A_M r³/√S dr − ∫A_M r/(r−2M) dr

INGREDIENTS (all explicit):
- A_M = ∂_Mφ_0 = 𝒜^M/√S + Σ_{k=0}^4 c_k^M U_k  (reduction of ∂_M F, ALGEBRAIC c_k^M;
  at E=7/5: c_k^M=[−0.8715,−0.3597,1.0354,0.2129,−0.1894]).
- v(r) = E U_3 + (r−r0) + 2M ln((r−2M)/(r0−2M))   (advanced-time clock).
- U_k = ∫r^k/√S dr  (abelian integrals of the first/second kind, on the Schwarzschild curve S at a=0).

DECOMPOSITION into NAMED functions (expanding the two integrals):
  ∫A_M r³/√S = [∫𝒜^M r³/S elementary] + Σc_k^M(½U_k U_3 + ½W_{3k}) − const·U_3
  ∫A_M r/(r−2M) = [∫𝒜^M r/((r−2M)√S) = L_{2M}, third kind at the horizon]
                 + Σc_k^M·2M·D_k + Σc_k^M(rU_k−U_{k+1}) + elementary
- W_{jk}=∫(U_j dU_k−U_k dU_j): genus-2 polylog at ∞ (as in Thakurta-Kerr).
- **D_k=∫U_k/(r−2M) dr: genus-2 dilog at the HORIZON (NEW in Vaidya).**
- L_{2M}=∫𝒜^M r/((r−2M)√S): third kind at the horizon (θ-log, weight 1).

REUSE: given the curve S (= Schwarzschild, a=0), the c_k^M (algebraic) and the named
functions {U_k, W_jk, D_k, L_{2M}}, δφ_V can be evaluated for any slow M(v).
Transcendence: W_jk (polylog at ∞) + D_k (horizon dilog); algebraic coefficients.
Difference from Thakurta-Kerr: the dilog sits at the HORIZON (tortoise term of the advanced time),
not at ∞ (third-kind dipole).
