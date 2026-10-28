# Adiabatic φ(r,A) — semi-analytic closed form and polylog structure

Consolidated note of the programme: **closed form of the brachistochrones φ(r)
with a slowly varying conformal factor A(η)** (Thakurta–Kerr → Kerr through
E_eff=Ê/A). Branches t (coordinate time) and τ (proper time). All the results
were validated numerically.

Reference parameters: M=1, a=9/10 (or 0.9), E=Ê, J. Genus-2 curve of the branch:
y²=R₆ (t), y²=S_τ (τ).

---

## 1. The final form (goal reached at the level of structure)

```
φ(r,A) = φ₀(r; Ê/A)  +  (A'/A) [ Closed(r) + ψ(r) ]  +  O((A'/A)²)

Closed(r) = −½ Ê · ∂_E φ₀(r) · η(r)          [CLOSED]
ψ(r)      = ½ Ê ( ρ − ρ̃ )                     [genus-2 polylog]
            ρ = ∫ ∂_Eφ₀ · h dr,  ρ̃ = ∫ η · ∂_E F dr,  h = dη/dr
```

η(r) is the **clock** accumulated along the frozen orbit: coordinate time t(r)
(t branch) or proper time τ(r)=∫L_τ dr (τ branch). The τ branch only has, in addition, the
prefactor A⁻² on the timing, NOT on the shape; the shape of φ is common to the two branches.

**How φ(r) is obtained:** each piece is a special function on the same genus-2
curve → evaluate and add. See the table in §2.

---

## 2. The pieces of φ(r) — all identified

| piece | explicit form | class | status |
|-------|-----------------|--------|-------|
| φ₀(r) | ∫F dr, F=dφ/dr (K/√R₆ or K_τ/√S) | first + third kind (Kleinian σ,℘,ζ) | ✓ |
| ∂_Eφ₀(r) | A(r)/√R + Σ_{k=0}^4 c_k ∫r^k/√R (no third kind) | second kind | ✓ err=0 |
| η(r) | ∫(dη/dr)dr, clock (see §3) | second kind [τ] / second kind + cover [t] (superseded for t: see the note in §3) | ✓ |
| Closed | −½Ê ∂_Eφ₀ · η (product of the two above) | closed | ✓ |
| ψ(r) | ψ_ζ + ψ_ab + ψ_Li (THREE pieces, see §5) | weight-1 ζ + TWO weight-2 pieces | structure ✓ |

Closed φ(r) = [Kleinian pieces] + (A'/A)[Kleinian + one genus-2 polylog].

---

## 3. The clock differentials dη/dr — made EXPLICIT and validated

> **Historical section (note added in translation, October 2026).** Everything said
> below about a frame-dragging "cover" for the t branch, including the "Physics"
> paragraph and the t/τ "polylog dichotomy" proposed in §7 and §8, was later shown
> to be wrong: the t clock lies entirely on the curve. See
> `../../KerrScripts/tbranch_clock_oncurve.md`. The text is kept as it was written.

Prerequisite for Closed: is dη/dr on the curve? (test: (dη/dr)²·[sextic] is a square).

### τ branch (proper time) — purely ON-CURVE
- β cancels in Q; `(dτ/dr)²·S_τ = [r²(r−2M)]²` is a PERFECT SQUARE.
- **dτ/dr = r²(r−2M)/√S_τ** — a RATIONAL differential on y²=S_τ.
- Validation: `max|√(Q/w) − r²(r−2M)/√S_τ| = 8.9e-16`.
- ω_a=∂_E F_τ dr (second kind on the curve) + ω_b=dτ (on the curve)
  ⇒ **ψ_τ = PURE genus-2 hyperelliptic polylog.**

### t branch (coordinate time) — ON-CURVE + frame-dragging cover
- β cancels in Q_t; `(Q_t/w)·R₆ = E²r⁶ = (Er³)²` → √(Q_t/w)=Er³/√R₆ on the curve.
- BUT the term B/f brings in β=√(2Mr/(r²+a²)) (Doran river velocity), NOT in √R₆:
  ```
  dt/dr = ρ_t(r)/√R₆  +  c_β(r)·√(2Mr/(r²+a²)),   c_β=(1−2Ma²/(rΔ))/f
  ```
- Validation: `max|dt/dr − ρ_t/√R₆ − c_β·β| = 2.5e-14`.
  ⇒ **ψ_t = [genus-2 polylog on R₆] + [remainder on the frame-dragging cover β].**

  [Note added in translation, October 2026: this decomposition was later shown to be
  wrong, the t clock lies entirely on the curve; see `../../KerrScripts/tbranch_clock_oncurve.md`.]

### Physics
Proper time is frame-independent → a clean genus-2 curve.
Coordinate time drags β along (frame dragging) → an extra square-root cover.
**The t/τ polylog dichotomy = the frame-dragging dichotomy.**

---

## 4. Validation of the hybrid form (φ_hybrid == φ_full)

Hybrid = closed (analytic) pieces + numerical ψ. Comparison with the full adiabatic solution
(non-autonomous coefficient −Ê∫∂_E F·η dr). Identity by parts; residual = trapezoidal rule.

| A'/A | t branch (η=t) | τ branch (η=τ) |
|------|--------------|--------------|
| 0.005 | 1.4e-06 | 2.5e-08 |
| 0.02  | 5.6e-06 | 1.0e-07 |
| 0.06  | 1.7e-05 | 3.0e-07 |

Scripts: `kerr_adiabatic_phi_hybrid.py` (t), `kerr_adiabatic_phi_hybrid_tau.py` (τ).

---

## 5. The nature of ψ — decomposition by kind (FINAL, from the residues)

An honest history (two corrections): the paper said "polylog"; then I said "pure ζ"
(WRONG); the decomposition of the residues gives the precise answer.

**Residues computed** (`/tmp` on-curve work, rational parameters):
- `res_∞(ω_a=∂_E F) = 0`  → ω_a is purely of the **first + second kind** (no third kind; also β_±=δ_±=0
  at the horizons). ✓
- `res_∞(ω_b=dη) = 1.063` ≠ 0 → dη has a **third-kind dipole at ∞**
  (physics: dη~dr/r at large r, η~log r; time grows logarithmically).

**ψ = ½Ê ∫(A dB − B dA)**, A=∂_Eφ₀, B=η. EXACT ANALYTIC decomposition (no fit,
`KerrMetric/kerr_holo_component_check.sage`) gives THREE pieces:
```
ψ = ψ_ζ  (weight 1)  Kleinian ζ,σ              terms (2nd kind)×(1st/2nd kind)
  + ψ_ab (weight 2)  holomorphic×holomorphic   det·∫(u1 du2 − u2 du1)   [Beilinson]
  + ψ_Li (weight 2)  third kind / dilog         ½Êρ₀[∂_Eφ₀ L − 2𝓛₂]
```
- **c_k** (second-kind reduction of ∂_E F): EXACT, rational in E (at E=7/5:
  [-0.531,1.979,-0.812,-0.360,0.189]). Analytic, SymPy.
- **ρ₀ = M/(E²−1)^(3/2)** EXACT (third-kind residue of the clock; at these parameters 1.0631).
- **ψ_ab**: coefficient = det(holomorphic components) = b^A_0 b^B_1 − b^A_1 b^B_0 =
  −4.27−7.94i ≠ 0 (|det|/|bA||bB|=0.80). Holomorphic components from exact a-periods.
  NEW for genus ≥ 2: in genus 1 there is only one holomorphic differential → ψ_ab is ABSENT (the elliptic zζ−2logσ
  is clean). The FIT hid it; exact algebra reveals it.
- **ψ_Li**: L=log(θ-ratio) [Fay's theorem, NOT a fit], analytic ρ₀, 𝓛₂=∫L ∂_E F dr
  an explicit hyperelliptic dilog (a single integral of θ, like Li₂).
- HISTORY of the corrections: paper "polylog" → "pure ζ" (WRONG) → "ζ+dilog" (2 pieces,
  INCOMPLETE) → **"ζ + holomorphic×holomorphic + dilog" (3 pieces, from exact algebra)**.
- ψ_ζ (weight 1) still has to be reduced explicitly with Baker's ℘ identities.

## 5bis. Quasi-periods η — COMPUTED and VALIDATED (the wall of §7-iii is broken)

`KerrMetric/kerr_quasiperiods_bel.sage`. The second kind is now computed with Sage:
- **ODD** model (quintic through x=1/s, one point at ∞): clean second-kind differentials.
  The even model (degree 6, two points at ∞) gives x²dx/y of the third kind → wrong η.
- Canonical **Baker–Enolski–Leykin** second-kind differentials:
  `dr_1=(λ3 s+2λ4 s²+3λ5 s³)ds/4y`, `dr_2=λ5 s²ds/4y`.
- Periods through `matrix_of_integral_values` (integers: no Singular crash).
  Pairing du=[1,s]↔dr=[dr1,dr2].
- **VALIDATION**: κ=η ω⁻¹ **symmetric to 1.4e-12**; generalised Legendre relation
  `ω'ηᵀ−ωη'ᵀ = −iπ·I` (canonical). ⇒ η is correct.

The ingredients of the closed form of ψ are NOW all available: u(r) (forward Abel map), θ+derivatives
(abelfunctions), κ/η (here). NEXT (mechanical): σ=exp(−½uᵀκu)θ, ζ_i=∂_i log σ,
ψ(r) = a closed combination of ζ(u(r)); validate against the numerical ψ.

**PAPER TO BE CORRECTED:** `sec:adiabatic` says ψ is a "genus-two polylogarithm /
irreducible" — imprecise. ψ closes in Kleinian ζ,σ (second kind), it is not a polylog.

---

## 6. Genus-2 engine (level 3) — abelfunctions in Sage

- `abelfunctions` 0.2.0 compiled in SageMath 10.9 (forced x86_64 build, API patch
  `is_LaurentSeries`; see the memory note `abelfunctions-sage-install.md`).
- It exposes RiemannTheta (with derivatives), RiemannSurface, AbelMap, RiemannConstantVector.
- Bug: the RiemannSurface pipeline (holomorphic_differentials→Singular integralbasis) crashes on
  non-monic leading terms/large coefficients. WORKED AROUND: τ, A|B from the **Sage** RiemannSurface;
  forward Abel map u(r) by direct integration; θ+derivatives from abelfunctions.
- Checks: A⁻¹B−τ=1e-16; u(r₀)=0; ∇θ(0)=1.6e-16 (θ even); symmetric non-zero Hessian.
- Script: `KerrMetric/kerr_psi_forward_abel.sage`.

---

## 7. In the paper (main.tex + main_prd_revtex.tex)

- Subsection `sec:adiabatic` "Semi-analytic first-order adiabatic orbit shape":
  hybrid equation, second-kind reduction, irreducible ψ, validation table for t+τ.
- Honest framing: ψ is NOT reducible to depth 1 (shuffle) but it is a genus-2 polylog.
- Acknowledgements + software bibliography: SageMath, abelfunctions (with link),
  NumPy/SciPy/SymPy/Matplotlib (with DOIs).

**TO BE ADDED (proposed; withdrawn, see the note at the start of §3):** the t/τ polylog dichotomy (§3) — a clean result:
proper time = clean curve, coordinate time = frame-dragging cover.

---

## 8. What remains

1. (optional) Figure of the x-y curves φ(r,A) assembled from the analytic pieces, t and τ.
2. Genus-1 separatrix: ψ closes in the tabulated elliptic dilogarithm (closed win).
3. Named genus-2 Kronecker–Eisenstein kernel for a generic ψ (research grade).
4. Write the t/τ polylog dichotomy into the paper. (Withdrawn: see the note at the start of §3.)

Related files: `../progress.md` (§8,§9), `../KerrMetric/doranTau.md`,
`../KerrMetric/doranT.md`, scripts in `../KerrMetric/` and `.`.
