# Brachistochrones in Thakurta (conformal Schwarzschild) — results

**Script**: `thakurta_brachistochrone_sympy.py` (all checks pass).
Context: `../FLRWmetric/FLRWresults.md` (conformal move),
`../SdSMetric/SdSresults.md` (static limit with Λ),
`../VaidyaMetric/VaidyaResults.md` (PMP Hamiltonians).

```
g = a(η)² · ĝ_Schw :   ds² = a(η)² [ −f dη² + dr²/f + r²dΩ² ] ,  f = 1−2M/r
```

The exact marriage of the two building blocks: the conformal move of FLRW on a static
Schwarzschild base. As a physical model of a cosmological black hole it is debated
(radial heat flux, PBH debate); as a geometric laboratory it is
the minimal case with **non-trivial optical geometry AND a non-trivial clock**.

## R1. Symmetry
`£_{∂η} g = (2a′/a) g`: `∂_η` is a **conformal Killing field for every a(η)**
(verified); genuine ⇔ `a′ = 0`. No stationarity: Perlick does not
apply directly.

## R2. Kinematics of the conformal rail `Ê = −u_η`
```
dτ/dη = a²f/Ê ,    γ = Ê/(a√f) ,    v² = 1 − a²f/Ê²
```
**Freezing surface** `a(η)²f(r) = Ê²` — a curve in (η, r) that
unifies the two known barriers:
- FLRW limit (f=1): the cosmic instant `a = Ê`;
- static limit (a=1): the radius `f = E²` (the w=0 barrier of SdS/Schwarzschild).

With growing `a` the accessible region `f < Ê²/a²` **shrinks towards
the horizon**: the expansion pushes the rail particles onto the black hole.

## R3. Indicatrix and Hamiltonians (exact)
An ellipse without wind (static base, `′ = d/dη`):
```
r′ = f v cosθ ,    φ′ = (√f·v/r) sinθ
H_η = v√f·√( f p_r² + J²/r² ) − 1
H_τ = v√f·√( f p_r² + J²/r² ) − a²f/Ê
```
`p_φ = J` conserved; `dH/dη ∝ a′` (conformal forcing); transversality
`H = 0` at arrival; static ⇒ H conserved ≡ 0.

## R4. Fermat-type reduction — where one goes BEYOND Perlick
Along every rail worldline: `dη = dl_opt / v` with `dl_opt` the arc length
of the **Schwarzschild optical metric** (`dr²/f² + r²dφ²/f`). The η
branch is Fermat with index
```
n(η, r) = 1/v = Ê / √(Ê² − a(η)² f(r))
```
- `a = 1`: `n = E/√(E²−f)` — **Perlick 1991 recovered** exactly;
- generic `a(η)`: `n` depends **jointly** on η and r ⇒ no
  global clock (the λ(η) trick of FLRW requires f=1, homogeneity):
  a non-separable Fermat principle — beyond Perlick in an essential way;
- τ branch: the cost `a²f/Ê` is weighted in space AND time: never reducible.

## R5. Static validation (a=1)
Hamilton flow (Schwarzschild coordinates, not EF) against the closed forms:
`dφ/dr` to 10 digits at r = 4, 6, 9 for **both** branches
(`𝒦_τ = J`, `𝒦_η = fJ/E`).

## R6. Capture by expansion (conformal de Sitter, a = −1/(Hη), H=0.02)
An orbit launched **outwards to escape** from r=4 (J=1.3, Ê=1.2, a(η₀)=1):
- it rises up to r = 8.74 while the freezing surface comes down
  (∞ → 8.74): **contact** at η = −36.59 with `v → 0`;
- after the contact the indicatrix degenerates to a point: the particle
  **rides the surface**, which drags it to `r_freeze → 2M`
  (3.21 → 2.004 as η → 0).

**Capture by expansion**: a phenomenon absent both in FLRW (no
horizon) and in Schwarzschild (no expansion) — the first genuinely
"cosmological black hole" result of the programme.

> **[PARTLY SUPERSEDED, September 2026]** The contact with the freezing surface stands; the
> continuation "rides the surface" is not a proven dynamics. At freezing the indicatrix
> collapses to a point and the control problem stops being well posed: the papers treat
> freezing as a boundary of the domain (Paper I, Lemma I.A(ii); Paper II, `lem:compact`).

## R7. ORIGINAL Thakurta: conformal Kerr (`thakurta_kerr_sympy.py`)
Thakurta's paper (1981) is "Kerr in an expanding universe":
`g = A(η)²·g_Kerr` — the minimal **rotating** cosmological black hole. Equatorial,
rail `Ê = −u_η` (`v̄² = 1−A²f/Ê²`, `P̄ = P + A²(2Ms/r)²/Ê²`):

**Indicatrix = ellipse with an angular wind** (exact):
```
centre: φ′₀ = (2Ms/r)v̄²/P̄  (frame dragging) ,  r′₀ = 0
R² = f v̄² + P̄φ′₀²  =  Δ·v̄²/P̄        [identity: fP + (2Ms/r)² = Δ]
```

**Hamiltonians in closed form**:
```
H_η = p_φφ′₀ + R√( (Δ/r²)p_r² + p_φ²/P̄ ) − 1
H_τ = p̃_φφ′₀ + R√( (Δ/r²)p_r² + p̃_φ²/P̄ ) − A²f/Ê ,
      p̃_φ = p_φ − 2MsA²/(rÊ)     (CONFORMAL gravitomagnetic shift)
```

Validations: s→0 = Thakurta-Schwarzschild (numerical, 1e-12); A=1 =
equatorial Kerr closed forms doranTau/doranT (`𝒦_τ = J`,
`𝒦_t = (fJ+2Ms/r)/E`) to 10 digits at r=4,6,9 (spin 0.9, E=1.2, J=1.3).

**R7a — the two critical surfaces separate sharply:**
- `Δ = 0` (horizon `r_+ = M+√(M²−s²)`): a **conformally
  invariant** degeneration — it depends neither on A(η) nor on Ê (verified: the root of
  R² is identical to 1e-13 for every A). The indicatrix stays regular
  across the ergosphere (at f=0: R² = v̄²Δ/P̄ > 0), like the constrained
  worldline in Doran coordinates;
- `v̄ = 0` (conformal freezing `A²f = Ê²`): it **breathes** — it comes down
  from infinity when A exceeds Ê (A=1.25→r=25.5; A=3→r=2.38) and
  compresses the accessible region into the horizon–freezing shell:
  the rotating version of the capture by expansion.

> **[PARTLY SUPERSEDED, September 2026]** The algebraic regularity of the indicatrix beyond the stationary
> surface is not physical admissibility: the control domain adopted in the papers is exterior to
> r = 2M on the equator (`lem:compact`), and what lies beyond is analytic continuation (type (A)).

## R8. Conformal-time-dependent optical metric for Thakurta-Kerr

> **[Note, September 2026]** The same distinction as in R7a applies: the formulas stay regular beyond
> r = 2M, but the physical domain is exterior. With a varying A the optical metric depends on η
> (∂_η a = 2Ê²AA′h₀/(Ê² − A²f₀)², `eq:aT-running`), and it does not reduce to an autonomous metric.
(`thakurta_kerr_optical_sympy.py`)

Solving the indicatrix for `dη` at a given spatial displacement — **non-autonomous Randers
in closed form**:

```
dη = n(η,r) · α_K + β_K

α_K² = r²dr²/(fΔ) + Δdφ²/f²          (Riemannian)
β_K  = −(2Ms/(rf)) dφ                 (gravitomagnetic one-form)
n(η,r) = 1/v̄ = Ê/√(Ê² − A(η)²f)     (index)
```

**Theorem (verified, residuals ~1e-14):** `α_K` and `β_K` are *exactly*
the Randers data of the **null** Fermat principle of Kerr (Gibbons–Werner). All the
physics of mass, rail and expansion sits in the scalar index `n(η,r)`,
which multiplies only the Riemannian part; the frame-dragging one-form is
**rigid**: conformally invariant and independent of Ê.

Limits (verified): `Ê→∞` ⇒ n=1 (null optics of Kerr); `A=1` ⇒
`n = E/√w` (stationary Perlick on Kerr; the Euler–Lagrange equation of F reproduces doranT to 10
digits); `s=0` ⇒ Thakurta-Schwarzschild; `M=s=0` ⇒ `n = 1/v` FLRW.

Singularities: `f=0` (ergosphere) — α, β diverge: the *optical reduction*
dies at the ergosphere (an intrinsic 1/F, as in static Kerr) while
the Hamiltonian/indicatrix stays regular down to the horizon (R7a);
`v̄=0` — n→∞ at the conformal freezing (a moving surface).

It is the "non-autonomous Finsler" anticipated in `nonStationaryBrachi.md` §1.2,
now explicit: **Zermelo with a time-dependent landscape, in which only the metric
term breathes, not the wind**.

### R8a. The three branches from the same ellipse (different costs, same indicatrix)

| branch | cost L₀ | Hamiltonian | optical metric |
|---|---|---|---|
| η | 1 | H_η (R7) | `dη = n·α_K + β_K` (Randers) |
| cosmic t | A(η) | `H_t = H_η form − A(η)` | **the same curves as η**: `t(η)=∫A dη` monotone ⇒ optical clock (degenerate η,t pair as in FLRW) |
| τ | `A²𝔉/Ê` (control-dependent) | H_τ with `p̃_φ` (R7) | `dτ = (A²f/(Êv̄))·α_K` — **PURELY Riemannian** |

Cancellation verified (1e-12 at random points of the ellipse):
substituting `dη = nα+β` into `dτ = (A²/Ê)(f dη + (2Ms/r)dφ)`, the term
`fβ` cancels the gravitomagnetic one exactly ⇒ the τ branch has no wind —
the conformal version of the "purely Riemannian τ branch" of doranTau, which therefore
survives the expansion. The wind is a phenomenon of the arrival-time
branches only.

## Figures
`make_figures_thakurta.py` → `Thakurtafigures/`: capture by expansion
(orbit against surface + accessible region), critical surfaces of
Thakurta-Kerr (rigid horizon against breathing freezing surface), residuals of the
flow against the Kerr closed forms at A=1 (1e-13). In `../PaperFigures/`:
`fig_indicatrices` (conceptual: the indicatrices in the three cases).

`make_figure_thakurta_branches.py` → `fig_thakurta_kerr_branches` (R8a as a
figure): (a) with |J|=1.3, arrival r=8: the τ branch BOUNCES in a mirror
pair (`Δφ = ±2.039`, periastra equal to **1e-11** — pure Riemannianity
is not imposed, it emerges from the flow under p_φ → −p_φ); the η branch SPIRALS
onto the horizon for both signs (self-tuned marginal crossing,
|J| < J_t^c ≈ 3.4) and **both spirals end up
prograde** (`Δφ = +3.45` and `+2.55`): frame dragging forces
co-rotation (the doranT mechanism). (b) the two conformal scalars `n` and `k_τ`
on the rigid α_K, with their asymptotes at freezing.

## R9. Inversion of the τ branch and trichotomy at the ergosphere
(`thakurta_kerr_plunge_inversion.py`, figure `fig_thakurta_kerr_plunge_map`)

Terminological NOTE: the trichotomy (doranTau) classifies the behaviour
at the light surface/ERGOSPHERE `r_e = 2M`, not at the horizon.
Capture (crossing + spiral onto r_+) happens ONLY on the
separatrix `J = J_c`, with measure zero; below J_c there is a cusp
reflection at r_e, not a plunge.

At constant A, `g = A²g_Kerr` is rescaled Kerr: the rail Ê ⇔ the Kerr problem
with `E_eff = Ê/A`, `J_eff = J/A` (T_τ̄ = A·T_τ ⇒ p_φ → A·p_φ).
**Analytic inversion radius** (rescaled trichotomy):

```
J > J_c(A):  smooth periastron, root of Δ(r) − (J/A)²[(Ê/A)² − f(r)] = 0
J < J_c(A) = s·A²/Ê:  CUSP reflection at the ergosphere, r_min = r_e = 2M
J = J_c(A):  crosses the ergosphere and spirals onto r_+ (capture, measure zero)
```

**Static colormap validation** (25×25 grid in (J, A), Hamilton flow
against the analytic result, cusp included): relative residual **max 7e-12**
— analytic contours indistinguishable from the numerical colormap; the line
`J_c(A)` is the capture separatrix. Verified `H_τ(p_r=0, r_inv)=0`.

**Dynamical (a = −1/(Hη), H=0.02, launch from r₀=8, A₀ ∈ [1, 1.3])**:
the non-adiabatic shift of the inversion `(r_min^num − r_qs)/r_qs` is
**always positive, +49% ÷ +120% (median +73%)**: the inversion is
**anticipated** relative to the quasi-static one — the approaching freezing
surface decelerates the particle globally and makes it turn much
earlier; no captured cell (the separatrix keeps measure ~zero on
this grid). Unlike Vaidya, adiabaticity fails at
first order even for small H: the drift `dH/dη ∝ A′` is cumulative
along the infall. Open question: does the dynamics widen the capture separatrix
at the ergosphere into a finite window (the analogue of Vaidya R8)? This needs
a zoom around `J_c(A)` with a cusp/crossing classification.

## R10. Plunge inversion t vs τ (Kerr protocol) and its conformal control
(`thakurta_kerr_plunge_t_tau.py`, figure `fig_thakurta_kerr_plunge_t_tau`;
complementary: `thakurta_kerr_inversion_t_tau.py`)

Protocol of `KerrMetric/evidence_t_plunges_deeper.txt`: the same initial
conditions `(r₀, p_r0, J)` in the flows of the QUADRATIC optical Hamiltonians
of the two branches, integration down to `p_r = 0`, comparison of `r_min` relative
to the singularity. The Thakurta-Kerr quadratics (from R8/R8a):

```
H_τ = (1/2k²)[(fΔ/r²)p_r² + (f²/Δ)J²]
H_t = (1/2n²)[(fΔ/r²)p_r² + (f²/Δ)X²],   X = J + 2Ms/(rf)
```

at A=1 they coincide EXACTLY with those of the Kerr project
(`c_p = FΔ·DE/(E²Σ³)` etc., verified term by term).

**Ratio of the Hamiltonians** (centrifugal potentials at p_r=0):

```
V_t/V_τ = ρ² ,   ρ = A²[ f + 2Ms/(rJ) ]/Ê
```

**The inversion plane is (a, J)** — black-hole spin against the z angular
momentum of the particle, as in the Kerr study; A is the parameter that moves
the locus.

- **Quantitative reproduction of the evidence** (A=1, a=0.4, r₀=10,
  p_r=−0.5): Δr(J=0.5) = **+0.0616** (evidence: +0.06), Δr(J=4.5) =
  **−0.417** (evidence: −0.42) — the conformal Hamiltonians reproduce the
  numbers of the independent Kerr code.
- **Colormap (J, a) at A=1**: algebraic r_min (H conserved), ODE to
  **7e-11**; exact inversion locus J_inv(a) from 0.32 (a→0) to 1.75
  (a=0.98), growing with the spin (frame dragging).
- **Corrected analytic formula**: the inversion condition is

  ```
  ρ(r*)² = H_t0 / H_τ0        [NOT ρ = 1]
  ```

  because at equal (r₀, p_r0, J) the two branches start with DIFFERENT
  Hamiltonian values (coefficients 1/n² against 1/k², X against J): the ratio
  of the potentials at the turn must be normalised to the launch ratio.
  In the far field `H_t0/H_τ0 ≈ ρ(r₀)²` ⇒ compact reading:
  **inversion where ρ(r*) = ρ(r₀)** (the profile ρ returns to its launch
  value). Verified against the exact locus: deviation **8.5e-13**.
  (The naive version ρ=1 assumed H_t0 = H_τ0 and was shifted.)
- **Numerical/analytic overlay**: colormap and contour Δr=0
  redone with PURE ODEs (28×28×2 flows); numerical inversion locus
  by bisection on the flows against the formula: point-by-point agreement,
  **max |ΔJ| = 2.4e-9** on a ∈ [0.24, 0.96] (for a ≲ 0.2 J_inv drops
  below the scan floor J=0.3: frame dragging vanishes and
  the inversion moves to J→0).
- **Conformal drift**: A grows ⇒ the inversion locus moves to higher
  J (J_inv(a=0.9): 1.58 → 1.68 → 1.93 for A = 1 → 1.1 → 1.2):
  the "τ deeper" region expands with the expansion; in the far
  field ρ → A²/Ê: for **A > √Ê** the inversion becomes generic
  (E_eff = Ê/A below ~1: time dilation wins, a TOV-like
  crossover, now driven by the scale factor).

### R10b. Theorem: the conformal inversion exists and is reachable
(`conformal_inversion_AJ.py`, `fig_thakurta_kerr_inversion_AJ`)

**CORRECTION** of an earlier WRONG note ("the inversion does not enter
the accessible window, 0/484 cells"). The count 0/484 was right
for THAT window (`J ∈ [3.9, 5.5]`), but the conclusion was false: the
window fell **between** two inversions. As a function of J (at `a=0.9`,
given `A`) `Δr = r_t − r_τ` changes sign TWICE:
- at **small J** (`J ~ 1.6–1.9`): rotational inversion (R10, plane
  `(a,J)`, `ρ(r*)²=H_t0/H_τ0`);
- at **large J** (`J ~ 16–24`, for `A ≳ 1.2`): **conformal inversion**
  (this result).
The scan `[3.9, 5.5]` lay in the middle (t deeper) → 0 cells, but
the inversion is certainly accessible on both sides.

**Inversion condition in closed form** (r_min^τ = r_min^t, i.e.
both branches turn at the same r; from the turning conditions `p_r=0, H=0`):
```
cond(r,A) = A² b Q − P̄ (Ê − A² f) = 0 ,   Q = b v̄² + √(v̄² Δ)
J_inv = P̄/Q
```
**Theorem (existence + reachability)**. At the boundaries in A:
- `cond(A=1) < 0` (Kerr, "t deeper", no inversion);
- `cond(A_freeze) = Ê r Δ (Ê−1)/(r−2M) > 0` for r>2M, Ê>1 (EXACT:
  at `v̄²=0` one has Q=0 and the bracket `Ê−A²f → Ê−Ê² = −Ê(Ê−1)`, with the
  factor `1/(r−2M)`).

By the intermediate value theorem there is **`A_inv(r) ∈ (1, A_freeze(r))`**
for every turning radius r>2M: the conformal inversion **always exists
and lies BELOW freezing (reachable)**. Far field `A_inv → √Ê`
(consistent with R10: A>√Ê makes the inversion generic).

**Checks**:
- direct integration of the two branches at `(A_inv, J_inv)`: `r_min^τ =
  r_min^t` to ~1e-4 (e.g. r*=6: A=1.324, J=16.76; r*=9: A=1.236, J=24.27);
- extended `(J, A)` map (J up to 30): 1029/3271 cells in inversion;
  numerical contour `Δr=0` ≡ analytic curve `cond=0` point by point.

Honest magnitude: in the red region τ is deeper only slightly (`Δr ~ +0.1÷0.3`), at
large J (weakly bound orbits, turning at r*~6–9 outside the ergosphere):
a genuine but "gentle" inversion, not a deep plunge.

### R10c. The inversion depends on the PROTOCOL: it disappears with fixed endpoints

> **[SUPERSEDED, September 2026 — see R13]** The "theorem" below does not hold: the pointwise ratio
> n_t/n_τ > 1 does not order the minimum radii of two different functionals. Current status: Lemma (A) (dominant t
> integrand at equal turning radius, proved, symbolic M) plus Lemma (B) on the regions of
> `prop:lemmaB`, at frozen A and for the symmetric single-periastron family; CAP certificates at
> seven static configurations; conjecture elsewhere. The numerical checks below remain as evidence.
(`KerrMetric/colormap_spin_fixed_endpoints.py` → `fig_colormap_spin_fixed_endpoints`;
`KerrMetric/colormap_conformal_fixed_endpoints.py` → `fig_colormap_conformal_fixed_endpoints`;
BVP: `KerrMetric/bvp_kerr_simmetrico.py`, `KerrMetric/bvp_conforme_inversione.py`,
`KerrMetric/bvp_fixed_endpoints.py`)

R10/R10b compare the two branches at **the same launch** (same r₀, same
instant): there the inversion exists, both rotational and conformal. But the
physically clean comparison is with **fixed endpoints**: A=(r₀,−Φ), B=(r₀,+Φ)
identical for t and τ. In that protocol the inversion does **NOT happen**.

**Theorem (conformal, Schwarzschild through the transfer E_eff=Ê/A).** With the
optical indices `n_τ=√(f/(E²−f))`, `n_t=E/√(f(E²−f))` one has
```
n_t/n_τ = E/f > 1   for every r>2M, E>1,   and decreasing in r.
```
The Randers/Beltrami geodesic gives `r_min` as the outer root of
`N(r)·r = J`; at fixed Φ, n_t is larger everywhere ⟹ the t branch turns
higher ⟹ `Δr = r_min^t − r_min^τ > 0` **always**. No root.

**Numerical checks (all Δr>0, no sign change):**
- **spin**, Kerr Hamilton flow, grid (J,a) 18×20, a∈[0.05,0.95]:
  `Δr ∈ [+0.446, +2.064]`; a map with vertical bands, gradient in J
  (endpoint angle), **flat in a** → the spin hardly matters;
- **conformal**, (J,A) through E_eff, grid 44×40: `Δr ∈ [+0.002, +3.711]`,
  no contour Δr=0 (confirming the theorem);
- **direct BVP** (integration of the branches, root shooting on J at Δφ=2Φ):
  symmetric `r_τ=3.617 < r_t=4.735`; asymmetric `r_τ=4.00 < r_t=5.29`.
  With fixed endpoints it is **τ** that goes lower (the opposite of the same-launch case).

**ASYMMETRIC endpoints** (`colormap_spin_asymmetry.py`,
`colormap_conformal_asymmetry.py`; A=(r_A=10,0), B=(r_B=6,Δφ)). The
verdict does not depend on the symmetry: full columns on the (J,a) and
(J,A) grids again give `Δr>0` everywhere —
spin `Δr ∈ [+0.762, +2.636]`, conformal `Δr ∈ [+0.161, +3.244]`. Maps with
vertical bands (gradient in J, flat in a/A), consistent with the theorem
`n_t/n_τ=E/f`, which is **pointwise in r** and hence independent of the endpoints.

**Protocol-dependence picture:**

| protocol | spin a | conformal A |
|---|---|---|
| **same launch** | inverts (R10) | inverts (R10b) |
| **fixed endpoints** | NO (Δr>0) | NO (theorem n_t/n_τ=E/f>1) |

The plunge inversion is therefore an effect of the **comparison
protocol**, not an invariant geometric property of the two
brachistochrones: it requires freedom in the arrival point. With both
endpoints fixed, t always stays above τ.

## R11. Quasi-constants: universal capture at the wall and the polar programme
(`thakurta_kerr_quasi_constants.py`, figure `fig_thakurta_kerr_quasicostante`)

Attempt: adiabatic invariant `I_r = (1/2π)∮p_r dr` for equatorial radial
librations in the shell (bound regime `E_eff = Ê/A < 1`), with the
natural small parameter `ε = A′/A` per period (the expansion
of the non-autonomous Hamiltonian).

**Negative but exact result — universal capture**: the radial
potential `∝ w_A/Δ` is monotone ⇒ the outer turning point is ALWAYS the freezing
wall; there the orbital speed vanishes (`∝√v̄²`) while the
wall comes down at finite speed: at the first contact one measures
`|ṙ_orbit| = 8×10⁻⁴` against `|ṙ_wall| = 3×10⁻²` (**38×**) ⇒ whatever touches the
wall is captured, radial libration does not exist, `I_r` cannot be
defined. (This strengthens R6/R9: capture by expansion has no escape
on the equatorial τ branch.)

**Where the quasi-constant lives** (programme for the next session,
the direct analogue of the Carter quasi-constants of Kerr):
1. `I_θ = ∮p_θ dθ` — the POLAR action off the equator: the motion in θ is a
   genuine libration (double centrifugal turn), robust at the wall;
   it requires the 3D Hamiltonians (ellipsoidal indicatrix, the same
   construction);
2. **conformal quasi-Carter**: `Q_eff = Q_Carter(E → Ê/A(η))` — the Killing
   tensor of Kerr is conformal in Thakurta-Kerr (exact for light);
   for the massive rail the drift `dQ/dη ∝ A′` can be computed at first order
   and tested with the drift-colormap workflow of the Kerr project.

### R11a. Transfer of the expansions K_t, K_τ at sixth order
The O(a⁶) expansions of the Kerr project (`K_tau_expansion.txt`,
`risultati_paper_summary.txt` §1–2: `K_τ = Q_std + a²f₂p_θ² + a⁴f₄p_θ²
+ a⁶f₆p_θ²`, with the numerical verdict: NLO optimal for τ, bare Q_std
optimal for t) transfer to Thakurta-Kerr **quasi-statically
by substitution** `E → E_eff = Ê/A(η)` in the coefficients (exact conformal
rescaling at constant A). Two remarks:

1. **The freezing wall appears as a pole of the expansion**:
   `DE₀(r) = (E²−1)r² + 2Mr = r²·w(r)` exactly ⇒ the denominators
   `D₀·DE₀`, `Δ₀ = D₀²DE₀` etc. vanish where `w_eff = 0`, that is on the
   wall `r_w = 2M/(1−E_eff²)` — which in Kerr (E>1) does not exist (w>0
   everywhere, only the pole r=2M) but in Thakurta-Kerr (E_eff<1) enters from the
   far field and ADVANCES towards the black hole with A(η): the radius of
   convergence of the quasi-constant series is eaten by freezing.
2. The total drift has two small parameters: `O(a²)` (Kerr residual,
   unchanged) + `O(A′/A)` (conformal drift): a hierarchy that can be tested with the
   drift colormaps. The saturation at NLO and the divergence for a ≳ 0.6
   (Kerr verdicts) must be re-mapped as functions of E_eff.

Consistency note: the plunge-inversion criterion of
`risultati_paper_summary.txt` §4, `(f_τ/E)²(H₀τ/H₀t)((J−A_φ)/J)² > 1`,
coincides exactly with our R10 `ρ(r*)² = H_t0/H_τ0` at A=1.

### R11b. Check of the transfer (`thakurta_kerr_K_expansion.py`)
Off-equatorial 3D H built; `K_NLO(E→E_eff)` tested:

- **bracket**: `{K_NLO, H}` stays O(a²) also at E_eff<1 (scaling
  ratio exactly 4.00) — consistent with the non-integrability in section 3 of the
  Kerr summary: the NLO term reduces the drift along the trajectories, not the pointwise
  order of the bracket. The transfer by substitution is formally
  exact (conformal rescaling, constant A).
- **3D drift test** (24 inner trajectories r₀∈[4,7], the same protocol
  in the two regimes): Kerr control (E=1.2, a=0.4): NLO wins 19/24;
  **Thakurta-Kerr (E_eff=0.96, A=1.25): NLO wins only 9/24, ratio
  0.999 — the NLO correction LOSES effectiveness in the conformal regime
  E_eff<1** (a verdict of the t-branch type: bare Q_std is better). Plausible:
  `N₂ = (E²−1)r² + 4Mr − 4M²` changes its sign structure for E<1.
- **pole at the wall verified**: f₂ grows by 100× approaching
  `r_w = 2M/(1−E_eff²) = 25.51` (E_eff=0.96).

Open: carry the ORIGINAL Kerr numtest (the same ensembles/drift
metrics) into the regime E_eff<1 for a 1:1 comparison with their 0.916; understand
whether an NLO correction adapted to the conformal regime exists (re-optimise
the ansatz with the pole at the wall).

### R11c. Structural identity: where it holds and where it breaks
**Constant A — THEOREM** (2 lines, off the equator): `u_TK = u_K/A` ⇒
the constraint Ê on ∂_η ≡ the constraint E_eff = Ê/A on ∂_t; `dτ_TK = A dτ_K` ⇒
`T_τ^{TK} = A·T_τ^K(E_eff)` as FUNCTIONALS ⇒ `H_τ^{TK} = A⁻²H_τ^K(E_eff)`
and (t branch, without prefactor) `n²_{TK} = n²_K(E_eff)` identically.
The perturbative construction commutes with the substitution: the
quasi-constants CANNOT be structurally different.

**Dynamical A(η) — it breaks at the mixed order a²ε (ε = A′/A)**:
`dK/dλ = {K,H} + ∂K/∂η`; `Q₀` does not depend on E ⇒ no counterterm
at O(ε); the first new term:
`∂K_NLO/∂η = ε·a²(2E_eff²cos²θ − E_eff·∂_E f₂·p_θ²) + …`.
The counterterm `h` solves the transport equation
`{h, H₀} = −∂K/∂η|_{a²ε}` along the unperturbed flow: a generic structure
OUTSIDE the Kerr family (flow integrals + pole at the wall):

```
K^{TK} = K^{Kerr}(E_eff(η)) + (A′/A)·a²·h(r,θ,p) + O(a²ε², a⁴ε)
```

`h` = the first object of the quasi-constant programme with no analogue in
Kerr; a probable explanation of the verdict R11b.

### R11d. The computations for non-constant A(η) (`thakurta_kerr_dynamic_K.py`)
**Dynamical 3D PMP Hamiltonian** (new, off the equator, generic Σ):

```
H_τ = p̃_φφ′₀ + R·√( (Δ/Σ)p_r² + p_θ²/Σ + p̃_φ²/Ḡ ) − A²f_Σ/Ê
f_Σ = 1−2Mr/Σ,  b = 2Mar sin²θ/Σ,  Ḡ = G + A²b²/Ê²,
φ′₀ = b·v̄_Σ²/Ḡ,  p̃_φ = J − A²b/Ê,  R² = v̄_Σ²·Δsin²θ/Ḡ
```

with the key identity `f_Σ·G + b² = Δsin²θ` (determinant of the t-φ block, exact in
SymPy). Equatorial reduction ≡ R7 (0.0 at random points).

**Dynamical checks** (`a(η) = −1/(H_cη)`):
- a=0: `L² = p_θ² + J²/sin²θ` conserved to **6e-13** along the non-autonomous
  flow (exact spatial Noether theorem — validation of the 3D machinery);
- hierarchy of the candidates (a=0.4, A: 1→1.25): the running `K(E_eff(η))`
  clearly beats `K(frozen E)` (6.8e-3 against 9.6e-3): **the
  running substitution is the zeroth-order counterterm**; NLO ≈ Q_std
  in the dynamical case (consistent with R11b).

**Source of the counterterm (compact form, verified in SymPy)**:
using `DE₀ − N₂ = −2M·D₀/r` ⇒ `∂_E f₂ = −2MEcos2θ/(r·DE₀²)`:

```
dK/dη|_expl = ε·a²·S₁ ,   ε = A′/A
S₁ = 2E_eff²·[ cos²θ + M·cos2θ·p_θ²/(r·DE₀²) ]
```

`DE₀ = r²w_eff` ⇒ **double pole on the freezing wall**: the
counterterm `h` (solution of `{h, H₀^{Thak-Schw}} = −S₁`) inherits the
pole — a structural confirmation that h lies outside the Kerr family.

### R11e. Transport solved: closed averages, counterterm, physical limits
(`thakurta_kerr_h_counterterm.py`)

**Angular averages in closed form** (orbital plane, `cosθ = sin i·sinψ`,
`k = 1−J²/L²`) — verified by quadrature to **1e-16**:

```
⟨cos²θ⟩ = k/2 ,  ⟨1/sin²θ⟩ = L/|J| ,  ⟨p_θ²⟩ = L(L−|J|) ,
⟨cos2θ·p_θ²⟩ = |J|(L−|J|)
⟹ ⟨S₁⟩ = E_eff²(1−J²/L²) + 2ME_eff²|J|(L−|J|)/(r·DE₀²)
```

**Structural result (extends R11)**: for `E_eff < 1` the τ branch has
no bound orbits EVEN at constant A: `W_eff(r)` decreases
monotonically towards the wall ⇒ **the wall is a global attractor** (endpoint
at finite η with `p_r → ∞`, a coordinate singularity of the reduction).
The dynamical quasi-constants therefore live on ARCS (scattering,
E_eff>1) or require the t branch.

**Counterterm and numerical verdict** (scattering arcs, a=0.15,
ε≈0.002): `K_comp = K(E_eff(η)) − a²∫εS₁dη` removes the explicit drift
*exactly by construction*; on realistic arcs, however, the conformal drift
is **subdominant to the end-to-end Kerr residual**
(6.1e-3 against 6.7e-3 static: the same order). The secular law ⟨S₁⟩ holds
as a torus average: on short arcs (1–3 polar oscillations) the
exact/prediction ratio is ~0.50 (spread 0.36–0.83) — incomplete
mixing, an expected factor.

**Drift map (θ, r)** (`make_drift_map.py` →
`fig_thakurta_kerr_drift_map`): relative pointwise drift
`|{K,H₃D} + εa²S₁|/|K|` on a (θ, log r) grid, s=0.9, A=1.25,
ε=0.025, (p_r, p_θ, J) = (−0.3, 1.0, 1.2). Structure: a red band at the
pole r=2M (D₀ of f₂) next to the ergosphere; a green quasi-equatorial
basin at moderate r; a north–south asymmetry from p_θ>0; green
filaments = loci of cancellation; region beyond the wall masked.
Range: 1e-9 → 2e2, median 7e-3. Variant in the **(A, r)** plane at θ=60°
(`make_drift_map_A.py` → `fig_thakurta_kerr_drift_map_A`): the wall
`r_w(A)` comes down from the right (mask boundary ≡ analytic curve, a free
check), the pole r=2M at the bottom, the green central basin compressed between the
two poles as A grows — the portrait of the "radius of convergence
eaten by freezing" (R11a).

**Summary of the computations for non-constant A(η)**: the dynamical quasi-constant
is `K^{TK}(η) = K^{Kerr}(E_eff(η)) − a²∫ε S₁ dη` with a closed `S₁` (R11d) and
closed averages; the limits of its use are physical, not technical: the
bound regime E_eff<1 does not exist (attracting wall) and on the arcs the
conformal term is a fine correction to the Kerr residual.

## R12. Algebraic structure of the equatorial orbits (elliptic integrals?)
Rationalised radicand of `dφ/dr` (verified in SymPy):

| object | degree of the radicand | genus | closed form |
|---|---|---|---|
| equatorial geodesics | 4 | 1 | elliptic (classical) |
| τ brachistochrone (Schw. AND Kerr) | 5, squarefree | **2** | genus-2 theta (Hackmann–Lämmerzahl) |
| t brachistochrone | 7 | **3** | beyond |

> **[SUPERSEDED — see R12h]** The row "t brachistochrone: degree 7, genus 3" is wrong: the rationalised
> radicand is the sextic R₆ = r·Q₂·[(E²−1)r+2M], genus 2.
| **separatrix J = J_c = a/E** | 4 (doranTau factorisation `Δ−J_c²w = f(r²+a²/E²)`, residual 0) | 1 | **elliptic (Weierstrass ℘)** |
| null limit | 4 | 1 | elliptic (photons) |
| thresholds (double roots) | degenerate | 0 | elementary (consistent with doranT) |

The rail factor `√(wf)` raises the genus by 1 compared with the
geodesics; the rational 𝒦(r) of the t branch by one more. The only
classically integrable brachistochrone is the penetrating separatrix.
Conformal (constant A): `E → Ê/A` ⇒ an elliptic island at `J_c(A) = sA²/Ê`,
crossed dynamically when `A(η) = √(JÊ/s)`. The times T_τ, T_t
share the radicand ⇒ the same classification.

### R12a. φ(r) of the separatrix in closed form (Weierstrass)
(`../KerrMetric/kerr_separatrix_weierstrass.py`, validation
`kerr_separatrix_validation.py`: reduction 5.6e-17, ℘ identity 9e-28,
WW antiderivative 20.53 exact to 30 digits)

Structural surprise: `P₂ = r²+c²` DIVIDES Q₄ ⇒ the residues at the complex
poles vanish: only the poles at the horizons r± remain. With
`c = a/E`, `Q₄ = r((E²−1)r+2M)(r²+c²)`:

```
A = Ma²/(2E²) ,  B = a²(E²−1)/(12E²)
g₂ = a⁴(E²−1)²/(12E⁴) − a²M²/E²
g₃ = a⁴[36E²M²(1−E²) − a²(E²−1)³]/(216E⁶)

z(r) = ℘⁻¹(A/r + B; g₂, g₃)      [= ∫₀^r dr′/√Q₄]

α± = r±((E²−1)r± + 2M)/(r± − r∓)
λ± = α±/√Q₄(r±) ,   v± = z(r±) ,  ℘(v±) = A/r± + B
Λ₀ = (E²−1) − α₊/r₊ − α₋/r₋ + 2λ₊ζ(v₊) + 2λ₋ζ(v₋)

φ(r) = (a/E)[ Λ₀·z(r) + λ₊ ln(σ(z−v₊)/σ(z+v₊))
                       + λ₋ ln(σ(z−v₋)/σ(z+v₋)) ] + const
```

Linear in z + two ratios of Weierstrass sigma functions (Legendre
equivalent: 1 F + 2 Π with characteristics at the horizons). Conformal:
E → Ê/A(η) in all the parameters.

### R12b. Domain of validity and Doran coordinates for Thakurta-Kerr
- **The formula R12a holds ACROSS the ergosphere**: the doranTau
  cancellation `Δ−J_c²w = f(r²+c²)` removes every trace of r_e from the
  radicand (`Q₄(2M) = 4M²(4E²M²+a²) > 0`) and from the poles (which sit at the
  horizons r±, inside the ergosphere). This is exactly what makes the
  separatrix special. The sub-critical orbits (J<J_c), instead, have
  the cusp at r_e, where the optical REDUCTION degenerates (the worldline is
  regular, as in doranTau §5).
- **The breakdown is at the HORIZON**, and it is the known BL artefact: the
  third-kind terms diverge logarithmically as r→r₊ (infinite azimuthal
  winding), as for geodesics in BL coordinates.
- **Doran coordinates for Thakurta-Kerr exist by conformal lifting.**
  `g_TK = A(t_D − F(r))²·g_Kerr-Doran` with the Doran transformation
  of the base Kerr metric (t_D = η + F(r), φ_D = φ + G(r)): the conformal Killing field `∂_η ≡ ∂_{t_D}`
  is the SAME vector field (unchanged by the change of chart), the
  rail constraint is coordinate-free, and the metric is regular
  at the horizon; the only price is that A = A(t_D − F(r)) depends also
  on r in the new chart. The separatrix in Doran coordinates adds the shift
  `∫aβ/Δ dr` with a CUBIC radicand `2Mr(r²+a²)` (verified): a second
  elliptic system (genus 1) — φ_D stays in closed form, a combination
  of two elliptic curves (Q₄ and the Doran cubic).

### R12c. Doran shift in closed form: the complete φ_D
(`../KerrMetric/kerr_doran_shift_weierstrass.py`; V1-V2 exact symbolic,
V3 reduction to 1.1e-16, σ/ζ antiderivative already validated to 30
digits)

`aβ/Δ·dr = 2Mar·dr/(Δ√C₃)` with `C₃ = 2Mr(r²+a²)` (the factor (r²+a²)
cancels: poles only at the horizons). A second elliptic curve in
CLEAN Weierstrass form:

```
r(z) = (2/M)·℘(z; g₂, g₃) ,   g₂ = −M²a² ,  g₃ = 0   (quasi-lemniscatic)
z(r) = ℘⁻¹(Mr/2)
at the horizons: r±²+a² = 2Mr± ⇒ √C₃(r±) = 2Mr±, ℘′(v±) = ±M²r±

φ_shift(r) = [a/(r₊−r₋)]·( B₊(z) − B₋(z) )
B_k(z) = 2ζ(v_k)z + ln[σ(z−v_k)/σ(z+v_k)] ,   ℘(v±) = Mr±/2
```

**φ_D(r) = φ_BL(r) [R12a, curve Q₄] + φ_shift(r) [cubic C₃]**: the
Doran separatrix in closed form on TWO elliptic curves; the
log divergences of the two pieces at r₊ cancel (the analytic content
of Doran regularity at the horizon). Conformal: E → Ê/A only in the
BL part (the shift does not contain E).

### R12d. End-to-end check of the closed form and the PMP limit at the ergosphere
(`../KerrMetric/kerr_separatrix_trajectories.py`,
`fig_separatrix_3trajectories`)

Three trajectories of the separatrix, with the same initial conditions (r₀=9,
a=0.9, E=1.2), down into the ergosphere:
- **T3, the Weierstrass form actually evaluated** (σ,ζ,℘ through Jacobi theta functions on the
  TRUE lattice: rhombic, τ = 1/2 + 0.4639i, ω₁ = 3.1305; residual of the
  ℘ identity 8.4e-15, uniformisation 1.5e-11) against **T2, quadrature**:
  **max |Δφ| = 2.8e-10 across the ergosphere** — the formula verified
  end to end, special functions included;
- **T1, Hamilton ODE**: it follows the others to ~1e-4/1e-5 and grazes
  `r_min = 2.000000` (EXACT marginal turn at r_e) — it cannot
  cross: **new result**: inside the ergoregion `H(p) > 0`
  for every p (the cost `𝔉 = f + bφ′` changes sign on the spurious branch
  of the ellipse): the PMP-η description of the τ branch ENDS at r_e — the
  exact Hamiltonian version of the "1/F singularity intrinsic to the
  reduction" of doranTau. Inside, the closed form follows the orbit of the
  constrained worldline (anchored by the Doran validation of doranTau,
  max|Δφ| = 4.4e-6).

### R12e. Gallery: the 3-method overlay is general, not special
(`../KerrMetric/kerr_separatrix_gallery.py`, `fig_separatrix_gallery`)

The R12d test repeated for FOUR parameter pairs (a,E), each with ITS OWN
lattice (different g₂,g₃,τ). Automatic search of the lattice type
(rhombic Re τ=1/2 / rectangular Re τ=0) with validation through the
uniformisation ℘(z(r))=A/r+B and golden-section minimisation of the
℘ residual.

| a | E | lattice τ | ℘ residual | unif. | \|W−quad\| | graze | \|ODE−q\| |
|---|---|---|---|---|---|---|---|
| 0.9 | 1.20 | 0.5+0.4639i | 6.6e-14 | 1.5e-11 | 2.8e-10 | 2.000000 | 1.5e-4 |
| 0.5 | 1.20 | 0.5+0.4795i | 2.7e-15 | 2.5e-11 | 1.2e-10 | 2.000000 | 2.6e-4 |
| 0.9 | 1.05 | 0.5+0.4901i | 7.0e-14 | 2.1e-11 | 2.5e-10 | 2.000000 | 1.5e-4 |
| 0.7 | 1.50 | 0.5+0.4387i | 3.8e-15 | 6.0e-11 | 2.8e-10 | 2.000000 | 2.2e-4 |

- **Evaluated Weierstrass form ≡ quadrature to ~1e-10 in all cases**: the
  σ/ζ closed form (R12a) is not an artefact of a particular lattice —
  it holds on different lattices. The ODE/quadrature/℘ overlay is NOT
  a lucky case of (0.9,1.2).
- **The lattice is always rhombic** (Re τ=1/2): for E>1 the discriminant is
  <0 in all four cases. The τ separatrix lives in the rhombic regime; the
  rectangular branch exists in the code but is not selected here.
- **The ODE always grazes r=2.000000M = r_e** (equatorial ergosphere): the
  PMP-η termination at the ergosphere (R12d) is confirmed in every case —
  a structural phenomenon, not a parametric one. The residual \|ODE−quad\|~2e-4
  lies entirely in the turning region near r_e; outside it is at the level of the
  quadrature.

### R12f. Retrograde trichotomy: the smallest negative J for capture
(`equatorial_trichotomy.py`, `trichotomy_figure.py`,
`fig_thakurta_trichotomy_Jneg`)

Capture band (the ingoing orbit REACHES the ergosphere `r_e=2M`;
by R12d the τ branch ends there) as J varies, including J<0, for the two
equatorial branches. M=1, a=0.9, Ê=1.2, r0=8.

- **τ branch**: **SYMMETRIC** band `J ∈ [−sA²/Ê, +sA²/Ê]`. Edges:
  A=1 → `±0.7511` (against the closed separatrix `J_c=sA²/Ê=0.75`, diff 1e-3);
  A=1.1 → `±0.9088` (against 0.9075). The symmetry is not accidental: the
  doranTau factorisation `Δ − J_c²w = f(r²+c²)` depends on **`J²`**,
  so prograde and retrograde orbits have the SAME threshold despite
  frame dragging. ⇒ **`J_neg^τ = −sA²/Ê`**.
- **t (η) branch**: a **wide and ASYMMETRIC** band. Capture down to
  `J≈−8` (strongly retrograde); prograde capture beyond `+3`. Negative
  edges: A=1 → `−8.05`, A=1.1 → `−8.37`. Below `J≈−15` the orbit cannot
  start (no ingoing root).

Physics: the t branch is enormously more *capture-prone* (consistent with "t
goes deeper", cf. Vaidya R13/R15). A retrograde orbit with `J=−5`
is **captured by t but scattered by τ**. Frame dragging breaks the
symmetry only for the t branch; the `J²` structure of doranTau protects the
symmetry of the τ branch. `A>1` slightly widens both thresholds.

**Analytic form of the t thresholds (new):** the turning function of the t branch
(from H_η=0 at p_r=0) is `N(r)=(P̄−Jbv̄²)²−J²Δv̄²`, `b=2Ma/r`. At r_e:
`N(r_e)=P̄_e(P̄_e−2aJ)`, `P̄_e=4M²+2a²+a²/E²` ⇒ **a single** crossing
(prograde) at **`J_+^t=P̄_e/(2a)`** (=3.435 for a=0.9,E=1.2; ≡ "beyond +3").
No retrograde crossing at r_e. The retrograde edge is the **marginal
capture** (double root `N=N'=0` at r*>r_e, the analogue of the retrograde
photon sphere): **`J_c^-=−P̄(r*)/[v̄(√Δ−bv̄)]|_{r*}`** = −8.054 (r*=3.514),
an exact match with the numerics. t band = `[J_c^-, J_+^t]`; conformal:
E→Ê/A, a→s, J→J/A.

### R12g. The reflection at the ergosphere is a CUSP (proof)
(`ergosphere_cusp.py`, `fig_thakurta_ergosphere_cusp`)

Analytic proof that an equatorial τ brachistochrone with
`|J| ≠ J_c` that reaches the ergosphere is reflected there with a **cusp**
(not a smooth periastron). From
```
dφ/dr = J √(w f) r / (Δ √(Δ − J² w)) ,   w = E² − f
```
near `r_e=2M` (where `f→0`, `f ≈ (r−r_e)/2M`, `w→E²`, `Δ→a²`):
```
dφ/dr ~ K √(r−r_e) ,  K = J E r_e/(a²√(a²−J²E²)√(2M))
⟹ φ − φ_e ~ (2/3) K (r−r_e)^{3/2}    (power 3/2 = CUSP)
```
The cusp comes from `√f→0` at the ergosphere while the turning radicand
`√(a²−J²E²)` stays finite for `|J|<J_c`. **Radial** tangent (`dφ/dr→0`
from both sides) ⟹ a cusp, not a corner.

**Trichotomy at the ergosphere** (sign of `a²−J²E²`):
- `|J|>J_c` (`a²−J²E²<0`): turn BEFORE `r_e`, smooth periastron,
  classical exponent 1/2 (`dφ/dr→∞`);
- `|J|=J_c` (`a²−J²E²=0`, with `J_c=a/E`): `√f` and `√(Δ−J²w)` vanish
  together ⟹ FINITE `dφ/dr`, smooth crossing (separatrix);
- `|J|<J_c` (`a²−J²E²>0`): `dφ/dr~√(r−r_e)` ⟹ CUSP.

Numerical checks: symbolic coefficient `K` ≡ closed formula;
exponent of `dφ/dr` at `r_e`: J=0.5→0.4999, J=0.7→0.4997 (cusp),
J=J_c→−0.0001 (smooth), J=0.9>J_c→turn at r=2.141>r_e; and
`Δφ/(r−r_e)^{1.5} → 2K/3` (power 3/2). Thakurta-Kerr through `E=Ê/A`,
`J→J/A`, `J_c=sA²/Ê` (r_e conformally invariant). It puts on a rigorous basis
the "cusp reflection" of the trichotomy (§ doranTau, R12f).

> **[Clarification, September 2026]** The two marginal cases are not equivalent (`prop:classification`):
> J = +J_c crosses r_e at finite speed (a corner, not a cusp), J = −J_c reaches it only
> asymptotically. "Smooth crossing" above applies to the +J_c branch.

**Crossing in Doran coordinates** (figure, green branch J_c). The
`√(wf)` becomes imaginary inside the ergosphere (`f<0`), but on the separatrix
`√f` **cancels**: `dφ_BL/dr = J_c√w·r/(Δ√(r²+c²))`, `c=a/E` — regular
across `r_e`. To regularise the HORIZON one passes to Doran coordinates:
`φ_D = φ_BL − φ_shift`, `dφ_shift/dr = a√(2Mr/(r²+a²))/Δ` (R12c); at `r_+`
the two `1/Δ` divergences cancel (verified: BL−shift = +0.036
finite, BL+shift diverges). So the `J_c` orbit **crosses `r_e` and
goes down to the horizon** with finite winding — the green curve appears "on the other
side" of the ergosphere (impossible in BL coordinates, which stop there).

### R12h. Separatrix of the t branch: genus 2 in general, genus 1 at the saddle-node
(`../KerrMetric/kerr_tbranch_separatrix_weierstrass.py`)

**CORRECTION + extension of R12.** The earlier statement "the t branch has a
rational 𝒦 ⟹ no closed form" was **imprecise about the cause**: a rational 𝒦
gives only integrals of the **third kind**, which close anyway in ℘,ζ,σ.
What matters is the **genus of the radicand**.

- **Generic t branch**: `dφ/dr = 𝒦/√R₆`, `R₆ = r·Q₂(r)·[(E²−1)r+2M]`
  **sextic ⟹ genus 2** (hyperelliptic). No doranTau cancellation
  (the gravitomagnetic wind does not cancel). Verified symbolically
  (degree 6) and numerically (`dφ/dr=𝒦/√R₆`).
- **Dynamical separatrix of the t branch** = the **retrograde saddle-node**
  `J_c^-` (R12f): there `R₆=R₆'=0` ⟹ a **double root** `(r−r*)²` at
  `r*` (marginal). Hence `R₆=(r−r*)²Q₄`, `Q₄` quartic ⟹ **genus 1**.
  A **different** mechanism from the τ separatrix (doranTau `(r²+c²)|Q₄`),
  the same elliptic outcome.

**Explicit closed form** (the same structure as eq. A3 of R12a):
```
z(r)   = ℘⁻¹(A/r + B; g₂, g₃)     [invariants of Q₄]
φ(r)   = Λ₀ z + Σ_k λ_k ln(σ(z−v_k)/σ(z+v_k))
poles  c_k ∈ {r₊, r₋, r*}   (THREE real ones; NO pole at r=0, NO pair ±ic
             — cleaner than the τ separatrix, which has 4 poles)
v_k=℘⁻¹(A/c_k+B),  λ_k=−A α_k/(c_k² ℘'(v_k)),
Λ₀=Λ_poly − Σα_k/c_k + Σ 2ζ(v_k)λ_k
```
Sample M=1, a=0.9, E=1.2 (r*=3.5139, J_c^-=−8.0535): g₂=−11.650,
g₃=−2.578, A=12.985, B=2.639; residues α: r₊→+13.23, r₋→−6.59, r*→+21.13;
Λ_poly=+3.544. **Check**: φ from the Hamiltonian flow against the decomposition
= +1.557010139, **diff 2×10⁻¹⁴** (rational-fit RMS 3×10⁻¹³).

**Final picture**: both separatrices (τ, t) have **genus 1** but through
different mechanisms (cancellation against a double root); only the
**non-separatrix** t orbits are genuinely genus 2 (Kleinian). The prograde
crossing `J_+^t` has a **simple** root at r_e ⟹ it stays genus 2.

### R13. Results of the audit of the proofs (September 2026)

- **Non-inversion with fixed endpoints (Lemma B):** it holds where Φ_τ decreases on [r_min^t, r₀) (frozen A,
  symmetric single-periastron family). Proved regions: (i) a=0, V ≥ c*V₀, c* = sech²U*,
  U* tanh U* = 1 (0.30518, optimal); (ii) a=0, r ≥ 9M/4, V ≥ V₀/4; (iii) r ≥ 5M/2, V ≥ c*V₀; (iv)
  r ≥ 3M, V ≥ V₀/4; (v) E² ≥ 3/2: a single peak. W′ < 0 for r ≥ 5M/2. False in general: the quarter without
  a radial margin (b=1, r₀=2.01) and unimodality (b=1/100, r₀=10⁴). R_* withdrawn. CAP RUN 6: 7
  configurations with matching at V₀/2. `verify_lemma_B.py`, `no_inversion_schwarzschild_CAP_grid.py`.
- **Proof A in symbolic M:** J* = [E²r(a²+r²)+2Ma²]/(aD_E), J** = [E²r²(r−M)+E²a²r+Ma²]/(aD_E),
  N_G = (r−r_min)·4MrΔ/[(r_min−2M)DE_min]·W; it also holds at a=0. `no_inversion_reduction.py`.
- **Peak of Φ_τ:** a fold of the map, not a circular τ orbit (∂_r[rΔ − J²D_E] = D_E V′ > 0 on the shell).
- **Costates:** J_eff = p_φ/A for the t and τ branches (factorisation of the Hamiltonian); J_η = J_t/A;
  (J_τ/A)/J_η = f/E at a=0.
- **Optical metric with varying A:** ∂_η a = 2Ê²AA′h₀/(Ê² − A²f₀)² ≠ 0; the foliation remains, the cost
  does not reduce to a single autonomous metric.
- **Not Einstein:** Ric_{ηr} = 2(A′/A)Γ̄^η_{ηr}, Γ̄^η_{ηr} = M(r²+a²)/(r²Δ) on the equator; hence W
  is not a selector of an (unweighted) Ricci soliton for varying A.
- **Bakry–Émery identity (30 September 2026):** generalised m-Bakry–Émery tensor, defined
  algebraically, at m = 2−n = −2 (negative parameter; no quasi-Einstein result imported),
  φ = 2 ln A. It is not the conformal rescaling of the metric of Lecian's sources. Ric_φ^m = μ g,
  μ = −ḡ^{ηη}(ln A)″/A² (Ricci-flat seed, harmonic BL time where the chart is regular). Then
  ½L_W g + Ric_φ^m = λ g, λ = A′/A + μ: an almost-soliton-type equation in this convention. For M > 0,
  λ constant on spacetime ⇔ A = A₀e^{cη} (λ = c), because on the equator
  ḡ^{ηη} = −1 − 2M(r²+a²)/(rΔ) depends on r; the flat seed is excluded. Geodesic drift
  λ + Ric_φ^m(u,u) = A′/A; on the controlled rail g(∇_u u, W) = A′/A. It characterises the selector,
  it does not derive it. `verify_weighted_soliton.wls` (complete Ricci tensor of Thakurta–Kerr and of the case a = 0,
  an unweighted control that must fail).
- **Tension of the projection (Lorentzian):** τ^r = (M/r²)(2Ê²−3f)/(Ê²−f) (a=0, A=1), it changes sign for
  Ê² < 3/2; the Riemannian assembly μ + ∇ ln Λ is not the tension.
- **Optical curvature spectrum (frozen Schwarzschild):** K_t, K_t, K_tan; R_opt = 2(2K_t + K_tan);
  out of the plane cos²χ K_t + sin²χ K_tan. Null metric: K = −M(2r−3M)/r⁴. `verify_myers_optical.wls`.
- **Taş surfaces:** a regular family of extremals at frozen A generates a timelike ruled
  surface (not immersed at the common launch); with L = F_T²/2 the magnetic separation is the normal
  component of the linearised equation of the generators. Kerr and the non-stationary case are open in Taş.
- **W_sup:** an algebraic question (√(ΔD_E/r) remains); ingoing branch with σ_r = −1 in the radial action.

## Open questions
1. Dynamics on the freezing surface (after the contact the
   Hamiltonian formalism degenerates: a description constrained to the
   boundary is needed — an analogue of the J<J_c cusp of Kerr?).
2. Escape map: (J, Ê, η₀) → escape against capture by expansion; threshold
   Ê_c(r₀, η₀).
3. Sultana–Dyer (a ∝ η², exact source) as a physical variant.
4. Quantitative comparison with asymptotic McVittie (Λ>0 → SdS).
