# Brachistochrones in Vaidya — results

**SymPy scripts** (all checks pass):
- `vaidya_metric_sympy.py` — study of the metric, its source, horizons
- `vaidya_brachistochrone_sympy.py` — equations of motion of the constrained worldline (parameter r)
- `vaidya_brachistochrone_vparam.py` — equations of motion in the parameter v, bounce, timing
- `vaidya_penetration_map.py` — penetration map J_c(v₀)
- `vaidya_hamiltonians_sympy.py` — PMP/Zermelo Hamiltonians of the two branches
- `vaidya_optical_metric_sympy.py` — optical metric through the homothety (m = μv)
- `make_figures_vaidya.py` — validation figures (`Vaidyafigures/`)

Context: `../FLRWmetric/nonStationaryBrachi.md` (plan), `FLRWresults.md`
(FLRW case), `../KerrMetric/doranTau.md`, `doranT.md` (stationary case).

---

## 1. The metric and what it describes

**Ingoing Vaidya** (advanced Eddington–Finkelstein coordinates, c=G=1):

```
ds² = −(1 − 2m(v)/r) dv² + 2 dv dr + r² dΩ²
```

Verified: `R = 0`, the only component is `G_vv = 2m′(v)/r²`, that is

```
T_μν = [m′(v)/(4π r²)] l_μ l_ν ,   l_μ = −(dv)_μ ,  l·l = 0,  ∇_l l = 0
```

**ingoing radial null dust** (a spherical rain of radiation).

| variant | energy | describes |
|---|---|---|
| ingoing, `m′>0` | NEC/WEC satisfied (`ρ ∝ m′`) | **black hole accreting** radiation |
| ingoing, `m′<0` | violates the NEC | **toy model of Hawking evaporation** |
| outgoing (`u`, `m′(u)<0`) | positive outgoing radiation | **radiating star** (Vaidya 1951) |

Key properties (verified):
- `m′=0` ⇒ vacuum ⇒ Schwarzschild in advanced EF coordinates; coordinates regular at the horizon.
- `(£_∂v g)_vv = 2m′/r ≠ 0`: **no timelike Killing field** — the intended non-stationary regime.
- **Apparent horizon**: outgoing expansion `θ₊ = (r−2m(v))/r²` ⇒ `r_AH = 2m(v)`,
  with velocity `ṙ_AH = 2m′(v)`.
- **Event horizon ≠ apparent horizon** (numerical, `m = 1+0.02v` on `v∈[0,50]`, then constant):
  during accretion `r_EH > r_AH` (gap ~0.2–0.27), then they converge to the
  final Schwarzschild. The EH is teleological: it swells *before* the mass falls in.

---

## 2. Equations of motion of the brachistochrone (non-autonomous constrained worldline)

**Rail constraint**: `−u_v = E` held fixed. `∂_v` is not Killing: the rail imposes
the constraint with a force whose component `f_v ≠ 0` (the analogue of `f_η = aH` in FLRW).

**Reduction** (parameter `λ = r`, on a stretch with increasing r, `′ = d/dr`).
From the constraint `Λ = dτ/dr = (f v′ − 1)/E`; eliminating Λ:

```
φ′² = W(v, v′, r) = [ f v′² − 2v′ − (f v′ − 1)²/E² ] / r²

τ branch:  F_τ = (f v′ − 1)/E − J √W        (T_τ = ∫ Λ dr)
t branch:  F_t = v′ − J √W                   (T_v = ∫ v′ dr = Δv)
```

with `J` = Fermat momentum (multiplier for Δφ). Equations of motion: Euler–Lagrange in `v(r)`,
second order, non-autonomous.

### R1. What is conserved and what is not
- **`J` conserved**: spherical symmetry survives (φ cyclic). What is lost is `E`
  (which is indeed imposed by the rail), not `J` — the mirror image of
  what one might expect: the Killing field that is lost is the temporal one.
- **Forcing ∝ m′ (exact)**: `F` depends on `v` only through `m(v)` ⇒
  `∂F/∂v = (∂F/∂m)·m′(v)` (chain rule; `∂F/∂m` in closed form in the scripts).

### R2. p_v and the teleology of the rail
`p_v = ∂F/∂v′`. Spatial arrival with a free instant `v` ⇒ transversality
`p_v(r_arrival) = 0`.

- Stationary (`m′=0`): `p_v` conserved ⇒ `p_v ≡ 0` on the orbit.
- Non-stationary: `dp_v/dr ∝ m′` ⇒ `p_v(r) ≠ 0` along the path even with
  `p_v(r₁)=0`: the brachistochrone **anticipates the future accretion**
  (the same teleology as the event horizon, here in the dynamics).

### R3. Validation in the Schwarzschild limit
With constant `m` and `p_v = 0`, `dφ/dr` reproduces the Kerr closed forms at a=0
(numerical match to 12 digits, both branches):

```
dφ/dr = 𝒦 r √(w f) / ( Δ √(Δ − 𝒦² w) ) ,   Δ = r²f ,  w = E² − f
𝒦_τ = J           (τ branch)
𝒦_t = f J / E     (t branch)
```

### R4. Numerical demo (linear Vaidya)
Ascent `r: 3.5 → 10`, arrival at `v = 40`, `E=1.2`, `J=1.3`,
`m(v) = 1 + 0.01v` against Schwarzschild `m=1` (backward integration
from the arrival with `p_v(r₁)=0`; the periastron of the orbit is at `r≈3.01`, where the
parametrisation in `r` degenerates — the start lies above it):

| | Schwarzschild | Vaidya μ=0.01 |
|---|---|---|
| v at start | 16.103 | 13.044 |
| Δφ | 0.4014 | 0.5293 |
| T_v = Δv | 23.897 | 26.956 |
| T_τ | 7.757 | 7.480 |
| p_v at start | 0 (exact) | **0.0188** |

Validation figures in `Vaidyafigures/` (`make_figures_vaidya.py`):
EH/AH horizons and self-similar rays, residuals against the Kerr closed form at a=0,
orbits and `p_v(r)`, direct variational minimum of `T_τ` at fixed `Δφ`.

---

## 3. Optical metric: is it possible?

**Generic Vaidya: NO.** A timelike Killing field is needed (Perlick 1991); there is none, and
there is not even a conformal Killing field: `(£_ξ g − 2g)_vv = 2(v m′ − m)/r` with
`ξ = v∂_v + r∂_r`.

**Linear Vaidya `m = μv`: YES, through the homothety.** The residual vanishes ⇔
`m = μv`: `ξ` is a **homothetic Killing field** (`£_ξ g = 2g`). Self-similar
coordinates `s = ln v`, `x = r/v`:

```
g = e^{2s} · ĝ(x) ,    ĝ = (2x − f) ds² + 2 ds dx + x² dφ² ,   f = 1 − 2μ/x
```

`∂_s` is a **genuine** Killing field of `ĝ` — the same structure as FLRW (`g = a² ĝ`),
with `e^s` in place of `a`.

### R5. Randers optical metric (null brachistochrone in v)
For null rays on `ĝ` (`A = −ĝ_ss = f − 2x`):

```
ds_arrival = [ dx + √( dx² + A x² dφ² ) ] / A
⇒  F = α + β ,   α² = dx²/A² + x² dφ²/A ,   β = dx/A
```

(checked against the null condition of `ĝ`). Minimising
`s_arrival = ln v_arrival` = minimising the arrival time `v`. Since
`v = e^s` is monotone, the **massive t branch** with the homothetic rail
`Ê = −u·ξ` also reduces to a **static** problem on `ĝ`. The τ branch does not:
`∫dτ = ∫e^s dτ̂` is not invariant.

### R6. Light surfaces = self-similar horizons, critical accretion
`A = 0 ⇔ 2x² − x + 2μ = 0`:

```
x± = (1 ± √(1 − 16μ)) / 4 ,     real ⇔ μ ≤ μ_c = 1/16
```

- `x₋` is the self-similar **event horizon** (`r_EH = x₋·v`;
  the null hypersurface `x=const` ⇔ the same quadratic; numerical bisection:
  8-digit match for μ=0.02: `x₋ = 0.04384472`).
- `x_AH = 2μ < x₋`: the AH lies inside the EH, as it must under accretion.
- `x₊` = outer escape surface (attractor of the outgoing rays).
- `μ_c = 1/16` is the critical value known from the self-similar literature
  (Papapetrou / naked singularities): for `μ > 1/16` there are no surfaces `x±`.

The `1/F` singularity of the optical reduction falls on the **homothetic
horizons** `x±` — the analogue of the light-surface/ergosphere pair of
Kerr, but here the inner surface IS the event horizon.

---

## 3b. Trajectories in the parameter v and the time window
(`vaidya_brachistochrone_vparam.py`, figures `fig_vaidya_bounce`,
`fig_vaidya_timing`)

The parametrisation in `r` degenerates at the periastron (`dv/dr → ∞`). With `v`
as the parameter (`′ = d/dv`):

```
Λ = (f − r′)/E ,   φ′² = W_v = [ f − 2r′ − (f−r′)²/E² ] / r²
F_τ = (f − r′)/E − J√W_v
```

The same functional (exact change of variable: `(fv′−1)dr = (f−r′)dv`),
the same extremal curves, but the Euler–Lagrange equation is regular at `r′ = 0`: one integrates
**through the bounce**. Transversality is taken from the r-parametrised form
(`p_v(r₁)=0` ⇒ `r′(v₁) = 1/v′(r₁)`).

**Validations (E=1.2, J=1.3, r₁=10):**
- cross-check with the r-parametrisation at the passage r=3.5: `v`, `Δφ`, `T_τ` agree to
  ~1e-7 (Schwarzschild and Vaidya μ=0.01);
- static: `r_min = 2.72713516` = exact root of `r²f − J²(E²−f) = 0`
  to **9e-12**;
- Vaidya μ=0.01: periastron `r_min = 3.010492` at `v = 7.29` — exactly
  where the parametrisation in r broke down (3.0104915).

**R7. Time window (scan of the arrival instant v₁ ∈ [25,70]):**

| μ | r_min(v₁=25) | r_min(v₁=70) | r_min/2m(v_peri) |
|---|---|---|---|
| −0.01 (evaporating) | 2.68 | **1.89** | 1.27 → **1.63** (moves away) |
| 0 (static) | 2.727 (flat) | 2.727 | 1.36 constant |
| +0.01 (accreting) | 2.76 | **3.51** | 1.48 → **1.30** (moves closer) |

The periastron depends on *when* one arrives — a genuinely non-stationary
phenomenon, absent in Kerr. **Counter-intuitive**: under evaporation the
absolute periastron goes down but the horizon recedes faster ⇒ the
*relative* penetration gets worse; under accretion the periastron rises but
the horizon grows faster ⇒ late arrivals graze the relative
horizon. The orbit–horizon "race" has opposite winners in absolute
and in relative terms.

## 3c. Penetration map (J, v₀)
(`vaidya_penetration_map.py`, figure `fig_vaidya_penetration_map`)

Set-up: ingoing launch from `r₀=10` at the instant `v₀`, locally optimal aim
(`p_v=0` at frozen mass, solved in closed form; time reflection
`u_in = 2/f − u_out` — the radial roots of W add up to `2/f`). Outcomes:
reflects (`r′=0` outside) or captured (crosses `r=2m(v)`). Threshold
`J_c = sup{J : capture}` by bisection. DOP853 integration with rtol 1e-12
(RK45 with 1e-10 is NOT enough for grazing turns: r_min wrong by 4e-4).

**Static validation**: r_min reproduces `2m + J²E²/(2m)` at all the J
tested (0.05→2.001797 against 2.001800; 0.005→2.000019 against 2.000018);
`J_c = 0` — penetration of measure zero, as in the Kerr trichotomy at a=0.

**R8. Main result — accretion opens a finite window:**

| case | J_c(v₀=0) | J_c(v₀=60) |
|---|---|---|
| static | 0 (measure zero) | 0 |
| evaporating μ=−0.01 | 0 | 0 |
| **accreting μ=+0.01** | **0.919** | **1.286** (grows with m(v₀)) |

- Under accretion the τ branch acquires a **finite capture interval**
  `J ∈ (0, J_c]` — the structure that in Kerr belonged only to the t branch
  (dichotomy with a threshold) emerges here from the **dynamics**: grazing orbits
  linger near the periastron and the growing horizon closes the gap
  (a chase; see the figure: two almost identical orbits on either side of J_c
  separate only in the grazing phase).
- Under evaporation the horizon recedes faster than the orbit
  approaches: `J_c = 0` as in the static case (accretion/evaporation
  asymmetry, consistent with R7).
- Empirical scaling: `J_c ≈ (8÷10)·√μ` for μ ∈ [0.0025, 0.02] (v₀=0):
  consistent with the estimate `δ_p = J²E²/2m ~ growth of the horizon during
  the grazing phase`.
- Technical caveat: for `J ≲ 0.01` the reduction degenerates
  (`∂²F/∂r′² ∝ J`): the Euler–Lagrange equation is stiff and the outcomes are not reliable — the probes start
  from J=0.01.

## 3d. Hamiltonians of the t and τ branches (`vaidya_hamiltonians_sympy.py`)

The constraint `−u_v = E` makes the indicatrix (admissible velocities at fixed (v,r))
an **ellipse** — exact identity, `w = E²−f`:

```
r′ = (f − E²) + E√w cosθ ,     φ′ = (√w/r) sinθ
```

Zermelo navigation: **ingoing radial wind** `f − E²` (always < 0 for
E>1) + an ellipse with semi-axes `E√w` (radial) and `√w/r` (angular). The Pontryagin
maximum over θ gives the Hamiltonians in **closed form** (p_φ = J):

```
H_v  = p_r(f − E²) − 1 + √w · √( E²p_r² + J²/r² )
H_τ  = p_r(f − E²) − E + √w · √( (Ep_r + 1)² + J²/r² )
```

(τ branch = v branch with `p_r → p_r + 1/E` and cost E: the cost dτ/dv depends
on the control). Properties:
- `p_φ = J` conserved; `dp_r/dv = −∂H/∂r ≠ 0` (radial inhomogeneity);
- `dH/dv = ∂H/∂v ∝ m′(v)` — the non-autonomous forcing in Hamiltonian form;
- transversality: `H = 0` at arrival; static ⇒ `H ≡ 0` (H = −p_v).

**Validations**: the Hamilton flow of H_τ reproduces the v-parametrised Euler–Lagrange bounce
exactly (r_min, Δφ, T_τ, v_start to 6 digits, μ=0 and 0.01); static
`max|H| = 4e-12`; dynamical `H(start) = ∫∂H/∂v dv` to 6 digits (−0.062126).

**R9 — self-tuning of the v branch, generalised dynamically**:
with the same (E=1.2, J=1.3, arrival v=40) the v branch penetrates down to
`r_min = 2.000000` (static: it grazes the horizon EXACTLY — the
doranT mechanism `𝒦(r_e) = fJ/E → 0` survives at a=0) and under
accretion down to `r_min = 1.04` = grazing `2m(v)` at the epoch of the
periastron (reached at early v, when the black hole was small). The v branch
"tunes itself" to the horizon **of its epoch**.

## 3e. Non-autonomous optical forms for generic m(v)
(`vaidya_optical_nonautonomous_sympy.py`)

Inverting the indicatrix (inverse Legendre transform of the Hamiltonians of §3d), for
**any** m(v) — with `f = f(v,r)`, `w = E²−f`:

```
dτ = √( dr² + f r²dφ² ) / √w                      (Riemann–Jacobi)
dv = dr/f + (E/(f√w)) √( dr² + f r²dφ² )          (Randers: drift + norm)
```

Checked against the constraint to ~1e-14; with `m′=0` they are **exactly the
two Perlick (1991) metrics** in EF coordinates (τ: identical to `L_τ = √((B²+Af)/w)`
of doranTau with B=1; v: the Euler–Lagrange equation reproduces doranT at a=0 to 10 digits): Perlick
generalised *verbatim* with `f → f(v,r)`.

**The point**: the forms always exist, but without a conformal Killing field nothing is rigid —
drift and norm both breathe. "Optical" only in a **non-autonomous**
sense (Lagrangian dual of the Hamiltonians), with no underlying static
geometry (except for `m = μv`, the homothety, R5-R6). Hierarchy:

| case | optical structure |
|---|---|
| FLRW | separable index → clock, trivial geometry |
| Thakurta-Schw | static geometry × index n(η,r) |
| Thakurta-Kerr | static geometry + **rigid** wind × index |
| Vaidya μv | static geometry on the self-similar quotient |
| Vaidya m(v) | everything breathes: weak sense only |

## 3f. The rail W is the Kodama vector
(`vaidya_kodama_sympy.py`)

`K^a = ε^{ab}∇_b r` (on the 2-plane normal to the spheres): in advanced EF coordinates
**K = ∂_v exactly** — our W. Verified:
- `K·K = −f` (timelike outside the apparent horizon, like a Killing field);
- `∇·K = 0` exactly for every m(v);
- **Kodama miracle**: `J^μ = G^{μν}K_ν = (0, 2m′/r², 0, 0)` with
  `∇·J = 0` for any m(v) — a current conserved WITHOUT a Killing field;
- the associated charge = Misner–Sharp mass `= m(v)`.

⇒ The constraint `−u_v = E` is `−u·K = E`: the **Kodama energy**. In every
dynamical spherically symmetric spacetime W has a canonical
definition (it is not a convention of foliation); Kodama → Killing in the
static limit, so the hierarchy of legitimacy of §4b of the kick-off note
is refined: Killing (theorem) → conformal Killing field (conformal selection) → **Kodama
(quasi-local selection, dynamical spherical case)** → declared convention
(no structure).

### R12. Kodama energy conserved along the brachistochrone
(`kodama_conservation.py`, `fig_kodama_conservation`)

Numerical confirmation + exact identity of §3f:
- **[A] symbolic**: `K = ∂_v`, `∇·K = 0`, `K_a = (−f,1,0,0)`, Misner–Sharp
  charge `= m(v)` — all exact for any m(v).
- **[B] rail identity**: from the constraint `W` alone it follows that
  `−g(T,T) = (f u^v−1)²/E²` ⟹ **`−u·K = E` identically**. The Kodama
  energy IS the rail invariant; it is not imposed by hand.
- **[C] numerical check** (accreting Vaidya `m=1+0.02v`, `E=1.2`):
  `|−u·K| = 1.200000000000` constant, `|−u·K − E| < 6.7×10⁻¹⁶`.
- **[C'] non-triviality**: for a geodesic (rail switched off) `−u·K`
  would drift at the rate `−m′(u^v)²/r ∝ m′` (here 0.02–0.14, cumulative 0.33):
  it is the missing symmetry that the rail compensates. It is NOT Noether
  conservation — it is **control** that pays a cost `∝ m′`.

**Discriminating test** (why Kodama and not something else): in the static limit
only the NON-normalised Kodama vector `K=∂_v` has zero cost. With `K̂=K/√f`
one would have `−u·K̂ = E_t/√f`, which drifts even in static Schwarzschild
(where a true Killing field exists) — a sign that it is the wrong choice. Hence
`W=∂_v` is fixed by the criterion **static ⟺ zero cost**.

### R13. r_min(t) − r_min(τ) has a fixed sign (non-rotating)
At the turning points (`p_r=0`, shell `H=0`), at `a=0`:
- τ: `(E²−f)J² = f r²`
- t: `(E²−f)J² = E² r²/f`

Exact local ratio `J²_t(r)/J²_τ(r) = E²/f(r)² > 1` (exterior,
`f<1≤E`). With increasing `J²_branch(r)` ⟹ **`r_min^t < r_min^τ` in the
STATIC limit** (`m′=0`): the t branch goes deeper, with a locked sign.
The plunge inversion of Kerr/Thakurta-Kerr is the effect of the frame-dragging
term `2Mas/r`, which deforms `E²/f²` until it crosses 1.
Verified on a grid in `(E,J)` (Cardano/Ferrari for the exterior roots of the
cubic/quartic). The sign "t deeper" is ROBUST also in the dynamical case:
physical evaporation narrows the gap but does NOT invert it (R16). To
invert it ROTATION is needed (the term `2Mas/r`) — there is no physical
evaporative inversion.

### R14. Direct validation of the minimum principle (three tests)
Three independent variational validations, all by perturbation with
fixed endpoints, with an exact minimum on the solution:
- **FLRW** (`fig_flrw_variational`, R2): the comoving straight line minimises
  BOTH the t and τ functionals (`argmin` at ε=0 for both).
- **Non-autonomous Vaidya** (`fig_vaidya_variational`): the non-autonomous Euler–Lagrange
  solution minimises `T_τ` at fixed `Δφ` (`T_min` against `T_EL` to 1e-8).
- **static t against τ** (`verify_brachistochrone_minimum.py`,
  `fig_brachistochrone_minimum_check`): each branch, built INDEPENDENTLY
  from its own Beltrami ODE, is the minimum of ITS OWN time (`argmin` at ε=0)
  and NOT of the other one (the other functional is monotone in the perturbation).
  A variational confirmation of the t/τ dichotomy (the same one that gives R13).

⇒ The curves found are indeed brachistochrones (true minima,
positive second order), and t and τ are genuinely distinct curves.

### R15. Law of the minimum plunge radius, t against τ, in Vaidya
(`plunge_vaidya_t_tau.py`)

t branch = minimises the advanced time `v` (with fixed radial endpoints this
is equivalent to minimising t: `v = t + r_*`, a constant difference); it reduces
to the Schwarzschild t branch in the static limit. Local law of the turn
(periastron event `p_r=0`), in the form of R13:
```
τ:  (E²−f_*) J² = f_* r_min²
t:  (E²−f_*) J² = E² r_min²/f_* ,   f_* = 1 − 2m(v_peri)/r_min
```
with `J` the optical/Fermat momentum and `E` the Kodama energy (R12).

- **Consistency at m′=0 (Schwarzschild)**: at μ=0 the integration of the two branches
  REPRODUCES the static roots R13 (residual ~1e-9): `r_τ=7.7573`,
  `r_t=4.7011` (E=1.2, J=8, m=1). Continuous limit, no
  discontinuity. Note: for small J the t branch does NOT turn (capture,
  consistent with "t deeper") — a large J is needed for both to turn.
- **Robust sign**: `r_min^t < r_min^τ` for every μ∈[−0.02,+0.02].
  Accretion WIDENS the gap (τ moves out, t goes deeper: gap 1.62→3.06→5.07
  from μ=−0.02 to +0.02), evaporation narrows it — but **does not
  invert it** in the physical regime (R16: the linear inversion at μ≈−0.058 is
  an artefact of `m<0`; with `m(v)>0` the sign stays locked).
- **Dynamical correction**: the law is NOT the static one evaluated at
  `m(v_peri)` — the residual of that form grows ~linearly in m′
  (`res/|μ|` ≈ 475–495 for τ, ~constant): it is the non-autonomous correction
  from the teleological momentum `p_v` (memory `∝ m′`, §3f/R12). The turn
  remembers the whole accretion history, not only `m(v_peri)`.

Summary: the same FORM of law, exact Schwarzschild limit, robust
sign, plus a correction `∝ m′` that is the signature of
non-stationarity.

### R16. There is NO physical evaporative plunge inversion [CORRECTED]
(`physical_inversion.py`, `fig_vaidya_no_inversion_evaporation`)
It replaces the previous version (`../exploratory/VaidyaMetric/evaporation_inversion.py`,
`fig_vaidya_inversion_evaporation`), which was WRONG.

**History**: the dynamical correction `∝ m′` of R15 narrows the gap for
`m′<0`; a first analysis (linear model `m=1+μv`, r1=15, v1=40)
seemed to show an INVERSION at `μ_inv=−0.0577` (a second, evaporative
mechanism). **It is an ARTEFACT of the negative-mass regime.**

- With `m=1+μv` the mass becomes NEGATIVE for `v>1/|μ|`. At `μ_inv=−0.0577` the
  turn was at `m≈0.67` (positive), but the transversality anchor
  at `v1=40` was at `m=−1.31` (`m=0` at v=17.3): a boundary condition in an
  unphysical region (naked singularity) that contaminates the trajectory.
- The perturbative series `gap(μ)=gap₀+gap₁μ+…` was asymptotic **with a
  radius of convergence `R~|μ_inv|`**: the singularity IS precisely
  the entry into `m≤0`. The cube root (−0.053) extrapolates the physical
  small-μ response into the `m<0` regime.

**Control (two models, tracking m_min along the orbit)**:
- **Linear, 2D scan (v1, μ)**: EVERY case with `gap<0` has `m<0`; EVERY
  case with `m>0` everywhere has `gap>0`. No overlap ⟹ no
  physical inversion. Anchoring at `v<0` (m>1) the gap stays largely
  positive (+1÷+2.7).
- **Exponential `m=e^{−λv}` (m>0 ALWAYS)**: the gap narrows towards
  a plateau (~0.85) but **NEVER inverts**, not even at `λ=1` (m(40)=4e-18).

**Physical reason**: the inversion requires a strongly evaporated geometry
(small m) at the turn; the linear ramp cannot be "locally very
evaporated but globally m>0" (it keeps going below zero), and
the exponential, keeping m>0, evaporates only asymptotically while the
turn stays at m~O(1). In both cases the turn happens where m is
still substantial ⟹ the static ordering is preserved.

**Conclusion**: **physical evaporation NARROWS the gap but does not
invert it.** The ordering `r_min^t < r_min^τ` (t deeper, R13) is
robust for all physical dynamical evaporation. There is **ONLY ONE
inversion mechanism, the ROTATIONAL one** of Kerr/Thakurta-Kerr
(`ρ(r*)²=H_t0/H_τ0`, R12c, closed because stationary). The statement
"rotation is needed" of R13/R15 is CORRECT.

### R17. Boundary term per branch and time dependence of the optical metric (audit, September 2026)

- **Branch Hamiltonian with cost:** H = max_θ[p·ẋ − ℓ]. H_v = h − 1; H_τ = h(p_r + 1/Ê, J) − f/Ê,
  with ℓ_τ = dτ/dv = (f − ∂_{p_r}H_τ)/Ê. Free transversality: h = 1 for the v branch.
- **Euler identity with cost:** (J∂_J + p_r∂_{p_r})H = H + ℓ. With the degree-zero self-similarity (which holds for
  both Hamiltonians) one obtains on the shell ΘH = d(r p_r)/dλ − ℓ, hence
  **S_D = [r p_r] − Δv on the v branch and S_D = [r p_r] − Δτ on the τ branch** (the formula −λ held only for H_v).
  Counterexample to −λ for τ: m=1, r=6, Ê=7/5, J=0, p_r = −5√291/97 (ingoing), ℓ_τ = 7/5 + √291/15.
- **Degeneracies of the 11×11 system:** besides |J| = J_c also J = 0 (S = r³(r−2m)²D_E, rank 8) and
  E = 1; disc_r S ∝ Ê¹²J¹⁰m¹⁶ q(J²).
- **Letters:** V_k has residues at the points at infinity (dV₂: −1/√a4); U_2, U_3 on the sextic have
  residues −1/√s6 (U_2) and s5/(2 s6^{3/2}) (U_3), on the sheet at infinity chosen in the verifier
  (on the other sheet they change sign). Abel origin at a branch point, antiderivatives from the launch.
- **Non-autonomous optical metric:** in the threading with respect to ∂_v, h = dr²/f + r²dΩ²,
  a_rr = Ê²/[f²(Ê²−f)], a_φφ = Ê²r²/[f(Ê²−f)]; ∂_v a_rr ∝ m′(2Ê² − 3f), ∂_v a_φφ ∝ m′(Ê² − 2f),
  never zero together on f > 0. The non-Killing character of W alone would not suffice: a calibrated metric
  −f(t)dt² + [f(E₀²−f)/E₀²]δ has a_ij = δ_ij even though W is not conformal Killing.
- **Thrust:** it diverges at freezing (Ê² − |W|² → 0) if the drift numerator stays away from zero; not
  for |W| → 0 alone (ingoing radial rail: |a| → |m′|/(4rÊ³)).
Verification: Hamiltonians, Euler identity, boundary term, counterexample ℓ_τ, degeneracies and
residues in `paper1/verification/verify_paper1_core.py` (exact checks). Derivations in the text,
rechecked in SymPy but not in an archived verifier: the derivatives of the optical metric, the
non-conformal-Killing counterexample, the discussion of the thrust.

## 4. Next steps

1. ~~Penetration phenomenology: thresholds in J / windows in v~~ — DONE
   (§3b timing R7, §3c map R8). The evaporating sub-case "the
   penetrating orbit can escape" remains (orbits launched inside that come out when
   the horizon recedes; the v branch with r_min ≈ 2m(v_epoch) in R9 is
   the hint that the regime exists).
2. Brachistochrones in the self-similar sector: use the optical reduction R5
   for the massive t branch with the rail `Ê = −u·ξ`; trichotomy/dichotomy at `x±`
   in analogy with `J_c` and `J_t^c` of Kerr; connection with the map R8.
3. Quasi-static limit `μ → 0`: perturbative expansion around
   Schwarzschild (the analogue of the expansion in `a` for Kerr); the scaling
   `J_c ≈ (8÷10)√μ` of R8 as the first datum to be explained analytically.
4. Penetration map for the v branch (self-tuned, R9): a different threshold
   definition (marginal grazing instead of an interval?).
5. Paper: FLRW (exact basis) + Vaidya (dynamical phenomenology), with a unified
   Hamiltonian PMP section with transversality H=0.
