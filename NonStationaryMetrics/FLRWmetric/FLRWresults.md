# t- and τ-brachistochrones in flat FLRW — results

**SymPy scripts** (all checks pass):
- `flrw_brachistochrone_sympy.py` — constrained worldline, equations of motion, rail force
- `flrw_optical_metric_sympy.py` — optical metric and clock λ(η) (R8)
- `flrw_hamiltonians_sympy.py` — PMP Hamiltonians of the t and τ branches (R9)
- `make_figures_flrw.py` — validation figures (`FLRWfigures/`)

Formalism: *covariantly constrained worldline* (from `../KerrMetric/doranTau.md`),
extended to the non-stationary case following the proposal of
`nonStationaryBrachi.md` §2: the constraint is placed on the conformal quantity.

## Setup

Flat FLRW metric (k=0), conformal time η (`dt = a dη`), motion in a plane:

```
ds² = a(η)² [ −dη² + dx² + dy² ]
```

`∂_η` is a **conformal** Killing field of the physical metric. Rail constraint:
`Ê = −u_η` held fixed (the conformal analogue of Kerr's `E = −u_t`).

Lagrangians (worldline, parameter λ, `Λ = dτ/dλ = a√(η′²−x′²−y′²)`):

```
τ branch:  L_τ = (1 + μÊ) Λ − μ a² η′
t branch:  L_t = a η′ + μ (Ê Λ − a² η′)        [T_t = ∫ dt = ∫ a dη]
```

## Results (checked symbolically)

### R1. The constraint fixes the kinematics; the rail is redshifted
From the constraint `−u_η = Ê`:

```
γ = dt/dτ = Ê/a          (E_phys = γ ∝ 1/a: redshift of the rail)
v(η) = √(1 − a²/Ê²)      (local physical speed relative to the comoving observers)
```

`v` depends **on time only** (homogeneity), not on position.

### R2. Degeneracy: the t, τ and η branches coincide as curves
`x, y` cyclic ⇒ `p_x, p_y` conserved ⇒ a **comoving straight line** for both
branches. On the solution:

```
dτ/dη = a²/Ê ,    dt/dη = a
```

Both functionals increase with the arrival η ⇒ the optimal curve is the
same (the straight line, that is the spatial 3D geodesic). The two branches differ only
in their multipliers:

```
μ_τ(η) = C/(Ê² v) − 1/Ê ,    μ_t(η) = C_t/(Ê² v)
```

**Answer to question 2 of the note (case k=0)**: the t-, η- and τ-brachistochrones are
the **same curve**; they differ only in the value of the travel time.
The argument uses homogeneity only ⇒ it holds for every FLRW (also k=±1, with
"straight line" → geodesic of the 3-geometry). In FLRW the physics of the brachistochrone
is entirely **temporal**; the shape dynamics of the orbit (Kerr) reappears only
when homogeneity is broken (→ Vaidya).

### R3. Rail force
Geodesic residual `f_μ = du_μ/dτ − ½∂_μ g_{αβ} u^α u^β` on the solution:

```
f_η = a′/a = aH     (conformal Hubble rate: the rail "sustains" Ê against the redshift)
f_x = −a′/(a v)
f·u = 0             (no proper work — exact)
du_η/dτ = 0         (constraint maintained)
```

The η equation holds identically on the ansatz (gauge identity):
the worldline is consistent. Comparison with Kerr: there `f_t = 0` (E is conserved by the
Killing field); here `f_η = aH ≠ 0` is the quantitative signature of non-stationarity.

### R4. Freezing: a temporal barrier, not a spatial one
Motion exists only for `a(t) < Ê`. At `a = Ê`: `v → 0` smoothly, and the particle
**freezes** at finite τ. It is the dynamical analogue of the light surface of
Kerr (`F=0`), but the barrier is in **time** instead of space.

### R5. de Sitter (closed form), `a(η) = −1/(Hη)`, η ∈ (−∞,0)
Start at `a=1` (η₀ = −1/H):

```
freezing:  a = Ê  at  t_f = ln(Ê)/H
Δx_max · H = ( √(Ê²−1) + asin(1/Ê) − π/2 ) / Ê
T_t(η₁)  = ln a(η₁) / H
T_τ(η₁)  = (−Hη₁ − 1)/(Ê H² η₁) ;   T_τ(freezing) = (Ê−1)/(ÊH)  finite
```

Null limit `Ê → ∞`: `Δx_max·H → 1` — the null brachistochrone reaches
**exactly** the comoving Hubble horizon (`r_H = 1/H` at η₀).
Conformal consistency with Kovner–Perlick checked (question 5 of the note: the
massive brachistochrone does NOT penetrate the horizon; only the null limit saturates it).

### R6. Conformal rail Ê against physical rail γ₀ (question 1 of the note)
Alternative constraint `−u·n = γ₀` (n = comoving observers):

- `γ₀` fixed ⇒ `v = √(1−1/γ₀²)` **constant**, no freezing;
  `−u_η = a γ₀` grows with the expansion ⇒ the rail must **pump energy**
  against the redshift (active rail).
- `Ê` fixed ⇒ passive/conservative rail (`f·u = 0`, fixed conformal component
  `f_η = aH`), with redshift `E_phys ∝ 1/a` and freezing.

They coincide only for `a = const` (Minkowski). The choice is physical, not
mathematical: `Ê` is the natural generalisation of the ideal Kerr rail
(no net work), `γ₀` models a powered rail.

> **[Clarification, September 2026]** f·u = 0 says only that the force is orthogonal to the four-velocity
> (a consequence of the normalisation), not that the rail is passive. The work against the
> selector, g(a,W), must be distinguished from the control needed to maintain the charge: with W conformal Killing the charge
> drifts at the rate a′/a (conformal time) along a free geodesic, and the rail must compensate it (Paper I, Lemma I.D).
> Freezing (Ê = |W|) and the stationary surface (W null) are also distinct boundaries.

### R7. Control limits
- Minkowski `a=1`: `v = √(Ê²−1)/Ê`, `γ = Ê` — standard. ✓
- Null limit: conformally invariant, Minkowski straight line. ✓

### R8. Optical metric from the conformal Killing field (`flrw_optical_metric_sympy.py`)
On the static `ĝ` (`g = a² ĝ`) the optical reduction exists and is **exact**:

- **Light**: `ĝ_{ηi} = 0` (no Randers term, no rotation) ⇒ optical metric =
  **comoving 3-geometry** `dl_k² = dχ² + S_k(χ)² dφ²`, index `n ≡ 1`.
- **Massive (rail Ê)**: Fermat with a purely temporal index `n(η) = 1/v(η)`.
  **Optical clock** `λ(η) = ∫ v dη` ⇒ `dλ = dl_k`: unit-speed motion in λ.
  Minimising `η_arrival` ⇔ minimising the comoving length ⇒ optical
  metric = comoving geometry for **all three branches** (t, τ, η —
  monotonicity). All the non-stationarity sits in the clock, none in the
  geometry.
- de Sitter: `λ(u) = [√(Ê²−1) − √(u²−1) − acos(1/Ê) + acos(1/u)]/(ÊH)`,
  `λ(η_f) = Δx_max` (identity checked). k=+1, numerical: brachistochrone =
  great circle (radiation `a=η`, argmin at ε=0).
- **Optical degeneration**: `n → ∞` at freezing `a → Ê`. In Kerr the
  optical reduction degenerates on a **spatial surface** (`F→0`,
  ergosphere); in FLRW it degenerates at a **cosmic instant**. The same
  `1/F` pathology, rotated from space to time.

### R9. Hamiltonians of the t and τ branches (`flrw_hamiltonians_sympy.py`)
Optimal-control formulation (Pontryagin), cosmic time t, control =
direction û, dynamics `dx/dt = (v/a)û`:

```
H_t   = (v/a)|p| − 1
H_τ   = (v/a)|p| − a/Ê          v(t) = √(1 − a²/Ê²)
```

- `dp/dt = −∂H/∂x = 0` (homogeneity) ⇒ fixed direction ⇒ comoving straight line:
  the trajectories follow from Hamilton's equations (checked against the closed forms
  R5 in de Sitter, 12 digits).
- `H` is **not** conserved: `dH/dt = ∂H/∂t ∝ ȧ` — the Hamiltonian version of the
  non-autonomous forcing (in Vaidya: `dp_v/dr ∝ m′`).
- Transversality (spatial arrival, free final instant): `H(t₁) = 0` ⇒
  `|p|_t = a₁/v₁`, `|p|_τ = a₁²/(Êv₁)` — the exact analogue of `p_v(r₁)=0`
  in Vaidya; the multiplier `μ(η)` of the worldline is the rescaled costate
  (`C = Êv(1+μÊ)` constant).

### R10. The exponential closed case is not de Sitter (audit, September 2026)

With a = e^{Ht} on closed slices, Δ_max = [√(Ê²−1) + arcsin(1/Ê) − π/2]/(ÊH) remains correct, but the
curvature scalar is **R = 12H² + 6e^{−2Ht}/R_c²**, not constant: this is not de Sitter. In true de Sitter,
in the closed foliation with radius H^{-1}cosh Ht, even a null ray from t₀ sweeps
**π/2 − arctan sinh Ht₀ < π**: the antipode is not reachable. The threshold is in HR_c at fixed Ê; at the
cut point uniqueness is lost, and global minimality beyond it. Figure `fig_flrw_curvature` regenerated with
the corrected title. The branches coincide in FLRW because the speeds and the clock rates depend only on
time (monotonicity), not because the functionals are proportional.

## Validation figures

In `FLRWfigures/` (`make_figures_flrw.py`): kinematics of the constraint and
freezing; comoving worldlines against the Hubble horizon (analytic Δx_max against
quadrature, error ~2e-6); direct variational minimum of T_t and T_τ on the comoving
straight line (validation of R2).

## Consequence for the research programme

FLRW is the "all time, no shape" case: useful as a laboratory for the
non-stationary constraint (rail redshift, freezing, force `f_η=aH`),
but the phenomenology of penetration with a non-trivial orbit requires
**spatial inhomogeneity + non-stationarity** ⇒ next step: **Vaidya**
(moving apparent horizon `r = 2m(v)`).
