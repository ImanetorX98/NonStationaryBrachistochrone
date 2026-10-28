# Brachistochrones in NON-stationary spacetimes — FLRW and Vaidya

**Research kick-off note (handoff).** Goal: extend the formalism of the
*constrained worldline* (developed for Kerr, where it works thanks to
stationarity) to **dynamical** spacetimes, where no timelike Killing vector
exists. From here a new session must **derive the equations of motion and a
new variational principle** for the brachistochrone.

Originating project: `../KerrMetric/` (see `doranTau.md`, `doranT.md`,
`doranOffEq.md` for the complete stationary case).

---

## 0. What we carry over from the stationary case (Kerr)

The relativistic brachistochrone (Perlick 1991): between two **spatial** points, the
curve that minimises the travel time (proper `T_τ` or coordinate `T_t`)
for a particle at **fixed specific energy** `E = −u_t` (an ideal rail that
conserves the energy). Two equivalent formulations:

- **Optical/Randers reduction** (it compresses to a fixed 3D Finsler geometry on the
  quotient space). Elegant but **singular at the light surface** (ergosphere:
  `g_tt→0`), and it **exists only if the spacetime is stationary**.
- **Constrained worldline** (4D curve `x^μ(τ)` with a multiplier for `−u_t=E`).
  Regular in regular coordinates (Doran), integrable through the ergosphere.
  Lagrangian:
  ```
  L = (1 + μ E) Λ + μ g_{tν} ẋ^ν ,     Λ = √(−g_{αβ} ẋ^α ẋ^β) = dτ/dλ
  ```
  Rail force `f^μ` with `f·u=0` and `f_t=0` (⟹ `E` conserved because `∂_t` is Killing).

**Ingredients that depend on stationarity** (and must be replaced here):
- `∂_t` (timelike Killing field) → conserved energy `E = −u_t` (the rail constraint);
- `∂_φ` (axisymmetry Killing field) → conserved Fermat momentum `J`;
- the **quotient** by the flow of `∂_t` → canonical "space" and "travel time";
- the **optical reduction** to a fixed 3D Finsler geometry.

---

## 1. The general variational principle (invariant core)

The functional `∫dτ` is invariant; what must be reformulated is the **constraint**. Three known
routes, to be developed:

1. **Kovner–Perlick Fermat principle** (Kovner 1990, Perlick 1990, *Living Reviews* Perlick 2004):
   for light, the rays from a source to an **observer worldline** extremise
   the **arrival time** (the observer's proper time at reception). It holds in
   **arbitrary** spacetimes. → The **null** brachistochrone (limit `E→∞`) generalises
   fully.

2. **Zermelo navigation / non-autonomous Finsler geometry**: prescribe the local
   velocity field `v(x,t)` (instead of deriving it from the energy). Least-time path =
   Randers–Finsler geodesic, now **time-dependent**. It is the natural generalisation
   of the **massive** brachistochrone.

3. **Covariant constrained worldline**: the Lagrangian with a multiplier, but with the
   constraint referred to a **chosen observer congruence** `n^μ` (foliation
   `t=const`). Constraint: `−u·n = γ` prescribed (local speed relative to `n`).
   **Non-autonomous** equations; no static optical reduction; the 4D worldline is
   integrated directly.

**Table of what survives:**

| ingredient | stationary | non-stationary |
|---|---|---|
| `∫dτ` to be minimised | ✓ | ✓ (invariant) |
| conserved energy `E=−u_t` | ✓ | ✗ (needs a Killing field) |
| fixed 3D Finsler geometry | ✓ | ✗ (non-autonomous Finsler) |
| Fermat for light | ✓ | ✓ (Kovner–Perlick) |
| canonical "spatial points" | ✓ | ✗ (needs a choice of foliation) |

---

## 2. The FLRW case (expanding universe) — conformal structure

**Metric** (spatially flat k=0, generalisable to k=±1):
```
ds² = −dt² + a(t)² [ dr² + r² dΩ² ]
```
with `a(t)` the scale factor. **No timelike Killing field** (`a` depends on `t`). BUT:

**Conformal time** `η`, with `dt = a dη`:
```
ds² = a(η)² [ −dη² + dr² + r² dΩ² ] = a(η)² · η_μν^(Mink)
```
FLRW is **conformally flat**. `∂_η` is a **conformal Killing field** of the physical
metric (a genuine Killing field of the static conformal metric `ĝ = η_Mink`).

**Consequences to exploit:**
- **Light**: null geodesics are conformally invariant → the **null brachistochrone
  in FLRW = the null brachistochrone in Minkowski** (trivial) through the
  conformal transformation. A good consistency test.
- **Massive particles**: `Ê = −u_η` is conserved for the **conformal** metric `ĝ`, but the physical
  energy `E_phys` is **redshifted** (`E_phys ∝ 1/a`). The rail constraint can be
  placed on the **conformal quantity** `Ê` (conserved!) instead of on `E_phys`.

**Proposed variational principle (FLRW):** reformulate the brachistochrone in the
static conformal metric `ĝ` (where `∂_η` is a genuine Killing field), solve it there with the
stationary formalism, then **map it back** to cosmic time `t` through `dt=a dη`.
The local physical speed `v(x,t)` becomes time-dependent through `a(t)`.

**Open questions for the new session (FLRW):**
1. Derive the Lagrangian of the massive brachistochrone in FLRW, distinguishing the
   constraint on `Ê` (conformal, conserved) from that on `E_phys` (redshifted). Which is the
   physically correct "rail"?
2. The brachistochrone that minimises **cosmic time** `t`, the one that minimises
   **conformal time** `η`, and the one that minimises **proper time** `τ`: are they different? How?
3. Effect of the **expansion** (Hubble rate `H=ȧ/a`): is an orbit that starts "later"
   (larger universe) faster or slower? An analogue of θ-focusing?
4. Concrete cases of `a(t)`: radiation (`a∝t^{1/2}`), matter (`a∝t^{2/3}`),
   de Sitter (`a∝e^{Ht}`, which HAS an additional conformal Killing field). de Sitter is the
   most symmetric — a good starting point.
5. Boundary/horizon: is the **Hubble horizon** `r_H = 1/H` the dynamical analogue
   of the light surface? Does the brachistochrone penetrate it?

---

## 3. The Vaidya case (accreting/evaporating black hole) — null structure

**Metric** (ingoing Vaidya, advanced Eddington–Finkelstein coordinates):
```
ds² = −(1 − 2m(v)/r) dv² + 2 dv dr + r² dΩ²
```
with `v` = advanced time (`v=const` are ingoing light cones) and `m(v)` = **dynamical
mass**: `dm/dv > 0` accretion, `dm/dv < 0` evaporation (modelled Hawking
radiation). **No timelike Killing field** (if `m` depends on `v`).

**Useful structure:**
- The coordinate `v` is **regular across the horizon** (as Eddington–Finkelstein).
  No coordinate singularity at the dynamical Schwarzschild radius.
- **Apparent horizon**: `r = 2m(v)` (dynamical, it moves with `m`). Different
  from the event horizon (global, not local in Vaidya).
- Radial ingoing null geodesics: `v = const`. Outgoing ones: `dr/dv = (1−2m(v)/r)/2`.
- `g_vv = −(1−2m/r)` changes sign at the apparent horizon — the analogue of `g_tt` but
  **dependent on `v`**.

**Proposed variational principle (Vaidya):**
- **Light**: Kovner–Perlick Fermat principle — extremise the arrival time at an
  observer worldline. Direct and well posed.
- **Massive particles**: covariant constrained worldline. The constraint `−u_t=E` is NOT conserved
  (`m(v)` breaks the symmetry). Replace it by: local speed prescribed relative
  to the ingoing observers (`v=const` slicing) OR relative to the dynamical
  ZAMO/Doran observers. **Non-autonomous** equations in `v`.

**Open questions for the new session (Vaidya):**
1. Derive the equations of motion of the brachistochrone in Vaidya with the
   non-autonomous constrained worldline. Which constraint replaces `−u_t=E`?
2. **Moving apparent horizon**: if the black hole **accretes** (`dm/dv>0`),
   the horizon grows towards the particle — can it "swallow" it even if it was
   reflecting? If it **evaporates** (`dm/dv<0`), the horizon recedes — can the penetrating
   orbit "escape"? It is a **genuinely non-stationary** phenomenon, with no analogue in Kerr.
3. Does penetration depend on the **timing** of the arrival relative to `m(v)`: is there a
   temporal penetration window? A critical `v_c` beyond `J_c`?
4. Simple case: **linear Vaidya** `m(v)=m_0 + μ v` (constant accretion) or
   `m(v)=m_0(1 − v/v_evap)` (evaporation). Numerically integrable.
5. Quasi-static limit (`dm/dv → 0`): is the Schwarzschild case recovered? Perturbative
   expansion in `dm/dv` (an analogue of the expansion in `a` for Kerr?).

---

## 4. Recommended strategy for the new session

1. **Start from the covariant constrained worldline** (§1.3) — it is the right formalism,
   already validated on Kerr, and the only one that extends without stationarity.
2. **de Sitter FLRW** as the first case (maximal residual symmetry, explicit conformal
   Killing field) — derive the equations of motion, check the Minkowski limit, then generalise to `a(t)`.
3. **Linear Vaidya** as the second case (simple null structure, computable moving
   horizon) — derive the non-autonomous equations of motion, study penetration with a moving horizon.
4. **Expected strength**: **genuinely dynamical** phenomena with no stationary
   analogue — a horizon that chases or recedes, temporal penetration windows,
   redshift of the "rail". This is the new contribution.
5. **Tools**: SymPy for the metrics and equations of motion (as in KerrMetric), numerical integration
   of the non-autonomous worldlines, checks of the static limits.

---

## 4b. ANSWER TO §1 (July 2026): the general variational principle

Found and validated on FLRW, Vaidya, SdS and Thakurta(-Kerr). Three ingredients,
none of which requires symmetries:

1. a time function v (choice of foliation = what is minimised);
2. the rail constraint −u·W = E with W = ∂_v: LINEAR in the velocity;
3. the normalisation u·u = −1: quadratic.

Linear ∩ quadratic = **elliptic indicatrix** at every event ⇒ Randers/Zermelo
form `dv = drift + norm` (Legendre dual of the PMP Hamiltonians
with transversality H=0). It exists for ANY spacetime.

The "conserved quantity" is replaced by the **active constraint**: E is a
controlled invariant (Noether → control). Measurable cost: f·ξ of the
rail force (Kerr: f_t=0 for free; FLRW: f_η=aH; Vaidya: f_v ∝ m′).

Hierarchy of what symmetry buys:
- ∂_v Killing → autonomous F (Perlick 1991 recovered), H conserved;
- conformal Killing field / homothety → rigid skeleton (scalar index n(η,r), rigid wind
  in Thakurta-Kerr, self-similar quotient in linear Vaidya);
- ∂_φ Killing → Fermat momentum J (pure Noether): integrability of the ODE;
- none → F exists but "breathes" entirely; without ∂_φ one would need HJ (a PDE).

**Legitimacy of the constraint (status for the paper):**
1. well-posedness: μ exists, f·u = 0 identically (= covariant "ideal" constraint:
   preserving u·u = −1 is equivalent to f·u = 0) — proved;
2. existence of minimisers also in the non-autonomous case: FILIPPOV's theorem
   (compact convex indicatrix = our ellipse); the PMP does not
   require autonomy;
3. the null limit Ê→∞ = KOVNER-PERLICK Fermat principle, proved for
   arbitrary spacetimes (anchor in the literature);
4. equivalence with ZERMELO navigation (prescribed v(x,t)): the
   constraint is a symmetric derivation of the velocity field, not
   an ad hoc postulate.
Status of the choice of W: theorem (Killing) → canonical selection (conformal Killing field:
the only one that degenerates to Perlick + the null limit + minimal cost f·ξ) →
quasi-local selection (KODAMA, dynamical spherical case: -u·K = Kodama energy,
VaidyaResults 3f) → a convention to be declared (generic case;
classical analogue: a brachistochrone with a moving wire). Mapped physical
alternative: powered rail γ₀ (FLRWresults R6).

**Realisability of the control (not heuristic):**
- constructive: f := Du/dτ has f·u = 0 automatically (normalisation
  identity) and is of magnetic type (f^μ = F^{μν}u_ν with F
  antisymmetric): every admissible curve can be realised by an external
  force without proper work; maintaining -u·W = E is ONE linear
  condition on 3 force degrees of freedom ⇒ solvable away from the
  critical locus (v̄→0, horizon) — which coincides with the mapped critical
  surfaces; μ computed explicitly (finite).
- variational: VAKONOMIC/optimal-control reading (the rail is optimised,
  not a given rail simulated): PMP (no autonomy
  required) + existence of minimisers through FILIPPOV (compact convex
  indicatrix). Note: for constraints linear in the velocity the
  d'Alembert–Chetaev procedure differs from the vakonomic one in
  general; here it does not matter (a shape-optimisation problem) and the
  two collapse in the stationary limit (Perlick).

Details: `../VaidyaMetric/VaidyaResults.md`,
`../ThakurtaMetric/ThakurtaResults.md`, `../SdSMetric/SdSresults.md`,
`FLRWresults.md`.

## 5. Key references

- **Perlick, V.** (1991), *The brachistochrone problem in a stationary space-time*,
  J. Math. Phys. 32, 3148. [stationary basis]
- **Kovner, I.** (1990), *Fermat principle in arbitrary gravitational fields*,
  ApJ 351, 114. [non-stationary Fermat principle for light]
- **Perlick, V.** (2004), *Gravitational Lensing from a Spacetime Perspective*,
  Living Rev. Relativity 7, 9. [general Fermat principle, review]
- **Gibbons, Herdeiro, Warnick, Werner** (2009), *Stationary metrics and optical
  Zermelo–Randers–Finsler geometry*, PRD 79, 044022. [Zermelo/Randers]
- **Vaidya, P. C.** (1951), *The gravitational field of a radiating star*. [Vaidya metric]
- Originating Kerr project: `../KerrMetric/doranTau.md`, `doranT.md`, `doranOffEq.md`.

---

## 6. Context from the Kerr work (for continuity)

Stationary results already established (to be reused as limits or comparisons):
- **τ branch** (proper time, purely Riemannian): penetrates the ergosphere **only** at
  `J=J_c=a/E` (measure zero, prograde). Off the equator it reflects; θ-focusing towards
  ~65° = a manifestation of `â_θθ=F→0` (optical metric "soft" in θ).
- **t branch** (coordinate time, Randers): penetrates for a bounded **interval**
  `(J_t^{c,−}, J_t^{c,+})` that includes retrograde orbits (frame dragging forces co-rotation
  at `F=0`). Thresholds in **closed form** (elementary + saddle-node). Crossover of
  mechanism at `a* ≈ M/√2`.
- **Constrained worldline**: regular; the t branch reaches the ergosphere off the equator
  from every θ_0 with finite velocity/acceleration (`dr/dτ→0` at `F=0`).
- The `1/F` singularity is **intrinsic to the optical reduction** (light surface),
  not to the coordinates. No exact Carter constant off the equator.

The expectation: in FLRW/Vaidya the "light surface" becomes **dynamical**
(Hubble horizon / moving apparent horizon), and the penetration phenomenology
acquires a **temporal dimension** absent in the stationary case.
