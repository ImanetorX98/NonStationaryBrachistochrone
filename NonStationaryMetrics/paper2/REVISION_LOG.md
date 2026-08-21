# Paper II — CQG-116884, revision log

Internal working file. **Not** submitted to CQG. It is the source from which
`response_to_referees_CQG.tex` (the Author Response) is written, and it records
*why* each change was made so the two stay consistent.

Referee 1 = the technical report (`report_CQG_116884.pdf`), 12 major + 10 minor
comments + a required-reconstruction list.
Referee 2 = `CQGReviewerstask2026I.pdf`, bibliography + self-citation requests.

Status key: **DONE** / **WIP** / **TODO** / **DECLINED** (with argument).

---

## R1.4 — the $J=-J_c$ cusp/corner contradiction  — **DONE (analysis)**

> *"This is a decisive internal contradiction. Section 3.2, Figure 9, Table A1,
> and the later summary classify the retrograde marginal value $J=-J_c$ as a cusp
> at $r=2M$. Appendix A.1 instead says that at $|J|=J_c$ the double factor in the
> radicand is partially cancelled by a factor in the azimuthal numerator, so the
> retrograde marginal has finite $d\phi/dr$ and is a corner rather than a cusp."*

The referee also forbids resolving it by editing one sentence: the case must be
rederived from the Hamilton equations with one-sided expansions, signs, branch
orientation and numerical convergence tests.

### Outcome: **both** printed statements were wrong, and the corrected result is stronger.

Verification: `verification/verify_cusp_corner.wls` (22 checks, shape) and
`verification/verify_marginal_hamilton.wls` (23 checks, dynamics). All pass.

Notation: $r_e=2M$, $J_c=a/E$, $D_E(r)=(E^2-1)r+2M$, $\Delta=r^2-2Mr+a^2$,
$\bar P_e = 4M^2+2a^2+a^2/E^2$, frozen $A$, equatorial, $\tau$-branch.

**(a) The shape radical is even in $J$; the Hamilton equations are not.**
This single line is the whole asymmetry:

    ptilde_phi(r_e) = J - J_c                                        (linear in J)

`S(r)` enters only as `J^2`, so no shape argument can ever distinguish
$\pm J_c$. The conformal gravitomagnetic momentum can, and does.

**(b) Exact relation between the momentum and the shape radical:**

    p_r^2 = S(r) / [ D_E(r)^2 Delta(r)^2 ]        (exact, verified symbolically)

the prefactor being regular and nonzero at $r_e$ (value $1/(4a^4E^4M^2)$). So the
root structure of $p_r^2$ *is* that of $S$ — nothing is lost passing between them.

**(c) The bracket factorises exactly at the marginal value:**

    r*Delta - J_c^2 D_E(r) = (r - 2M)(r^2 + a^2/E^2)

so `S = r (r-2M)^2 D_E (r^2 + J_c^2)`: $r_e$ is a **double** root, for **both**
signs of $J$. Equivalently `S = r^4 f^2 w (r^2 + J_c^2)` with $w=\hat E^2-f>0$ —
a perfect square times a positive factor. **The "reality wall" $\sqrt{wf}$ is
therefore removed at $|J|=J_c$ for both signs**, not just the prograde one. The
old argument ("counter-rotation forbidden inside $r_e$") is not what excludes the
retrograde marginal, and is withdrawn.

**(d) What actually distinguishes the two signs is the rate.** With
$\mathrm{rhs}(r) = f/E - \tilde p_\varphi\varphi_0'$ the denominator of
$\dot r = R^2(\Delta/r^2)p_r/\mathrm{rhs}$,

    rhs(r_e) = - a (J - J_c) / Pbar_e      -> vanishes iff J = +J_c

so:

| case | $p_r^2$ at $r_e$ | rhs at $r_e$ | $\dot r$ near $r_e$ | $r_e$ attained? | shape at $r_e$ |
|---|---|---|---|---|---|
| $\|J\|>J_c$ | $>0$ | $\neq0$ | — | no | smooth periapsis at $r_{\min}>r_e$, exponent 2 |
| $\|J\|<J_c$ | simple zero | $\neq0$ | $\sim\sqrt{r-r_e}$ | **yes, finite parameter** | semicubical cusp, exponent 3/2 |
| $J=-J_c$ | double zero | $\neq0$ | $\sim C(r-r_e)$ | **no — asymptotic, log-divergent** | finite one-sided tangent, never attained |
| $J=+J_c$ | double zero | **simple zero** | finite $\neq0$ | **yes, and crosses** | regular crossing |

So the main text's "cusp at $J=-J_c$" is wrong (the exponent is not 3/2), and the
appendix's "reflects with a corner" is *half* wrong: the tangent is indeed finite,
but there is **no reflection** — $r_e$ is an asymptotic endpoint reached only as
the parameter diverges. It is a separatrix in the proper dynamical sense.

**(e) Closed forms** (all verified symbolically and against high-precision numerics):

    cusp amplitude, |J| < J_c:
        phi - phi_e = (2/3) K (r-r_e)^{3/2},
        K = (E J / a^2) sqrt( 2M / (a^2 - E^2 J^2) )
        -> K diverges as |J| -> J_c^- : the cusp family does NOT connect to the
           marginal case, which is why merging them was illegitimate.

    one-sided tangent at |J| = J_c:
        dphi/dr |_{r_e^+} = 2 M E J / ( a^2 sqrt(4M^2 + J_c^2) )

    retrograde marginal decay rate:
        dr/ds ~= C (r - r_e),   C = - sqrt(4M^2 + J_c^2) / (8 M^2)
        (depends on a and E only through J_c)
        elapsed parameter ~ |C|^{-1} ln[1/(r-r_e)]

    prograde marginal crossing rate:
        dr/ds |_{r_e} = - a^2 E^2 sqrt(4M^2+J_c^2) / ( 2M [ a^2(1+E^2) + 4E^2M^2 ] )

Numerics at $M=1,a=0.9,E=1.2$ ($J_c=0.75$): $C=-0.267000117$, matched by the ODE
to $10^{-6}$; the elapsed parameter gains $34.4956$ per four decades of
$r-r_e$, against $\ln(10^4)/|C| = 34.4955$. Prograde crossing rate
$-0.161020082$, closed form and numerics agreeing to $10^{-6}$. Sub-marginal
$J=-0.6$: exponent of $d\varphi/dr$ measured $0.4999999993$ against $1/2$;
arrival parameter finite, $5.31288649$, with the truncated tail scaling as
$\sqrt{\epsilon}$ (ratio $9999.99986$ against $10^4$).

**(f) Contrast with the $t$-branch, now stated explicitly.** Both retrograde
marginals are asymptotic endpoints, but of different kinds:

- $t$-branch at $J_c^-$: double root at $r_*=3.514>r_e$, and $\mathcal K$ does
  *not* vanish there, so $d\varphi/dr$ diverges — **infinite winding**
  (logarithmic $\varphi$-barrier). Already correctly described in the appendix.
- $\tau$-branch at $-J_c$: double root at $r_e$, and $\mathcal K_\tau$ vanishes
  linearly there too, so $d\varphi/dr$ stays finite — **finite total winding**,
  the orbit converging to a single point of the circle $r=r_e$ with a definite
  limiting tangent.

**(g) Aside, worth keeping.** Near $r_e$ the retrograde marginal has
$d\varphi/dr<0$ and $\dot r<0$, hence $\dot\varphi>0$: it is locally
*co-rotating* despite $J<0$. So the discarded "counter-rotation" story was wrong
in its detail as well as in its role.

### Consequent edits (all flow from the table above)
- [x] §3.2 rewritten as `sec:trichotomy`: four-case classification stated as
      `prop:classification`, derived from the Hamilton equations
- [x] §3.2 retrograde paragraph: counter-rotation argument withdrawn, replaced
      by the rate degeneracy `eq:T-at-re`
- [x] Table `tab:penetration` restructured — the old "cusp bounce" column could
      not hold $-J_c$; now keyed on root multiplicity and attainment
- [x] appendix `app:explicit` passage at $|J|=J_c$ (the "corner" sentence)
- [x] "Conformal trichotomy" paragraph -> "Conformal form of the classification"
- [x] caption of `fig_master_penetration_taut`
- [x] caption of `fig_thakurta_cuspide_ergosfera`
- [x] caption of `fig_atlas_tau`
- [x] abstract and introduction
- [ ] conclusions §4
- [ ] **figures themselves regenerated** — the three above are still the old
      renderings; `fig_atlas_tau` in particular draws a cusp in the $-J_c$ panel
      and must be redrawn as an asymptotic approach truncated at a parameter
      horizon. Generators: `ThakurtaMetric/cuspide_ergosfera.py`,
      `KerrMetric/fig_master_penetration_taut.py` (**note: `KerrMetric/` is
      untracked — this is exactly referee 1's major 10 about scripts living in
      sibling directories outside the archive**), atlas generator TBD.
- [ ] title still reads "ergosphere trichotomy" — see decision note below

### Decision required: the title
The title advertises *"ergosphere trichotomy"*. The trichotomy is gone (it is a
four-case law) and referee 1 major 5 asks that the interior trichotomy be
replaced as the flagship result by the exterior separatrix at
$r_d\simeq3.5139M$. Leaving the title would contradict the reconstruction.
Suggested: *"Brachistochrones in conformal Kerr spacetimes: control domain,
exterior separatrices, and adiabatic response"*. **Author decision.**

---

## R1.1 — exterior/contact/continuation discipline — **WIP**

New subsection `sec:domain`, "The control domain, and the status of continued
formulae", inserted after the Hamiltonians in §2. Contains:

- `lem:compact`, with proof: the rail slice
  $\{g(u,u)=-1,\ -g(u,W)=\hat E\}$ is nonempty and compact iff $W$ is timelike
  and $\hat E\ge|W|$; a single point iff $\hat E=|W|$; noncompact if $W$ is null
  or spacelike. This is the referee's own argument, made into a lemma.
- $g(W,W)=-A^2(1-2M/r)$ for $W=\partial_\eta$, hence the domain is $r>2M$, and
  Paper I's two scalar inequalities read $r>2M$ and $\bar v^2>0$ here.
- The explicit statement that Doran coordinates do not repair this, and that a
  Doran- or ZAMO-anchored selector would be a different problem.
- `prot:status`: every curve is (E) exterior extremal, (C) limiting contact, or
  (A) analytic continuation with no optimality claimed and a distinct line style
  (this also discharges referee 1 M6).
- A `remark` discharging R1.2: PMP is necessary not sufficient, so
  "controlled-rail extremal" is the default term.
- `prot:endpoint`: fixed launch event, fixed spatial target, free arrival clock,
  branch-dependent shooting; comparisons at the relabelled costate $J_t=AJ_\eta$
  (referee 1 M1 and part of major 3).

Still to do: propagate the vocabulary through the remaining sections and
figure captions, and relabel type-(A) curves in the figures themselves.

---

## Remaining referee 1 items — **TODO**

R1.1 exterior/contact/continuation discipline · R1.2 extremal vs minimizer ·
R1.3 endpoint protocol + costate relabelling · R1.5 exterior separatrix as the
flagship · R1.6 tiered fixed-endpoint status + `cap_full.py` · R1.7 adiabatic
scope + one traceable dataset · R1.8 neutral higher-genus terminology ·
R1.9 physical dictionary · R1.10 dependency-complete archive · R1.11 Paper I
dependencies · R1.12 positioning vs literature · M1–M10.

## Referee 2 — **TODO**

Accept: Perlick 1991; Giannoni–Masiello–Piccione 1997; Giannoni–Piccione–Tausk
1999 (Morse theory) and 2002; Taş 2026; Haws–Kiser 1995. The Morse-theory item
connects to R1.2 and to the conjugate-point/Maxwell material.

Decline, with argument, the requested characterisation of the results as "a
trivial application of theorem T2" of an unrefereed preprint, and the demand for
"experimental validation" of a mathematical result. See the Author Response.
