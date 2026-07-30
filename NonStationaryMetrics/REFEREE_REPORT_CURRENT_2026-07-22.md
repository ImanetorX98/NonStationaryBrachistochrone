# Technical referee report on the current manuscript

## Purpose and audit snapshot

This document is an operational handoff for revising the current version of the
paper. It is intended to be sufficiently self-contained that another model or a
human collaborator can understand:

- what the paper is trying to prove;
- which parts are already correct and should be preserved;
- which claims are presently stronger than the derivations support;
- what must be derived, recomputed, renamed, or downgraded;
- which tests should pass before an objection is considered closed.

Audit date: **2026-07-22**  
Repository: `NonStationaryMetrics`  
Branch: `separatrice-adiabatica-vaidya`  
Audited commit: `61967a7`  
Main manuscript: `paper/main.tex`  
PRD mirror: `paper/main_prd_revtex.tex`

Line numbers below refer to the audited commit and may move after editing. Search
for the quoted equations or phrases if the line number no longer matches.

No scientific source file was modified during the audit. The only new file is this
report.

---

## Executive verdict

**Recommendation: major revision.**

The project has a strong mathematical core and is not a project to abandon. Several
earlier objections have been addressed correctly, the main algebraic reductions are
reproducible, and important distinctions between physical separatrices and algebraic
degenerations are now present in parts of the manuscript.

The paper is nevertheless **not yet rigorous enough to support its principal
non-stationary theorem or its strongest physical conclusions**. The central problem
is not a small algebraic error: the first-order adiabatic correction is obtained by
differentiating a frozen stationary family, while the paper itself correctly states
that the true non-autonomous extremal is off that frozen zero-Hamiltonian shell in the
interior. The current validations reconstruct the same assumed frozen-source formula
and therefore do not close this gap.

The manuscript can become rigorous, but the following two issues are blocking and
must be handled before treating the paper as complete:

1. derive the first variation of the **actual non-autonomous Hamiltonian boundary
   value problem**, including endpoint transversality and the variation of the
   trajectory;
2. distinguish the conserved Pontryagin costate `p_phi` from the mechanical angular
   momentum used in the conformal-scaling argument, and derive their relation rather
   than identifying them by notation.

Until those steps are done, the exact genus-two reduction should be advertised as an
exact closed form for an explicitly defined **frozen-source adiabatic candidate**, not
as a proved first-order correction to the true non-stationary optimum.

### Severity legend

- **BLOCKER**: affects the main theorem or the physical interpretation of the main
  result.
- **MAJOR**: a theorem/claim is false, conditional, internally contradictory, or not
  established by the cited computation.
- **MINOR**: notation, wording, reproducibility, or layout issue that does not by
  itself invalidate the mathematical core.
- **VERIFIED / PRESERVE**: checked and should not be casually rewritten.

---

## What the manuscript is trying to establish

The paper develops a controlled relativistic brachistochrone in non-stationary
spacetimes. Its intended logical chain is:

1. impose the actively maintained rail constraint
   `-u.W = Ehat` relative to a selector field `W`;
2. derive a velocity indicatrix and the corresponding branch Hamiltonians using the
   Pontryagin maximum principle;
3. apply the construction to FLRW, ingoing/outgoing Vaidya, and conformal Kerr
   (Thakurta--Kerr);
4. derive stationary/frozen turning laws, penetration thresholds, cusp behaviour,
   and separatrix degenerations;
5. treat a slowly varying mass or conformal factor by differentiating the frozen
   shape integral;
6. reduce the resulting `S^{-3/2}` differential to a finite collection of Abelian
   integrals and depth-two iterated integrals on a genus-two curve;
7. interpret those corrections as a universal first-order adiabatic theorem for a
   breathing indicatrix.

Steps 4 and much of step 6 are currently the strongest parts. The logical gap lies
between steps 3 and 5: it has not been shown that differentiating the frozen shape
family produces the first variation of the true non-autonomous optimal-control
trajectory.

---

# Part I — Blocking mathematical issues

## 1. BLOCKER — The adiabatic correction is not derived from the true non-autonomous extremal

### Current statements

The foundational section correctly states in `paper/main.tex`, around lines 172--185,
that:

- `J = p_phi` is conserved by axial symmetry;
- free-arrival transversality imposes `H(s_f)=0` only at the endpoint;
- for a non-autonomous system, `dH/ds = partial_s H`;
- therefore `H` is generally nonzero in the interior;
- the globally constrained quantity is the extended Hamiltonian
  `H_ext = p_s + H = 0`.

This is a valuable correction and should be preserved.

However, the Vaidya adiabatic construction begins around lines 516--550 by taking a
frozen Schwarzschild shape integrand `F(r;m,E,J)` and writing a correction proportional
to `partial_m F` weighted by the frozen clock. The Thakurta--Kerr construction around
lines 910--945 similarly uses a frozen Kerr shape integrand and the Euler source in
`E` and `J`. The unified theorem around lines 1515--1585 then elevates this operation to
a rigorous first-order theorem.

The frozen orbit used to define `F` is obtained after solving the stationary
zero-Hamiltonian relation. The true non-autonomous orbit explicitly does **not** obey
that relation throughout its interior. Therefore the following implication has not
been proved:

```text
first variation of the true non-autonomous extremal
    = parameter derivative of the frozen H=0 shape relation.
```

### Why the exact reductions do not solve this problem

The identities of the form

```text
partial_lambda F = N_lambda / S^(3/2)
                 = d/dr(A/sqrt(S)) + sum_k c_k r^k/sqrt(S)
```

are important and, where checked, correct. They prove an exact reduction of the
chosen differential. They do not prove that this differential is the forcing term in
the Jacobi/variational equation of the actual optimal-control boundary value problem.

Similarly, the integration-by-parts and shuffle identities prove the closed-form
assembly of the assumed source. They do not supply the missing dynamical derivation.

### Existing validation is circular

`ThakurtaMetric/fig_phi_closed_validation.py`, around lines 70--76, defines

```python
dphi_full = -Ehat*cumulative_trapezoid(dEF(rg)*eta, rg)
```

and then compares this quantity with an integration-by-parts rearrangement of the
same integral. The plot labels the former as `numerical (ODE)`, but it is not obtained
by integrating the true time-dependent canonical equations. It is the same assumed
quadrature represented in a different form.

This script is useful as a check of the closed-form decomposition, but its labels and
the corresponding manuscript wording must not describe it as a validation against
the non-autonomous solution.

`ThakurtaMetric/rail_conservation.py` has the same structural limitation. Its
docstring, around lines 21--31, assumes that the brachistochrone is
"instantaneously an on-shell Kerr orbit." The function `pr_onshell`, around lines
138--150, solves `H2=0` at each radius by construction. The script then demonstrates
that the prescribed rail energy is constant on the constructed curve. It does not
show that the constructed curve satisfies the complete non-autonomous costate system
or its endpoint conditions.

### Independent diagnostic

As a diagnostic, the same branch Hamiltonian was integrated as a genuinely
non-autonomous canonical system on a smooth incoming arc, and the first-order
coefficient was estimated by a symmetric finite difference. For the representative
choice

```text
M=1, a=0.9, Ehat=1.4, J=6, incoming arc r=12 -> 5,
```

the canonical-flow coefficient was approximately `+1.244`, while the frozen formulas
gave approximately `-7.146` for the energy-only source and `-20.074` for the current
energy-plus-angular-momentum source.

This diagnostic is not itself the final physical free-arrival BVP: the launch was
anchored on the stationary `H=0` shell, in the same spirit as the present plotting
scripts, rather than re-shooting the full endpoint-transversality problem. It should
therefore not be promoted as a replacement theorem. It is nevertheless a strong
falsification of the claim that the existing quadrature has already been validated
as the generic first variation of the non-autonomous canonical flow.

### Required derivation

Introduce an explicit slow parameter, for example

```text
lambda(s) = lambda_0 + epsilon lambda_1(s) + O(epsilon^2),
z(s;epsilon) = z_0(s) + epsilon delta z(s) + O(epsilon^2),
```

where `z` contains all state and costate variables, including the relevant clock
coordinate and its conjugate momentum. Starting from the extended Hamiltonian, derive

```text
delta dot z
  = J_sympl Hess_z(H_0) delta z
  + J_sympl partial_lambda grad_z(H_0) lambda_1(s)
```

plus any terms caused by the chosen parametrization.

The boundary conditions must be linearized as well:

- fixed initial event/state constraints;
- target spatial point, target worldline, or other actual terminal manifold;
- free-arrival-time condition;
- endpoint transversality `H(s_f)=0`;
- variation `delta s_f` of the terminal parameter/time;
- any normalization or shooting condition used to select the normal extremal.

Then derive the first-order correction to the observable `phi(r)` or to the endpoint
angle. There are two possible outcomes:

1. the current frozen-source term reappears, accompanied by a proof that all
   trajectory/costate corrections cancel or are absorbed into boundary terms; or
2. additional Green-function/Jacobi-field or boundary terms survive, in which case
   the genus-two reduction must be applied to the corrected source.

### Required numerical validation

Implement two independent computations for the same boundary data:

1. solve the full non-autonomous BVP for `+epsilon` and `-epsilon`, re-shooting the
   endpoint/transversality conditions in each case;
2. solve the derived linear variational BVP once at `epsilon=0`.

Compare

```text
[Observable(+epsilon)-Observable(-epsilon)]/(2 epsilon)
```

with the variational solution. The error should scale as `O(epsilon^2)` for the
central difference and remain stable under tighter ODE/BVP tolerances and mesh
refinement.

Only after this comparison succeeds should the special-function closed form be
called the first-order adiabatic correction to the physical optimum.

### Acceptance criteria

- [ ] A derivation begins from the full extended Hamiltonian, not from the frozen
      shape relation alone.
- [ ] State, costate, terminal-time, and transversality variations are included.
- [ ] The relation between the derived source and `partial_lambda F` is explicit.
- [ ] A full non-autonomous BVP, not a reconstructed quadrature, is used as the
      numerical reference.
- [ ] Symmetric finite differences and the variational solution agree with the
      expected convergence order.
- [ ] Figure labels distinguish "quadrature identity check" from "non-autonomous
      BVP validation."
- [ ] Until all previous items pass, abstract, theorem, conclusions, and captions
      use conditional language.

### What should be preserved

Do not discard the exact polynomial reductions, the Abelian integral basis, or the
integration-by-parts assembly. They remain valuable. The required change is to put a
proved dynamical source in front of that machinery.

---

## 2. BLOCKER — `J` conflates a Pontryagin costate with mechanical angular momentum

### Current statements

Around line 172 of `paper/main.tex`, the manuscript defines

```text
J = p_phi,
```

which is conserved because the branch Hamiltonian is independent of `phi`. This is a
canonical/Pontryagin statement.

Around lines 418--432 and again in the Thakurta--Kerr adiabatic section, the manuscript
states that under `g_TK = A^2 g_Kerr` both Kerr charges have conformal weight `A^-1`,
so that

```text
E_eff = Ehat/A,
J_eff = J/A,
```

and consequently uses the source `(E partial_E + J partial_J)F`.

### Why the identification is not automatic

The mechanical spacetime angular momentum normally means a quantity such as

```text
L_mech = g(u, partial_phi).
```

The conserved `p_phi` in the time-optimal control Hamiltonian is the costate conjugate
to the coordinate `phi`. For a forced rail, these are not automatically identical.
The constraint force may exert a torque about the symmetry axis even though the
optimal-control Hamiltonian is axially symmetric and its canonical costate is
conserved.

Thus the following must be derived, not assumed:

```text
p_phi  <->  L_mech,
```

including any branch-dependent factors, clock parametrization factors, or Legendre
map.

The current scripts expose the ambiguity. `rail_conservation.py` says `J` is conserved
and passes a fixed `J` into the frozen Hamiltonian, whereas the paper's revised
adiabatic source treats the effective Kerr angular momentum as `J/A`.

### Required correction

Use separate notation immediately, for example:

```text
P_phi       = Pontryagin costate, conserved by axial symmetry;
L_mech      = g_TK(u, partial_phi);
ell_eff     = Kerr mechanical/angular parameter used by F.
```

Derive the constrained Lagrangian or Hamiltonian Legendre map and express `ell_eff`
in terms of `P_phi`, `Ehat`, `A`, and the branch data. Then differentiate the actual
relation when constructing the adiabatic source.

Possible outcomes include:

- `ell_eff=P_phi/A`, validating the existing Euler source;
- `ell_eff=P_phi` in the shape Hamiltonian, eliminating the `J partial_J` term;
- a more complicated relation producing additional terms.

The result cannot be chosen by conformal dimensional reasoning alone.

### Acceptance criteria

- [ ] Canonical costate and mechanical angular momentum have distinct symbols.
- [ ] Their relation is derived from the actual controlled problem.
- [ ] The conformal weight of each quantity follows from that derivation.
- [ ] The Thakurta--Kerr source is recomputed from the derived relation.
- [ ] `rail_conservation.py`, the manuscript, and the reduction scripts use the same
      convention.
- [ ] Any numerical result depending on the `J partial_J` term is regenerated.

---

# Part II — Major internal contradictions and physical interpretations

## 3. MAJOR — Algebraic degeneration, physical separatrix, ergospheric threshold, and penetration threshold are mixed

### What is already correct

The manuscript now contains the crucial observation, around lines 621--626, that for
the Schwarzschild/Vaidya frozen curve with `E=1.4`, the value

```text
|J| approximately 7.0266,  r_d approximately -3.3637
```

is a **negative-radius algebraic genus degeneration**, not a physical external
separatrix. This is confirmed by `separatrix_classification.py`.

The same script correctly finds representative Thakurta--Kerr/Kerr `t`-branch roots,
including:

- a physical external retrograde double root near
  `J=-8.05352`, `r_d=3.51391`;
- a prograde/ergospheric double root near `J=+2.936`, `r_d=1.512`;
- large-`|J|` negative-radius algebraic degenerations.

This classification should be preserved.

### Remaining contradictions

Elsewhere the same Vaidya value `J/m approximately 7.0266` is still called a
separatrix without qualification, including around lines 525 and 610--615. More
seriously, the unified theorem and appendix around lines 1573--1590 and 2290 onward
return to calling it a moving physical separatrix and infer:

```text
accretion -> capture,
evaporation -> escape.
```

That physical interpretation does not follow from a double root at negative radius.
The algebraic cancellation of a triple pole can remain valid, but it is not evidence
for an exterior capture/escape transition.

The symbol `J_c` is also used for several inequivalent objects:

1. Vaidya dynamical penetration threshold `J_c(v_0)`;
2. Vaidya/Schwarzschild negative-radius degeneration `J_deg`;
3. Thakurta--Kerr `tau` ergosphere threshold, approximately `a/E`;
4. Thakurta--Kerr/Kerr physical external `t`-branch separatrix;
5. further negative-radius roots in Appendix B.

### Required correction

Adopt a global notation table and use it consistently:

```text
J_ergo                 ergospheric/grazing threshold;
J_deg                  algebraic discriminant root with nonphysical r_d;
J_sep^phys             external physical separatrix;
J_pen(v_0)             dynamical Vaidya penetration threshold.
```

For every double root reported in the paper or a figure, state:

- branch (`t`, `tau`, advanced/retarded clock);
- parameter values;
- double-root radius;
- horizon/ergosphere radius;
- whether `r_d` lies in the physical domain;
- whether it changes the number of accessible real turning points;
- whether the result is algebraic only or has a capture/escape interpretation.

### Acceptance criteria

- [ ] No negative-radius double root is called a physical separatrix.
- [ ] The value `7.0266` is consistently labelled `J_deg/m` in the Vaidya example.
- [ ] Capture/escape language is used only after a physical-domain phase-portrait or
      BVP analysis.
- [ ] `J_c` no longer denotes more than one mathematical object.
- [ ] Abstract, unified theorem, appendices, captions, and conclusions agree with the
      classification script.

---

## 4. MAJOR — The physical Vaidya evaporation claim is not established

### Correct convention now present

The paper correctly states around lines 490--497 that physical evaporation is
described by outgoing Vaidya,

```text
ds^2 = -f du^2 - 2 du dr,
```

with retarded time `u`, rather than by simply taking `m'(v)<0` in the ingoing advanced
metric. This convention should be kept.

### Unsupported conclusion

The abstract and body state that there is no physical evaporative inversion. The
main numerical script used for this conclusion,
`VaidyaMetric/inversione_fisica.py`, still derives and integrates the **ingoing**
metric equations using `m(v)`. It scans positive linear/exponential decreasing mass
profiles, but does not derive or integrate the outgoing `u`-Hamiltonian/PMP boundary
value problem.

It can support only the narrower statement:

> In the tested ingoing toy models, keeping the mass positive removes the apparent
> inversion seen after extrapolating a linear mass law into `m<0`.

It does not prove the physical outgoing-Vaidya statement.

### Required correction

Derive the outgoing indicatrix and both branch Hamiltonians directly from

```text
ds^2 = -f(u,r) du^2 - 2 du dr + r^2 dOmega^2.
```

Then solve the same endpoint protocol used for the ingoing comparison and scan a
well-specified class of positive decreasing mass functions. State whether the result
is:

- an analytic ordering theorem;
- a validated numerical observation on an explicit compact parameter domain; or
- a conjecture suggested by examples.

Do not use "absence" or "never" without one of the first two levels of evidence.

### Concrete sign error in the asymmetry paragraph

The current definitions are

```text
delta_accr / q    = A+B,
delta_evap / (-q) = A-B,
q=|m_dot|.
```

Therefore the physical corrections are

```text
delta_accr =  q(A+B),
delta_evap = -q(A-B).
```

It follows algebraically that

```text
delta_accr - delta_evap = 2q A,
delta_accr + delta_evap = 2q B.
```

The manuscript around lines 701--720 and
`VaidyaMetric/vaidya_asymmetry.py` around lines 6--9 and 38--51 instead describe the
difference as `2qB`. That compares the coefficients after dividing by mass rates of
opposite sign, not the two physical corrections at equal `|m_dot|`.

For the script's current numbers,

```text
A approximately 8.7276,
B approximately 6.1127,
delta_accr/q approximately +14.840,
delta_evap/q approximately -2.615.
```

Hence the physical difference is approximately `17.455 q = 2qA`, while the physical
sum/deviation from naive antisymmetry is approximately `12.225 q = 2qB`.

The horizon term `B` can still be interpreted as controlling the **failure of exact
antisymmetry**, but not the ordinary physical difference as currently worded.

### Acceptance criteria

- [ ] The sign algebra is corrected in manuscript, script comments, printed output,
      and figures.
- [ ] "Difference", "sum", and "asymmetry" are explicitly defined.
- [ ] The outgoing Vaidya Hamiltonian is derived rather than inferred by a clock sign
      replacement alone.
- [ ] The no-inversion statement is proved or restricted to the numerical domain
      actually tested.
- [ ] Ingoing decreasing-mass toy results are not presented as outgoing evaporation.

---

## 5. MAJOR — Fixed-endpoint no-inversion is still conditional

### Current rigorous status at commit `61967a7`

The reduction of the fixed-endpoint problem into Lemmas A and B is useful and should
be preserved.

**Lemma A is proved.** At matched turning radius, the factorization

```text
N_G = (r-r_min) [4 r Delta/(r_min-2M)] W_G(r)
```

with `W_G` linear and positive proves the pointwise ordering of the branch integrands
in the stated exterior domain. The symbolic proof in
`no_inversion_reduction.py` is a genuine algebraic result.

**Lemma B is not proved on the entire claimed scattering regime.** The current status,
documented more honestly in `TODO_lemmaB_closure.md`, is:

- monotonicity is proved in a grazing regime;
- in the frozen Schwarzschild reduction it has been extended by an exact closed form
  for `Phi'` and an elementary estimate to `V_min >= V_0/4`;
- a single-crossing criterion labelled `(DAGGER)` has been derived;
- `(DAGGER)` has been checked numerically to high precision;
- the interval between the deflection peak `r_pk` and the elementary quarter-regime
  is still open;
- for frozen Kerr, the corresponding monotonicity remains numerical;
- the transition location is genuinely transcendental, so a simple polynomial
  certificate is not expected.

This is meaningful progress, but it does not justify the sentence in
`paper/main.tex` around lines 1470--1480 saying that fixed-endpoint no-inversion "is a
theorem in the scattering regime" while simultaneously stating that Lemma B is only
verified beyond the proved sub-regime.

The protocol table and captions around lines 1495--1510 also say "never" and attribute
the conformal case to the pointwise ratio `n_t/n_tau=E/f`. The preceding text correctly
states that this pointwise ratio alone does not order the minima of two different
functionals. The table therefore contradicts the derivation.

### Honest resolution options

Choose one of the following.

#### Option A — Full analytic proof

Prove `(DAGGER)` or an equivalent single-crossing/monotonicity result throughout the
remaining domain, then extend or separately prove the Kerr case.

#### Option B — Computer-assisted theorem on a compact domain

Specify compact ranges of `(E,a,r_0,r_min)` and use interval arithmetic with outward
rounding to certify the sign, including near the tight threshold. Document the
subdivision strategy and provide a reproducible certificate. This gives a rigorous
theorem on the certified domain, not for all parameters.

#### Option C — Conditional theorem plus numerical evidence

State:

> Assuming Lemma B (monotonicity/single crossing), fixed-endpoint no-inversion follows.
> Lemma A is analytic; Lemma B is analytic in the stated sub-regime and numerically
> supported elsewhere in the explored domain.

This is fully acceptable if claims, abstract, captions, table, and conclusions all
use the same status.

### Acceptance criteria

- [ ] The manuscript does not call Lemma B proved outside its proved domain.
- [ ] "Never" is removed unless supported by a global proof.
- [ ] The table no longer invokes a pointwise index ratio as a proof after the text
      explains why that implication is invalid.
- [ ] Parameter ranges and numerical precision are stated wherever evidence is
      numerical.
- [ ] If interval arithmetic is used, the certificate is reproducible from the repo.
- [ ] `TODO_lemmaB_closure.md`, `REFEREE_RESPONSE.md`, abstract, theorem, and conclusion
      report the same logical status.

---

# Part III — Foundational control-theory issues

## 6. MAJOR — Compact oval, convex control set, existence, and normality are conflated

### Boundary of an ellipse versus filled ellipse

Equation `(ellipse)` parametrizes the **boundary** of an ellipse by an angle. The
boundary is a compact, strictly convex curve in the geometric sense, but it is not a
convex set: the line segment joining two boundary velocities generally lies in the
interior and is not included.

Thus the statement around lines 160--166 that the admissible-velocity set "is already
compact and convex (the oval)" is false as written.

Two mathematically different models are possible:

1. fixed rail speed/mass-shell control: admissible velocities are the boundary oval;
2. relaxed variable-speed control: admissible velocities fill the ellipse.

If the second model is intended, it must be defined and one must prove that a
time-optimal solution saturates the boundary and that no relaxation gap changes the
physical problem. If the first is intended, do not invoke convexity of the velocity
image as though the full Filippov hypothesis were automatic.

### Unique support maximizer

Strict convexity gives a unique support point for a **nonzero** covector. At `p=0`,
every velocity maximizes `p.v`. PMP nontriviality prevents the entire multiplier pair
from vanishing, but it does not by itself justify silently dividing by or normalizing
the spatial costate at every point. Normality and nonvanishing of the relevant
costate component should be stated or proved in the domain where the closed branch
Hamiltonian is used.

### PMP is necessary, not an existence/global-optimality theorem

The manuscript calls the problem "well-posed" around lines 95--99. Pontryagin's
maximum principle supplies necessary conditions for an optimum under appropriate
hypotheses; it does not alone prove that an optimum exists, that the chosen extremal
is globally minimizing, or that it is unique.

A full well-posedness statement would require at least:

- precise admissible controls and regularity;
- existence of controlled trajectories on the interval/domain;
- compactness and the correct convexity/closure hypotheses, or a relaxation argument;
- reachability of the terminal set;
- lower semicontinuity/compactness for a minimizing sequence;
- treatment of abnormal extremals;
- a reason the selected extremal is a minimum rather than another stationary
  trajectory.

The standard Filippov existence discussion explicitly requires hypotheses beyond
compactness of an angular parameter; see
`https://liberzon.csl.illinois.edu/teaching/cvoc/node89.html`.

### Endpoint formulation

The introduction calls the problem the fastest rail worldline "between two events."
If both spacetime events are fixed, the coordinate-time separation is already fixed
and a free-arrival-time objective is not being minimized. The actual terminal data
appear to be an initial event plus a target spatial point, target orbit/worldline, or
terminal manifold with free arrival time. State this precisely and use the same
protocol in derivations and numerics.

### Freezing surface statement

Around lines 147--150 the collapse of the velocity oval is called "the
turning-point/separatrix locus." A radial turning point does not generally require the
entire velocity indicatrix to collapse; a trajectory can have zero radial velocity
and nonzero tangential velocity. A separatrix is also a global/dynamical phase-space
property, not generically identical to a local speed-collapse surface.

Separate these notions:

```text
freezing surface: local admissible-speed degeneration;
turning point: radial component of one trajectory vanishes;
separatrix: phase-space boundary/double-root structure in a specified physical domain.
```

### Ergosphere and selector replacement

The manuscript correctly notes that `W=partial_eta` becomes spacelike inside the Kerr
ergoregion and the compact-oval argument then fails. It next says that analytic
continuation is "equivalently" a re-anchoring to Doran/ZAMO and leaves the closed forms
unchanged.

Changing `W` changes the rail constraint `-u.W=Ehat` and therefore changes the
controlled physical problem. Equality of analytically continued expressions is not
by itself equivalence of the two optimization problems.

Use one of these formulations:

- treat the formulas inside the ergoregion as mathematical analytic continuations,
  without claiming they are physical optima for the original selector; or
- define a timelike selector globally, rederive its indicatrix/Hamiltonian, and prove
  the relation to the exterior problem.

### Acceptance criteria

- [ ] The control set is explicitly the boundary or the filled ellipse.
- [ ] Any convexification and absence of relaxation gap are proved or not claimed.
- [ ] Normality/nontriviality assumptions are stated.
- [ ] PMP is described as necessary unless a separate existence/minimality proof is
      supplied.
- [ ] The terminal manifold and free-time protocol are unambiguous.
- [ ] Freezing, turning, and separatrix are not identified with one another.
- [ ] Analytic continuation and physical selector replacement are distinguished.

---

# Part IV — Special-function and verification claims

## 7. MAJOR — The genus-two Kronecker--Eisenstein terminology is not established by the present derivation

### Current result that is supported

The repository supports the following conservative statement:

- the chosen frozen-source differential is reduced exactly to Abelian integrals on a
  genus-two hyperelliptic curve plus depth-two iterated Abelian integrals;
- for tested curves, period matrices and theta series give stable numerical
  evaluations;
- truncated lattice/nome sums converge rapidly in the tested examples;
- a five-dimensional numerical span appears stable to high precision in the selected
  data.

This is already a substantial result.

### Claims that exceed the evidence

Around lines 62--66, 1070--1090, and 2360--2395 the paper describes the weight-two
objects as a universal genus-two Kronecker--Eisenstein dilogarithm, asserts geometric
convergence, sometimes calls the basis irreducible, and states that every symbolic
identity was independently confirmed in Mathematica.

Problems:

1. Brown--Levin (`https://arxiv.org/abs/1110.6917`) is an elliptic/genus-one
   construction. It does not by itself identify an arbitrary genus-two theta-log
   iterated integral with the same named basis.
2. Current higher-genus constructions use additional structures and hypotheses; see
   `https://arxiv.org/abs/2306.08644` and `https://arxiv.org/abs/2406.10051`.
3. The repository scripts numerically truncate genus-two theta lattices and compare
   values. They do not currently prove a general convergence theorem or an identity
   with a canonical higher-genus Kronecker--Eisenstein polylogarithm.
4. Numerical rank five, even to very high precision, is strong evidence rather than a
   proof of exact dimension or motivic independence.
5. `VERIFICATION_STATUS.md`, lines 74--77, explicitly says the genus-two theta pieces
   were checked through Sage/`abelfunctions`, not Mathematica. This contradicts the
   abstract's "every symbolic identity" statement.
6. Around line 2385, `U_k` are called "third-kind elliptic integrals" in an
   off-separatrix genus-two discussion. They are hyperelliptic Abelian integrals; their
   first/second/third kind depends on the divisor/pole structure of each differential.

### Required correction

Unless a precise identification theorem is added, use terminology such as:

> depth-two iterated Abelian integrals on the genus-two frozen spectral curve, with a
> rapidly convergent numerical theta/nome representation for the tested period
> matrices.

If the Kronecker--Eisenstein name is retained, provide:

- the precise higher-genus definition being used;
- a map from each paper kernel/letter to that definition;
- divisor and monodromy data;
- convergence hypotheses and proof that the computed curves satisfy them;
- normalization conventions;
- a theorem identifying the paper's `W_kj` with the named objects.

The five-dimensional compression should remain explicitly "numerically observed" or
"strong numerical evidence" unless an independence/spanning proof is supplied.

### Acceptance criteria

- [ ] Genus-one references are not used as direct proofs of genus-two identities.
- [ ] The named higher-genus function class is defined precisely or the terminology is
      downgraded.
- [ ] Numerical convergence is distinguished from a general convergence theorem.
- [ ] Rank five and irreducibility are labelled conjectural/numerical unless proved.
- [ ] Mathematica and Sage verification claims match `VERIFICATION_STATUS.md`.
- [ ] `U_k` are classified using correct hyperelliptic differential terminology.

---

## 8. MAJOR/MINOR — Abstract and conclusion aggregate conditional results as theorems

The abstract currently includes all of the following as established main results:

- Kodama energy conservation;
- absence of physical evaporative inversion;
- an existence theorem for conformal plunge inversion;
- rigorous first-order adiabatic corrections;
- a universal genus-two KE dilogarithm;
- a separatrix-tracking theorem with physical moving-root interpretation;
- independent Mathematica confirmation of every symbolic identity.

These claims have different evidentiary status and should not be bundled together.

### Suggested status-correct replacements

Until the blocking derivations are complete, use wording along these lines:

```text
The rail actively maintains the Kodama/conformal-selector energy; this is not a
Noether conservation law for free fall.
```

```text
For slowly varying backgrounds we derive an exact special-function reduction of the
frozen-family source. Its identification with the first variation of the full
non-autonomous boundary-value problem remains conditional on the variational
derivation described in Section ... .
```

```text
The fixed-endpoint no-inversion follows from a proved algebraic comparison together
with a monotonicity lemma proved in a sub-regime and numerically supported elsewhere.
```

```text
The resulting depth-two genus-two iterated Abelian integrals admit rapidly convergent
theta-series evaluations in the tested examples.
```

```text
Algebraic coefficient identities and separatrix residues were cross-checked in
Mathematica; genus-two theta evaluations were checked with Sage/abelfunctions.
```

### Other overbroad statements to revise

- The geometric dictionary around lines 240--245 uses `if and only if`/uniqueness-like
  language connecting indicatrix shapes to spacetime classes. The three studied
  examples do not establish a global classification theorem. Restrict the statement
  to "in the models studied here."
- "Well-posed" should be replaced by "formulated as a normal PMP extremal problem"
  until existence/global minimality is proved.
- "Kodama energy conservation" in the abstract should say "active rail maintenance
  of the Kodama energy."
- "Airtight" around the unified theorem must be removed until Issue 1 is closed.

---

# Part V — Verified material that should be preserved

## 9. VERIFIED / PRESERVE — Foundational corrections already made successfully

The following earlier corrections are sound and should not be undone:

### Active rail invariant versus Noether conservation

The paper now correctly distinguishes an actively maintained quantity
`-u.W=Ehat` from a conserved Noether charge of free fall. The exact identity

```text
d(-u.W)/d tau = -a.W - u^a u^b nabla_(a W_b)
```

around lines 187--203 is correct. It also correctly distinguishes `a.u=0` (no kinetic
work on the particle/mass-shell preservation) from generally nonzero `a.W` (power
against the selector).

### Non-autonomous transversality

The statement that free-arrival transversality gives `H(s_f)=0`, not `H=0` throughout
the non-autonomous interior, is correct. The extended-Hamiltonian statement should be
retained.

### Vakonomic versus d'Alembert

The earlier false claim that convexity makes vakonomic and d'Alembert equations
coincide has been removed. Keep the current distinction. The remaining problem is
only the separate convex-set/existence wording discussed in Issue 6.

### Incoming versus outgoing Vaidya convention

The statement that physical radiation/emission uses outgoing retarded Vaidya is
correct. The required change is to derive and test that physical branch before making
universal conclusions.

### Doran azimuthal sign

The consistent relation

```text
phi_D = phi_BL - phi_shift
```

was fixed using horizon-log residue cancellation and should be preserved.

### Third-kind residue balance

The divisor accounting for `dr/(r-2m)`—residue `+2` at the horizon point and `-1` at
each point at infinity—is consistent with the residue theorem and should remain.

---

## 10. VERIFIED / PRESERVE — Algebraic and numerical reductions

The following checks passed at the audited commit or immediately preceding corrected
state:

### Generic reductions

Command:

```bash
python3 reproduce_reductions.py
```

Observed result:

- Thakurta--Kerr `tau` energy-source reduction: exact zero polynomial residual;
- Thakurta--Kerr `tau` angular-source reduction: exact zero polynomial residual;
- Thakurta--Kerr `t` energy/angular reductions: exact zero polynomial residuals;
- Vaidya mass-source reduction: exact zero polynomial residual;
- on-curve `t`-clock identity: exact.

These reductions should be reused after the dynamical source is corrected.

### Separatrix classification

Command:

```bash
python3 separatrix_classification.py
```

Observed result:

- correctly distinguishes physical exterior double roots from negative-radius
  algebraic degenerations;
- confirms that the Vaidya/Schwarzschild `7.0266` example is algebraic only;
- identifies the physical external Kerr `t`-branch root near `J=-8.0535`.

Treat this script as the authoritative regression test for terminology.

### Fixed-endpoint reductions

Command:

```bash
python3 no_inversion_reduction.py
```

Observed result:

- Lemma A symbolic factorization/proof passes;
- Lemma B partial/grazing checks pass;
- numerical evidence supports the remaining monotonicity domain but is not a proof.

Additional current progress is recorded in:

- `TODO_lemmaB_closure.md`;
- `no_inversion_schwarzschild_frozen.py`;
- `no_inversion_schwarzschild_closedform.py` if present in the working version.

### Vaidya algebraic decomposition

Command:

```bash
python3 VaidyaMetric/vaidya_asymmetry.py
```

Observed result:

- the quadrature identities `accr=A+B` and coefficient-level `evap=A-B` pass;
- the physical interpretation of their difference/sum has the sign-bookkeeping error
  described in Issue 4.

### Important limitation

Passing any of these algebraic tests does not validate the missing non-autonomous
variational step. Keep the regression tests, but do not cite them as proof of Issue 1.

---

# Part VI — Reproducibility, build, and editorial issues

## 11. Build status

The main CQG-style manuscript compiles successfully after full LaTeX passes:

```bash
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Observed output: approximately 61 pages, with no unresolved references or citations.

The PRD/REVTeX mirror also compiles, approximately 29 pages, but generates numerous
float-placement warnings.

### Layout defects observed

The main manuscript contains overfull boxes including approximately:

- `166 pt` near the fixed-endpoint proof/table region;
- `106 pt` and `68 pt` in later tables;
- `23 pt` in the appendix.

The PRD mirror has an overfull box around `93 pt` and several "float stuck" warnings.
These do not invalidate the mathematics but must be fixed before submission because
some equations/tables may be clipped or unreadable.

### Acceptance criteria

- [ ] Both manuscript variants compile from a clean auxiliary-file state.
- [ ] No undefined references/citations.
- [ ] No materially clipped overfull box.
- [ ] Figures and captions fit the target journal layout.
- [ ] The two manuscript variants contain the same scientific qualifications.

---

## 12. Reproducibility limitations

- There is no unified automated test suite; the repository relies on standalone
  scripts. Add a small regression runner that executes the fast exact checks and
  records expected tolerances.
- `ThakurtaMetric/conformal_source_Jterm.py` is extremely slow under the current
  SymPy environment. Its key coefficient identities are reproduced more quickly by
  `reproduce_reductions.py`; document which script is authoritative.
- Mathematica was installed in the audited environment but did not return even a
  trivial test expression during the audit. Existing `.wl` files could be inspected,
  but the full Mathematica suite was not rerun. Do not represent this audit as a new
  independent Mathematica execution.
- `VERIFICATION_STATUS.md` states that genus-two theta pieces require
  Sage/`abelfunctions`. Ensure the environment lockfile or installation instructions
  permit another researcher to reproduce those calculations.
- If numerical plots are used as evidence for global sign claims, save the parameter
  domain, grid, solver tolerance, event handling, and worst-case margin in a machine-
  readable output file.

Suggested fast regression entry point:

```text
1. reproduce_reductions.py
2. separatrix_classification.py
3. no_inversion_reduction.py
4. VaidyaMetric/vaidya_asymmetry.py
5. LaTeX build with warning summary
```

Do not include long exploratory scripts in the default fast suite.

---

## 13. Minor notation and terminology cleanup

### Spin symbol

The Thakurta--Kerr formulas have historically mixed `a` and `s`. The response file
says this was unified, but inspect the final manuscript and scripts again: all
definitions of `Delta`, `P`, wind terms, roots, and plotting labels must use the same
spin symbol.

### Energies and angular quantities

Maintain a notation table distinguishing:

```text
Ehat             actively maintained selector/rail energy;
E_eff            frozen stationary energy parameter;
P_phi            Pontryagin angular costate;
L_mech           mechanical spacetime angular momentum;
ell_eff          frozen Kerr shape parameter, after its relation is derived.
```

### Genus-two integral kinds

Do not call all `U_k` third-kind elliptic integrals. State the curve genus and classify
the differential from its poles/residues.

### "Teleological"

If retained, define it narrowly as endpoint-sensitive dependence arising from the
two-point/free-future boundary value problem. Avoid suggesting event-horizon
teleology unless that is actually what is being analysed.

### Geometry-wide dictionary

Replace global `iff`, "unique", or exhaustive ladder statements with "for the three
models studied" unless a classification theorem is supplied.

---

# Part VII — Prioritized revision plan

## Phase 0 — Freeze claims before doing new algebra

1. Mark the unified adiabatic result as conditional/frozen-source in abstract, body,
   captions, and conclusion.
2. Mark fixed-endpoint no-inversion as conditional outside the proved sub-regime.
3. Remove physical capture/escape language from negative-radius degenerations.
4. Correct the Vaidya `2A`/`2B` sign bookkeeping.

This prevents later work from building on claims already known to be overstated.

## Phase 1 — Rebuild the dynamical foundation

1. Define the exact non-autonomous optimal-control BVP.
2. Derive the extended Hamiltonian and all state/costate equations in one notation.
3. Establish the Legendre relation between `P_phi` and mechanical/frozen angular
   quantities.
4. Linearize the full BVP in the slow parameter.
5. Derive the observable correction including boundary and trajectory terms.
6. Validate it against symmetric finite differences of a re-shot non-autonomous BVP.

**Stop gate:** do not call the result an adiabatic theorem until this phase passes.

## Phase 2 — Reapply the exact algebraic machinery

1. Express the newly derived source on the frozen spectral curve.
2. Determine whether the existing `N_lambda/S^(3/2)` basis remains sufficient.
3. If additional Jacobi/Green-function terms occur, reduce them systematically.
4. Regenerate all `c_k`, `Q_kj`, boundary terms, and plots.
5. Check behaviour near turning points and separatrices with matched asymptotics rather
   than a bare divergent expansion.

## Phase 3 — Repair physical classification

1. Create the global threshold notation table.
2. Regenerate every separatrix/degeneration example through
   `separatrix_classification.py`.
3. Derive the outgoing Vaidya problem explicitly.
4. Re-run physical inversion/penetration comparisons using the same endpoint protocol.
5. State theorem, bounded numerical result, or conjecture according to the evidence.

## Phase 4 — Decide how to close Lemma B

Choose analytic proof, interval-certified compact-domain theorem, or conditional
theorem. Do not continue to call the full scattering statement proved while the
`(DAGGER)` gap remains.

## Phase 5 — Normalize the special-function claims

1. Define the higher-genus function class precisely or use conservative iterated-
   Abelian-integral terminology.
2. Separate symbolic algebra checks from numerical period/theta evaluations.
3. Preserve rank-five as numerical evidence unless an exact dimension proof is added.
4. Align abstract, `VERIFICATION_STATUS.md`, and appendices.

## Phase 6 — Final journal pass

1. synchronize CQG and PRD variants;
2. fix notation, overfull boxes, and floats;
3. run the regression suite;
4. perform a final claim-by-claim audit of abstract and conclusion.

---

# Part VIII — Global acceptance checklist

The paper can reasonably be called correct and rigorous when all applicable boxes are
closed:

## Core dynamics

- [ ] Full non-autonomous BVP is stated.
- [ ] Its first variation is derived.
- [ ] Endpoint/time/transversality variations are included.
- [ ] `P_phi` versus mechanical angular momentum is resolved.
- [ ] Full BVP finite differences validate the analytic first-order result.

## Physical interpretation

- [ ] Outgoing Vaidya is derived and tested directly.
- [ ] Physical and coefficient-level accretion/evaporation comparisons have correct
      signs.
- [ ] Negative-radius degenerations carry no capture/escape interpretation.
- [ ] Every threshold symbol has exactly one definition.

## Theorem status

- [ ] Lemma B is proved on the claimed domain, rigorously interval-certified on a
      stated compact domain, or consistently labelled conditional.
- [ ] PMP necessary conditions are not presented as global existence/minimality
      without a separate proof.
- [ ] Ergosphere continuation is not confused with an equivalent physical selector.

## Special functions

- [ ] Higher-genus naming is defined and cited accurately.
- [ ] Numerical rank/convergence is not described as exact proof.
- [ ] Mathematica/Sage claims match what was actually verified.

## Reproducibility and presentation

- [ ] Fast exact regression scripts pass.
- [ ] Non-autonomous BVP validation passes under mesh/tolerance refinement.
- [ ] Both LaTeX variants compile without unresolved references or clipped content.
- [ ] Abstract, theorem statements, captions, tables, appendices, and conclusion all
      report identical logical status.

---

# Part IX — Instructions for Claude or another revising agent

Use this workflow when applying corrections:

1. **Do not perform a broad stylistic rewrite first.** Resolve the two blockers and
   the notation/physical classification before polishing prose.
2. **Preserve verified algebra.** Treat `reproduce_reductions.py`, the Lemma A
   factorization, the Doran sign, and the residue balance as regression-protected.
3. **Do not claim that a quadrature identity validates a non-autonomous ODE/BVP.**
   Rename the current validation figures until a genuinely independent solver exists.
4. **Keep logical status explicit.** Use the labels `proved`, `conditional`,
   `computer-assisted on domain D`, `numerically observed`, and `conjectural`
   consistently.
5. **Before changing `J partial_J`, derive what `J` means.** Do not remove or retain the
   term based only on previous prose.
6. **For every separatrix statement, print and inspect `r_d`.** A discriminant zero is
   not automatically a physical separatrix.
7. **After each scientific change, run the relevant exact regression script.**
8. **Synchronize both paper variants.** A correction existing only in `main.tex` is not
   complete.
9. **Report unresolved gaps rather than hiding them.** A narrower correct theorem is
   preferable to a universal statement supported only by scans.

Suggested order of files to inspect:

```text
paper/main.tex
paper/main_prd_revtex.tex
ThakurtaMetric/rail_conservation.py
ThakurtaMetric/fig_phi_closed_validation.py
reproduce_reductions.py
separatrix_classification.py
no_inversion_reduction.py
TODO_lemmaB_closure.md
VaidyaMetric/inversione_fisica.py
VaidyaMetric/vaidya_asymmetry.py
VERIFICATION_STATUS.md
REFEREE_RESPONSE.md
```

When a point is fixed, update this report or create a response table with four
columns:

```text
Issue | Change made | Proof/test | Remaining limitation
```

An issue is not closed by changing prose alone if it asks for a derivation, and it is
not closed by adding a derivation alone if the abstract, captions, scripts, and mirror
paper continue to make the old claim.

---

## Final assessment

The manuscript contains publishable-looking components: the controlled-invariant
identity, the frozen algebraic reductions, the separatrix classification machinery,
the Lemma A proof, and the special-function assembly are substantial. The main
scientific risk is that exact manipulation of the frozen source is currently being
mistaken for derivation of the non-autonomous optimum.

The correct strategy is therefore not to discard the genus-two machinery. It is to
derive the actual variational source first and then place that source into the already
developed reduction framework. If that step succeeds and the physical/conditional
claims are normalized as described above, the project will be much closer to a
rigorous paper.
