# Response to the Technical Referee Report (main8.pdf)

We thank the referee for a careful, competent report. We accept the central diagnosis:
the independent spot-checks pass, so the *calculations* are correct; the problem is a
mismatch between the strength of the claims and the level of proof supplied. Our revision
therefore (i) downgrades every overclaim to what is actually proven, (ii) states the domain
of each theorem, (iii) incorporates new closed-form and machine-checked results that were
completed after submission, and (iv) reorganizes the honesty of abstract/tables/conclusion.
Below, "DONE" marks a change already in the revised source; "PLANNED" marks a change
requiring a new derivation or a reproducibility artifact, scheduled for the next round.

## Spot checks (Sec. 3)
All PASS items are retained unchanged. The two "NOT ESTABLISHED" items (fixed-endpoint
no-inversion, ergosphere-as-PMP-optimum) and the "FAIL as truncation test" (Fig. 9) are
addressed in 4.9, 4.4, and 4.7 respectively.

## Major comments

**4.1 Endpoint problem / "between two events" — ACCEPTED. (PLANNED)**
Correct and important. In a non-stationary spacetime, curves sharing a spatial endpoint at
different $t,\eta,v$ are not the same event. We will add, per branch, an explicit statement
of: state variables; evolution parameter; admissible control set; initial event/hypersurface;
terminal manifold; which terminal coordinates are fixed; terminal cost; whether the final
selector-time is free; and the transversality conditions. "Fixed endpoints" will be replaced
by "fixed spatial endpoints" unless both time coordinates are genuinely fixed, and the $t$-
vs $\tau$-comparison will state whether the optimized curves terminate at the same event or
only the same spatial location.

**4.2 Convexity/existence in Sec. 2.1 — ACCEPTED. (DONE)**
The indicatrix is the compact strictly convex *closed curve* bounding a convex region; the
control is its direction $\theta\in S^1$. Rewritten: uniqueness of the support-function
maximizer follows from strict convexity; no Filippov relaxation is needed because $S^1$ is
compact and the maximizer single-valued — we no longer claim the velocity locus is a convex
subset of the plane, nor convexify it. Added exclusion of $p=0$ (nontriviality) and explicit
deferral of freezing/horizon degeneracies.

**4.3 Force/work interpretation — ACCEPTED. (DONE)**
Removed the "magnetic-type" framing. Now: $a\cdot u=0$ is *automatic* for any rest-mass-
preserving four-force (derivative of the shell), not special, and does not imply zero work
w.r.t. a physical observer. We distinguish (i) rest-mass preservation, (ii) observer-measured
power, (iii) the constraint force, (iv) non-stationarity-compensation power. The thrust/fuel
paragraph will be qualified with the observer/exhaust assumptions or removed. (PARTLY DONE.)

**4.4 Optimal-control problem inside the ergosphere — ACCEPTED. (PLANNED)**
We agree the Doran coordinate change and the re-anchoring of the selector to a ZAMO/Doran
congruence are inequivalent operations, and that closed forms cannot be asserted invariant
without derivation. Revision route chosen: **restrict all optimal-control theorems to the
timelike-$W$ region and label the ergosphere continuation as purely analytic**, clearly
separated from proved control results (the alternative — a globally timelike selector with
re-derived Hamiltonians and matching conditions — is noted as future work).

**4.5 "Horizon"/"ergosphere" in conformal Kerr — ACCEPTED. (PLANNED)**
Apparent/trapping horizons are not conformally invariant. We will distinguish seed null
surface, conformal Killing horizon, apparent/trapping horizon, and event horizon; compute or
cite the relevant null expansions; and rename $r=2M$ the conformal-stationary-limit surface
unless a stronger interpretation is demonstrated. "Outside the horizon"/"horizon-penetrating"
will be revised accordingly.

**4.6 Adiabatic correction not from the full non-autonomous $H$ — ACCEPTED (deep). (PLANNED)**
The frozen-shell assumption is not justified by $\hat E$-conservation alone; an off-shell
$O(\dot\lambda)$ term from costate transport ($p_\eta,p_v$) can arise because $H\ne0$ in the
interior. We will re-derive the first-order correction by a multiple-scale expansion of the
full extended system $H_{\rm ext}=p_s+H=0$, including endpoint conditions and costate
transport, and prove that the parameter-derivative formula is the complete first-order term
(or add the missing off-shell term), validating against the full non-autonomous BVP rather
than against re-integration of the same frozen integrand. This is directly connected to 4.7.

**4.7 Fig. 9 scaling contradicts first-order accuracy — ACCEPTED. (DONE + PLANNED)**
Confirmed: residual$/\varepsilon = 2.8\times10^{-4}$ is constant (linear in $\varepsilon$),
not $O(\varepsilon^2)$. Caption corrected (DONE): the linear growth measures the size of the
retained correction plus a numerical floor and is *not* an $O(\varepsilon^2)$ truncation test.
PLANNED: a log--log convergence study against the full non-autonomous flow (expected slope
two if the parameter-derivative formula is complete), separating discretization, quadrature,
and truncation errors.

**4.8 Turning-point claims mutually inconsistent — ACCEPTED. (DONE + PLANNED)**
We removed the claim of agreement "through the turning point"; Fig. 9 panels are now stated on
compact sub-arcs bounded away from the turning point (DONE). PLANNED: either display and
validate the matched local solution near the simple turning point (showing the divergent
$S^{-3/2}$ pieces cancel once the moving-root shift is included) and compare to the full ODE,
or restrict all plots to compact sub-arcs; the separatrix case is flagged for separate scaling
(divergent radial period).

**4.9 Fixed-endpoint no-inversion overclaimed — ACCEPTED. (DONE)**
This is the central correction. We downgraded every "never inverts" to the proven/tested
regimes (DONE, abstract-level captions and the theorem statement). Concretely, the revised
Lemma (B) now records the new results obtained after submission:
- an exact closed form for $\Phi_\tau'$ (Abel integration by parts), which proves
  monotonicity in the grazing ($V_0\le2V_{\min}$) and quarter ($V_0\le4V_{\min}$) regimes,
  hence for all $r_0\le R_\ast(E)$;
- in the static limit, the large-$r_0$ asymptotics
  $\Phi_\tau=\arctan\sqrt{V_0/V_{\min}-1}-b^{-3/2}\sqrt{1/V_{\min}-1/V_0}+\cdots$, a single
  *non-degenerate* maximum with $\Phi_\tau''(V_{\rm pk})=-\tfrac12 b^{9/4}V_0^{-5/4}<0$,
  proving (B) for all sufficiently large $r_0$;
- a validated interval-arithmetic (CAP) certificate of the single-crossing property,
  completed for $r_0=10M$ over the whole continuum $r_{\min}\in(2.6,6.05)M$.
The universal statement is now stated as a theorem *in these regimes*, machine-certified at
representative $r_0$, and a numerically supported conjecture in full generality. Every "never"
was audited. We agree the pointwise $n_t/n_\tau>1$ bound is insufficient and say so.

**4.10 $J$ is a PMP costate before a physical momentum — ACCEPTED. (PLANNED)**
We will derive, branch by branch, the relation between the PMP costate $p_\phi$, the
Finsler/Jacobi momentum, and $g(u,\partial_\phi)$, and state where they coincide (stationary
shell / up to a position-dependent factor) versus differ (accelerated constrained worldline,
different Legendre maps for $t,\tau,\eta$). Notation and conformal scaling of $J_{\rm eff}$
corrected accordingly.

**4.11 Evaporation sign flip vs outgoing model — ACCEPTED. (DONE + PLANNED)**
Caption corrected to label the negative-$m$ case as an *ingoing model* statement, not
universal evaporation (DONE). PLANNED: cleanly separate (i) ingoing Vaidya either sign as a
math model, (ii) ingoing accretion, (iii) outgoing decreasing-mass evaporation, and recompute
any penetration/timing claim in the outgoing metric with a positive-mass domain.

**4.12 Higher-polylogarithm framework — ACCEPTED. (PLANNED)**
We will separate three levels: exact algebraic reduction to length-two iterated Abelian
integrals (rigorous); a chosen numerical theta-series representation; and conjectural
statements about irreducibility / minimal basis dimension (the rank-five observation is
explicitly numerical and curve-specific, not a theorem). Absent a full higher-genus framework
(base points, paths, tangential regularization, monodromy, single-valued completion), we adopt
the neutral term "length-two iterated Abelian integrals" and drop "genus-two dilogarithm" as a
theorem-level claim. The "collapse to Bloch--Wigner" needs an explicit elliptic-polylog
formula and is downgraded.

**4.13 Separatrix / threshold / algebraic degeneration notation — ACCEPTED. (PLANNED)**
Distinct symbols throughout: $J_{\rm pen}(v_0)$ (dynamical penetration threshold), $J_{\rm sep}$
(physical capture/scatter separatrix with divergent period), $J_{\rm deg}$ (algebraic double
root at $r<0$). "Separatrix" reserved for the physical phase-space boundary.

**4.14 Auditable reproducibility — ACCEPTED. (PLANNED)**
We will archive a versioned release tied to the manuscript with exact symbolic expressions /
machine-readable coefficient files, the Mathematica sources used as *independent* checks
(commit hash, environment, precision), and automated tests that compare genuinely independent
formulations. (Post-submission we have added an independent Mathematica cross-check of the
frozen no-inversion derivations, `verify_lemmaB_mathematica.wls`.) For BVP claims we will
report solver tolerances, mesh refinement, branch-cut sensitivity, and working precision.

**4.15 Narrow or split — ACCEPTED in spirit. (PLANNED)**
We agree the scope is too broad for uniform rigor. Preferred split: (1) controlled-rail
formulation + FLRW/Vaidya with a rigorous PMP theorem and outgoing-Vaidya comparison; (2)
conformal Kerr / ergosphere outside the timelike-selector boundary and protocol dependence of
inversion; (3) the higher-genus adiabatic reduction as a separate mathematical paper after the
iterated-integral framework is formalized. If retained as one article, the abstract/conclusion
are shortened and every theorem carries explicit hypotheses.

## Minor / editorial (Sec. 5) — ACCEPTED
Signature/units and dimensions of $\hat E,J,A,m$; derivative notation $\dot m,A_{,\eta}$ instead
of $m0,A0$; "strictly convex closed curve" not "convex set" (DONE via 4.2); Eq. (3) derived
from the terminal manifold; replace "teleological" by "endpoint-sensitive costate memory";
drop "magnetic-type" (DONE); clarify "proper-time brachistochrone"; horizon coefficient
degeneration is chart-dependent; distinct symbols for FLRW scale factor vs Kerr spin ($a$ vs
$A$); "freezing surface" as a kinematic boundary; dimensionless log arguments; base
point/sheet/path for every Abelian integral; inverse-Weierstrass branch; Doran sign stated
once; fix "App. Appendix A.1" duplicated labels; Fig. 9 caption (DONE); number the Eq. (40)
threshold; enlarge multi-panel figures; exact expressions in Tables A2/A3; residual norm/interval/
precision stated; drop "geodesic costs no fuel"; AI-assistance statement per journal policy;
ref. [31] DOI; abstract to at most three proven results + one limitation; conclusion to
separate proved / asymptotic / numerical / conjectural in distinct sentences.

## Summary of changes already in the revised source
Downgraded no-inversion theorem statement and captions with the new closed-form + asymptotic +
CAP results (4.9); corrected convexity/existence (4.2); corrected force/work interpretation
(4.3); corrected Fig. 9 caption and turning-point wording (4.7, 4.8); qualified the ingoing
evaporation caption (4.11); intro rail-force wording (4.3). The paper recompiles (62 pp).
All remaining items are scoped above as PLANNED with a concrete route.
