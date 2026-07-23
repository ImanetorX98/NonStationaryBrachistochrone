# Response to the technical review (17 objections)

Point-by-point record of how each objection in `Obiezioni_tecniche_main6.pdf` was
addressed. All 17 were evaluated; every one was found (at least partly) founded and
acted upon. Papers: `paper/main.tex`, `paper/main_prd_revtex.tex`.

Legend: **[C]** critical, **[M]** major, **[m]** minor (referee's grading).

---

## Results-affecting corrections

### #5 [C] — Conformal adiabatic source incomplete (missing `d_J F`)
**Founded (verified numerically).** Under `g^TK=A^2 g^Kerr`, both Kerr charges carry
weight `A^-1`: `E_eff=Ehat/A` AND `J_eff=J/A`. Holding `J` fixed drifts the conserved
`J_TK=g^TK(u,d_phi)` by ~5.5 (numerics). The first-order source is the Euler operator
`(E d_E + J d_J)F`, not `E d_E F` alone. The `d_J F` piece reduces in the identical
second-kind basis, `c_k^tot=c_k^E+c_k^J` (no new special function); verified by full
block assembly (Weierstrass + genus-2 dilog) to `1e-8`. Vaidya (mass parameter, not
conformal) unaffected.
Commits `9e808ad`, `ff5dc74`. New: `ThakurtaMetric/conformal_source_Jterm.py`;
`reproduce_reductions.py` now outputs `c_k^J`.

### #7/#8 [C] — "Separatrix" is a negative-radius algebraic degeneration; `Q_2` not linear in `J`
**Founded (verified).** `separatrix_classification.py`: `Q_2` is quadratic in `J`
(terms `J^2`, `(J-a)^2`). The discriminant `Res_r(Q_2,d_r Q_2)=0` has six real roots.
The physical separatrix (external double root `r_d>r_+`, two real turning points
merging) is the retrograde `J_c^-=-8.054`, `r_d=+3.514>2M`; `J=+19.089,-18.671` sit at
`r_d=-6.62,-6.59<0` — algebraic genus-degenerations `J_deg`, NOT physical separatrices.
App. B.1 rewritten: physical tracking `dJc^-/dE=+2.44` (the earlier `-/+148-150` belonged
to the `J_deg` roots), physical `b_i` via the Euler source. `J_deg` distinguished from
the dynamical penetration threshold `J_c(v_0)`. Schwarzschild `tau` at `E=1.4` has no
external double root: its `r_d=-3.36` is `J_deg`.
Commit `bb7561c`. New: `separatrix_classification.py`.

### #11 [M] — Fixed-endpoint no-inversion is not a theorem
**Founded.** Reduced to two lemmas (parametrize each branch by turning radius `r_min`;
`Phi_br(r_min)=int_{r_min}^{r0} F_br(r;J_br(r_min))dr`):
- **(A)** `F_t>F_tau` at matched `r_min` => `Phi_t>Phi_tau`. **PROVEN in closed form**
  (algebraic): via the factorization `N_G=(r-r_min)[4r Delta/(r_min-2M)] W(r)` with `W`
  LINEAR in `r`, `W(r_min)>0` and slope `w1>0` following from `J_t<J*,J**` (downward
  parabola `Q_2`, negative vertex, `Q_2(r_min,J*),Q_2(r_min,J**)<0` manifestly).
- **(B)** `Phi_tau` decreasing. **Conditional/partial**: proven in closed form for the
  grazing sub-regime `V(r0)<=2V(r_min)` (and extended to `~4V(r_min)` via convexity of
  `W`, `W''>0` proven); the exact threshold is the transcendental deflection peak
  `r_pk~3M`, beyond elementary certificates.
The fixed-endpoint turning points lie in the scattering regime `r_min>r_pk`, so the
no-inversion holds there; it is a theorem modulo (B). At FIXED `J` the order is the
opposite (t deeper) — the fixed-endpoint result is a reversal via `J -> Phi`, confirming
the pointwise `n_t/n_tau=E/f>1` bound is insufficient.
Commits `9e422f4`,`163317f`,`882ca0b`,`0aeb22a`,`6afc7eb`. New:
`no_inversion_reduction.py`; open step in `TODO_lemmaB_closure.md`.

---

## Foundational rigor (control theory)

### #1 [C] — Constrained-energy derivative: missing four-force term
**Founded.** `eq:cost` is now the exact identity
`d(-u.W)/dtau = -a.W - u^a u^b nabla_(a W_b)`. The rail force is magnetic-type in the
kinetic sense `a.u=0` (no work on the particle), distinct from `a.W` (power against the
selector); the rail supplies `a.W=-uu nablaW` to hold `Ehat` fixed. Commit `626b161`.

### #2 [C] — Transversality `H=0` vs `dH/ds != 0`
**Founded.** `H=0` stated at the free endpoint `H|_{s_f}=0`; the extended Hamiltonian
`H_ext=p_s+H=0` is the global invariant; `H` drifts in the interior (consistent with
`p_v(r_1)=0` but `p_v!=0` on the path). Commit `626b161`.

### #3 [M] — Convexity does not imply vakonomic=d'Alembert
**Founded.** False claim removed. Strict convexity => unique maximizer (smooth control,
no Filippov relaxation, the set is already compact convex); the problem is time-optimal
control, so the Pontryagin (vakonomic-type) reduction is correct and d'Alembert is not
invoked. Commit `9e3b6bf`.

### #4 [C] — Domain of the elliptic indicatrix
**Founded.** New paragraph: the oval is compact strictly convex only where `W` is
timelike and nondegenerate. Freezing surface (`R->0`), horizon (radial coeff degenerate)
and Kerr ergosphere (`d_eta`, hence `W`, spacelike) treated separately; ergosphere via
analytic continuation or a timelike Doran/ZAMO re-anchoring. Commit `626b161`.

### #13 [M] — Third-kind differential residue balance
**Founded (verified).** `dr/(r-2m)` has residue `+2` at the horizon Weierstrass point
and `-1` at each of the two points at infinity (sum zero, residue theorem); the Fay
theta-ratio carries the full divisor `2(r=2m) - inf_+ - inf_-`. Commit `626b161`.

---

## Physical / modeling

### #9 [M] — Ingoing vs outgoing Vaidya
**Founded (partly pre-existing).** The depth ordering `r_min^t<r_min^tau` is
branch-intrinsic (`Ehat/f>1`), independent of `sign(m')` and of the flux direction.
Physical evaporation is the OUTGOING metric `ds^2=-f du^2-2 du dr` (retarded `u`),
derived and compared in the asymmetry paragraph (now cross-referenced); the `m<=0`
linear-model "inversion" is an unphysical ingoing-model extrapolation. Commit `6d4335c`.

### #10 [M] — "Breathing-indicatrix theorem" incoherent
**Founded.** Restated precisely: in each geometry BOTH the centre and the semi-axes
breathe (`d_m c_r=-2/r`, `d_m R_r`, `d_m R_phi` all nonzero; Thakurta via
`d_A E_eff=-Ehat/A^2`). The geometries are distinguished by the wind ORIENTATION
(radial centre-wind vs azimuthal frame-drag) and the CLOCK singularity location
(horizon dlog vs infinity dipole), not by which datum varies. Commit `cdda868`.

### #12 [M] — Existence of the conformal inversion: missing hypotheses
**Founded.** Made conditional: `cond(A_freeze)>0` requires `Ehat>1`, `r>2M`, `Delta>0`
(exterior scattering); for `Ehat<=1` no root is guaranteed. Commit `cdda868`.

---

## Honesty / claims / notation

### #14 [M] — Special-function claims too strong
**Founded.** "Dimension exactly five" -> strong numerical evidence (stable rank 5, sv to
`1e-79`), not a proof; a motivic/Hodge count or explicit independent basis would be
needed. Irreducibility to weight one -> conjectural (theta-divisor obstruction). Commit
`cdda868`.

### #15 [M] — Inconsistent Doran sign
**Founded (fixed by computation).** Horizon log residues of `phi_BL` and `phi_shift` are
EQUAL (`1.0324=1.0324`), so they cancel under SUBTRACTION: the single consistent sign is
`phi_D=phi_BL-phi_shift` (`eq:doran-phi`). The two stray `+phi_shift` occurrences fixed in
each paper. Commit `d823792`.

### #16 [m] — Terminology
**Founded.** Kodama-energy caption: "held constant (actively rail-controlled, not a
Noether constant)". "teleological" defined once as "endpoint-sensitive, from the
free-future two-point problem---not causal". Commit `6d4335c`.

### #17 [m] — Minor mathematical form
**Founded.** degree-six -> degree-five numerator (`A=sum_0^5`); spin symbol unified to
`a` (was mixed `s`/`a`); reduction written as `d/dr(A/R)+sum c_k r^k/R`; `f2` is defined;
`E,Ehat,E_eff` and `J_deg` vs `J_c(v_0)` distinguished. Commits `ff5dc74`, `6d4335c`.

---

## Section 3 of the review (verified derivations, NOT to rewrite)

Confirmed correct and preserved: elliptic indicatrix + branch Hamiltonians; static
turning laws; FLRW degeneration; cusp `2/3` law; Weierstrass uniformization; the 11x11
reduction identity; the triple-pole cancellation (now noted to be an algebraic identity
at any double root, hence valid on the physical branch too).

## Outstanding

- **Lemma B full closure** (`#11`) beyond the grazing sub-regime — a genuine
  scattering-angle inequality tight at the transcendental `r_pk`. Routes and progress in
  `TODO_lemmaB_closure.md`. Everything else is closed or rigorously downgraded.
