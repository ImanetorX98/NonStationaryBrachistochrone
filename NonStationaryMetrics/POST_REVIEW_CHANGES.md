# What changed in Papers I and II after the CQG review

Manuscript **CQG-116884** (Paper II, conformal Kerr). Two reports:
**R1** = the technical report, 12 major + 10 minor comments; **R2** = the
bibliography/self-citation report.

This file is the durable record. `paper2/REVISION_LOG.md` holds the working
detail; the Author Response (`paper2/response_to_referees_CQG.tex`) is what CQG
receives. Where those three disagree, **this file and the response are the ones
to trust**, because they were written last.

---

## 1. Errors the review caught in our work

Five, all real, all fixed.

| # | What | Where it was | Comment |
|---|---|---|---|
| 1 | $J=-J_c$ classified as a cusp in the body and a corner in the appendix | §3.2, Fig. 9, Table A1, App. A.1 | R1 major 4 |
| 2 | "a PMP extremal is guaranteed only *locally* optimal" | **both** `paper1.tex` and `paper1_JMP.tex` | R1 major 2 |
| 3 | `cap_full.py` imported by three CAP runners but never archived | repository | R1 major 10 |
| 4 | printed SHA-256 hashes did not match the release | Table A3 | R1 major 10 |
| 5 | figure and table quoted different adiabatic slopes | Fig. 3 vs Table A3 | R1 major 7 |

**On (1): both printed statements were wrong.** The rederivation from the
Hamilton equations gives a four-case law, not a trichotomy — see §3 below.

**On (4) it was worse than the referee could see.** The same two files carried
*three* different SHA-256 prefixes: one printed in the manuscript, one in the
release the referee cloned, one in the working tree.

## 2. Errors we found ourselves while chasing the review

Four more, none visible from outside.

1. **Every figure generator failed on a fresh clone.** All thirteen inserted the
   repository root on `sys.path` while `paper_style.py` lives one directory
   down. Reproduced, then fixed in all thirteen.
2. **A silent numerical bug.** `r_turn()` returned the grid point just *inside*
   the forbidden region, so the $|J|>J_c$ curve was dropped from a published
   log–log panel with no error raised.
3. **A wrong page number** in `GiannoniPiccioneVerderesi1997`: we had 6178, the
   correct range is 6367–6381.
4. **A wrong count in a caption**: "six $\varepsilon$ values per curve" where
   there are five.

## 3. The mathematics that changed

### 3.1 Local classification at the conformal stationary limit (R1 major 4)

The shape radical $\mathcal S$ depends on $J$ only through $J^2$, so no argument
based on it can distinguish $\pm J_c$. The Hamilton equations can:
$\tilde p_\varphi(r_e)=J-J_c$ is **linear** in $J$. With the exact identities
$p_r^2=\mathcal S/(D_E^2\Delta^2)$ and
$r\Delta-J_c^2D_E=(r-2M)(r^2+J_c^2)$:

| case | root of $p_r^2$ | $\mathcal T(r_e)$ | $r_e$ attained | shape |
|---|---|---|---|---|
| $\|J\|>J_c$ | none | $\neq0$ | no | periapsis, exponent 2 |
| $\|J\|<J_c$ | simple | $\neq0$ | finite parameter | cusp, exponent 3/2 |
| $J=-J_c$ | double | $\neq0$ | **asymptotic only** | finite tangent, exponent 1 |
| $J=+J_c$ | double | simple | finite, crosses | regular crossing |

New closed forms: cusp amplitude $K=(EJ/a^2)\sqrt{2M/(a^2-E^2J^2)}$, which
**diverges** as $|J|\to J_c$ and so proves the two families do not connect;
retrograde decay rate $C=-\sqrt{4M^2+J_c^2}/(8M^2)$; prograde crossing rate.

The "reality wall" argument is **withdrawn**: at $|J|=J_c$ the wall falls for
both signs. The old counter-rotation story was wrong in detail too — near $r_e$
the retrograde marginal is locally *co*-rotating.

### 3.2 The control domain (R1 major 1)

New §2.1 with `lem:compact`, proved: the rail slice is nonempty and compact iff
$W$ is timelike and $\hat E\ge|W|$, a single point iff $\hat E=|W|$, noncompact
if $W$ is null or spacelike. Hence for $W=\partial_\eta$ the compact-control
problem **ends** at $r=2M$. Doran coordinates do not repair this.

### 3.3 Exterior retrograde separatrix as the flagship (R1 major 5)

The referee's numbers reproduce to 16 digits: $r_d=3.513905124011657M$,
$J_c^-=-8.053516003877019$. **Beyond what was asked**: both control-domain
inequalities hold strictly at $r_d$ ($g(W,W)=-0.430833$, $\bar v^2=0.700811$) and
on the whole approach arc, so it is a type-(E) object; and $r_d>2M$ persists
across $a\in[0.1,0.99]$, so the exterior character is structural.

### 3.4 The rail drift is the adiabatic parameter — **new**

Thakurta–Kerr is $g=A(\eta)^2\bar g$ with $\bar g$ the $\eta$-independent Kerr
seed, so $\mathcal L_Wg=2AA'\bar g=(2A'/A)g$: **$W=\partial_\eta$ is a conformal
Killing vector**, $\psi=A'/A$, $\nabla\!\cdot\!W=4\psi$. Hence

> $\dfrac{d\hat E}{d\tau}\Big|_{\rm geodesic}=\dfrac{A'}{A}=\varepsilon$,
> independently of the four-velocity.

The control effort and the expansion parameter are the same quantity. This is
the middle rung of the Killing → conformal-Killing → Kodama hierarchy made
quantitative, and the manuscript had never stated it.

A soliton selector would need $\mathrm{Ric}=(\lambda-\psi)g$, i.e. an Einstein
spacetime; TK is not, so the soliton class excludes this selector while the
conformal-Killing class contains it — which is why the hierarchy needs three
rungs.

### 3.5 Minimum thrust, and where it diverges — **new** (R1 major 9)

Holding the rail forces $g(a,W)=\varepsilon$. Cauchy–Schwarz in the rest space
of $u$ gives the **sharp** bound

> $|a|\ \ge\ \dfrac{\varepsilon}{\hat E\,\bar v}$, attained when $a\parallel W_\perp$.

Writing $A^2f=\hat E^2(1-\delta)$, the minimum thrust behaves as
$\varepsilon/(\hat E\sqrt\delta)$: it **diverges as an inverse square root at the
freezing surface** and is **finite at the conformal stationary limit**
($\to\varepsilon/\hat E$). The two boundaries of the control domain are of
different kinds — one a failure of compactness at finite cost, the other a
divergence of cost. The referee's prose caveat is now a theorem with a sharper
statement than the caveat had.

### 3.6 Curvature of the optical base — **new** (R2, and R1 major 2)

The optical projection is a **horizontally conformal submersion** with dilation
$\Lambda^2=\hat E^2/[f(\hat E^2-f)]$ — the Perlick factor itself. So the
Fuglede–Ishihara dichotomy (R2's theorem T2) is the abstract reason the optical
metric is a *conformal* rescaling of $h$. That is a real observation and §2.2
credits it.

Gauss curvature of the fixed-energy optical surface computed in closed form;
reduces as $\hat E\to\infty$ to the classical $-(2M/r^3)(1-3M/2r)$. Far field:

> $r^3K\ \longrightarrow\ -\dfrac{M(2\hat E^2-3)}{\hat E^2}$

so the distant optical geometry **focuses** for $\hat E^2<3/2$ and **defocuses**
for $\hat E^2>3/2$; critical rail energy $\hat E=\sqrt{3/2}$, asymptotic speed
$1/\sqrt3$. At the critical value $r^4K\to-23M^2/3$. Focusing produces conjugate
points, so this separates two regimes for distant extremals.

But $K\to0$ at large $r$ for every $\hat E$: **asymptotic flatness is
incompatible with the Myers hypothesis**, and positivity of the curvature
operator would force a space form. Those frameworks are therefore unavailable
here — now an argument with a computation behind it, not an omission.

### 3.7 Conjugate points and Maxwell points — **new** (R1 major 2)

On the frozen $\tau$-branch the optical geometry is Riemannian-pure, so the
classical second variation applies verbatim: $\ddot{\mathcal J}+K\mathcal J=0$.

> **$K(r)<0$ for every $r>2M$ if and only if $\hat E^2\ge3/2$.**
> Hence for $\hat E^2\ge3/2$ no exterior extremal carries a conjugate point, and
> every one is a **local minimiser** of arrival time.

Two facts fix it: the numerator of $K$ at the stationary limit is $-4\hat E^4M^4$,
negative for every $\hat E$; and its cubic coefficient is
$-(2\hat E^2-3)(\hat E^2-1)$. So the $\hat E^2=3/2$ threshold is **not merely
asymptotic** — it is global. Below it, positive curvature opens at large radius
whose inner edge runs to infinity as $\hat E^2\to(3/2)^-$: $r/M = 22.75,\,59.53,\,
576.7,\,5751.7$ at $\hat E^2=1.20,\,1.40,\,1.49,\,1.499$.

**But negative curvature does not give global minimality**, and not for a focusing
reason. The optical surface is $\{r>2M\}\times S^1$ — an **annulus**, not simply
connected. Cartan–Hadamard gives local minimality without uniqueness: the two
extremals exchanged by $\varphi\to-\varphi$ reconverge at swept angle $\pi$. So
for $\hat E^2\ge3/2$ the cut locus is **purely of Maxwell type**; global
minimality holds up to the first antipodal crossing and fails there by symmetry,
not by focusing. Below the threshold both mechanisms coexist.

Scope: established for the non-rotating frozen limit; the rotating case is not
claimed.

---

### 3.8 What the bibliographic audit produced — **new**

Reading the sources rather than their metadata found two errors of ours and four
results we then used. Ordered by what it changed:

**Bishop 1972 (Clairaut submersions).** We stated the theorem as "Clairaut iff
fibres totally umbilical". The theorem is: a Riemannian submersion **of complete
manifolds with connected fibres** is Clairaut iff the fibres are totally
umbilical **with mean curvature of gradient form**, $H=-\nabla\ln f$. Umbilicity
alone is not the condition. Corrected in the manuscript and the response, and
later confirmed from a primary source (Meena–Zawadzki's abstract).

**Giannoni–Masiello–Piccione 2001.** Described as "closest to the frozen-$A$
regime". Too generous: it concerns **null** geodesics and lensing, its machinery
rests systematically on the **Killing** property, and the conformally stationary
case enters as an extension under boundedness assumptions on the conformal
factor. Corrected. Reading it also gave a better citation for the
Cartan–Hadamard fact we rely on.

**Lecian, curvature-operator eigenvalues.** Produced the weighted-seed
correspondence of §3.4bis: with weight $e^{-2f^2}=A^{-2}$ the weighted metric
**is** the Kerr seed, $W$ is a genuine Killing field of it and
$\operatorname{div}_{\bar g}W=0$, matching their Eq. (32) verbatim. Their **T8**
(strictly parallel weight ⟺ Einstein) then draws the line: our weight
$f^2=\ln A$ is parallel iff $A$ is constant.

**Caponio–Corona–Giambò–Piccione 2024.** Their arrival time is differentiated
subject to the *semi-holonomic constraint defined by the Noether charge*. Our
rail is a semi-holonomic constraint of exactly that type — the structural reason
the theories look alike. Difference sharpened: there $K$ is a symmetry so the
constraint is preserved automatically; here $W$ is conformal Killing but not
Killing, and **a conformal Killing field conserves no charge along timelike
worldlines**.

**Caponio–Corona–Giambò–Piccione 2024, read in full.** Their treatment of the
semi-holonomic constraint is precise and marks the limit of what transfers.
Because $K$ is a symmetry, the Noether charge is automatically constant along any
critical point (their Lemma 3.3), so curves of constant charge are the only
candidates; that set is then a **$C^1$ closed submanifold** of the path space with
an explicit tangent space (Prop. 3.4), and critical points of the action are the
critical points of its restriction (Prop. 3.5). Their Lorentz–Finsler setting
also gives, for a *stationary Finsler spacetime*, the local product structure
$L_F=F^2(y,v)+\omega_y(v)\tau-\tfrac12\tau^2$ — the general form of which our
Randers reduction is an instance.

**And a risk this closed.** Their constraint set is a submanifold *because* the
charge is conserved. Ours is a single prescribed level $\{-g(u,W)=\hat E\}$, and
that argument is unavailable once the charge is not conserved. We were adjacent
to assuming the structure carries over. It does not, and the manuscript now says
so — which is precisely why Paper I's existence theorem goes by the direct method
on an augmented state (Grönwall + Filippov–Cesari) rather than by restricting a
functional to a constraint submanifold.

**Meena–Zawadzki 2024.** Their Theorem 5 has **three** conditions, not one; we
test only the necessary "dilation constant along fibres" and now say so, with a
caveat that their setting is Riemannian-on-Riemannian while ours has a Lorentzian
total space. Two consequences we had been asserting without support are now
theirs: **Prop. 1**, horizontal distribution totally geodesic ⟺ dilation
fibre-constant; and **Cor. 7**, the rescaling to a Clairaut *Riemannian*
submersion exists only when $\lambda$ descends to the base. That is the precise
reason the non-stationary optical structure does not live on a quotient.

**Axioms 15(4) 267 (Yamabe flow).** Its Eq. (1), $G\to e^{-2f}G$, is the
metrisation we had identified with the conformal factor — so that identification
is standard in that literature, not an ad hoc device. It also supplied a
classical framing for a result we had only verified: conformal Schwarzschild is
not Einstein for generic $A$, and **Birkhoff's theorem** says why — a spherically
symmetric vacuum solution is necessarily static, so no running conformal factor
attaches to Schwarzschild without sourcing matter.

**Lecian, 2-planes solitons.** Cites Gibbons on *energy-independent* Randers-type
Finslerian metrics and the Gauss curvature of equatorial Schwarzschild sections.
That sharpened our own threshold: in the energy-independent (null) case $K<0$
throughout the exterior and **no threshold exists**. Ours appears only because the
optical metric carries the rail energy — the leading far-field sign is that of
$-(2\hat E^2-3)$, which cannot change in a family with no $\hat E$ to vary.

**Giannoni–Masiello–Piccione 2000 (Morse for massive particles).** Their
variational data are ours minus the control: fixed launch event, arrival on a
timelike curve, arrival-time functional, pointwise velocity constraint. Their
Theorem 1.2 gives **Morse index = geometric index** at a nonconjugate endpoint,
hence index zero ⟹ local minimum — the step Proposition (no conjugate points)
needed. Their Morse relations also locate our Maxwell set: the count runs through
the Poincaré polynomial of the path space, and our annulus organises connecting
extremals by winding number. We do not verify their pseudo-coercivity and make no
global count.

**The five roles of $\psi=A'/A$.** Collected in §2.2: conformal factor of the
selector; rail drift along a geodesic; minimum thrust; failure of the dilation to
be fibre-constant (hence of Clairaut, of the horizontal distribution to be
totally geodesic, and of descent to a quotient); and failure of the weight to be
parallel (hence of Einstein). Five classically independent conditions coincide
here and degenerate together at frozen $A$.

## 4. Discipline imposed throughout

- **Protocol 1** (E)/(C)/(A): exterior extremal / limiting contact / analytic
  continuation with no optimality claimed. Applied in text *and figures*.
- **Protocol 2** endpoint: fixed launch, fixed target, free arrival; comparisons
  at $J_t=AJ_\eta$. §3.4 names its two protocols (P1) and (P2) explicitly.
- **Protocol 3** nomenclature: rail charge ≠ energy; $J$ costate ≠ angular
  momentum; $r_+$ = seed Kerr null surface; no capture/plunge/bounce inside $2M$;
  a table keeping $E,\hat E,E_{\rm eff},J,J_{\rm eff}$ distinct; units and error
  conventions.
- **Protocol 4** evidence levels: (I) exact symbolic identity = proof;
  (II) independent numerical agreement ≠ proof; (III) recovery of established
  limits; (IV) computer-assisted with stated grid. Plus an explicit remark that
  **no experimental validation exists** for this class of result, and why.
- Default term is **"controlled-rail extremal"**; minimiser/optimal/fastest
  reserved for statements with a certificate.
- Higher-genus terminology neutral: "genus-two dilogarithm" is a *descriptive
  shorthand*, never a classification.
- Title no longer says "trichotomy".

## 5. Reproducibility

- `cap_full.py` restored verbatim from the prototype, with a self-test. The CAP
  runner now certifies cells (`tag='S'`, 0.2 s/cell).
- `provenance/make_provenance.py`: **one command** runs the validation scripts,
  parses slopes and covariance $\sigma$, reads the $\varepsilon$-window out of
  the source, hashes every script, and emits the LaTeX macros the manuscript
  uses. Values can no longer drift.
- `provenance/make_crossrefs.py`: the response's cross-references are read out of
  `paper2.aux`. (A first draft cited §2.2/Protocol 2.2/Proposition 3.1 where the
  manuscript says 2.1/Protocol 1/Proposition 2.)
- `MANIFEST.tsv`: 16 artefacts, each → one command, one checksum. Environment
  lock verified against the live environment and extended with the CAS version.
- Sub-window refit demonstrating the finite-window bias: $2.118\to2.068\to2.043$.

**Verification suite** (`paper2/verification/`), all passing:

| script | checks | subject |
|---|---|---|
| `verify_cusp_corner.wls` | 22 | shape at $r_e$ |
| `verify_marginal_hamilton.wls` | 23 | dynamics at $r_e$ |
| `verify_exterior_separatrix.wls` | 14 | flagship separatrix |
| `verify_submersion_link.wls` | 14 | optical submersion, dilation |
| `verify_soliton_rail.wls` | 10 | soliton drift identity |
| `verify_tk_ckv_rail.wls` | 9 | TK conformal Killing, drift $=\varepsilon$ |
| `verify_myers_optical.wls` | 9 | optical curvature, Myers |
| `verify_thrust_bound.wls` | 13 | minimum thrust |
| `verify_conjugate_maxwell.wls` | 9 | conjugate points, Maxwell set |
| `verify_submersion_link.wls` | 14 | optical submersion, dilation |
| `verify_soliton_rail.wls` | 10 | soliton drift identity |
| `verify_tk_ckv_rail.wls` | 9 | TK conformal Killing, drift = epsilon |
| `verify_weighted_seed.wls` | 9 | the weighted seed correspondence |
| `verify_dilation_fibre.wls` | 4 | the Clairaut necessary condition |

## 6. Paper I

- The incorrect PMP sentence fixed in **both** variants.
- `paper1.tex`: Theorem I.2's admissible class stated in five parts (control
  domain, extended state and clock, admissible arcs and closure, endpoint
  protocol, excluded degeneracies).
- **The JMP variant is the current Paper I**; it already had hypotheses H1/H2/H3
  and a Filippov–Cesari proof.
- Vaidya boundary qualification added to **both** abstracts: the problem ends at
  $r=2m(v)$, the calculation gives an exterior approach and contact threshold,
  and $m'(v)<0$ ingoing is a formal continuation.

## 6bis. Paper I: what was ported back, and the frozen submitted text

About half the mathematics developed while answering this review is **general**
and belongs to Paper I, not to the conformal-Kerr manuscript where it happened to
be written. Two versions now exist and must not be confused.

| | |
|---|---|
| `paper1/submitted_JMP_2026-08/` | **frozen** — the text exactly as sent to JMP (from `91749fe`), with a README naming its one known defect and directing revisions elsewhere |
| `paper1/paper1_JMP.tex` | working copy, 54 pp, the version archived on Zenodo |

This lets the CQG response say that the conceptual corrections to the minimum
principle are public and citable, without disturbing the JMP process. If JMP
returns a report, revise the working copy.

**Ported**, as a new subsection after the selector hierarchy:

- **Lemma (rail drift and the control force).**
  $\frac{d}{d\tau}(-g(u,W))=-g(a,W)-\tfrac12(\mathcal L_Wg)_{\mu\nu}u^\mu u^\nu$.
  At the Killing rung the right side vanishes and the rail is free; below it the
  same identity **quantifies the cost of each rung**. Vaidya:
  $-m'(v)(u^v)^2/r$. Conformal-Killing rung: $\psi$ itself.
- **Sharp minimum-thrust bound**
  $|a|\ge|\tfrac12(\mathcal L_Wg)(u,u)|/\sqrt{\hat E^2-|W|^2}$, attained when
  $a\parallel W_\perp$. Diverges as an inverse square root at the freezing
  surface — **but not merely because $|W|\to0$**. The two boundaries of the
  control domain are of different kinds: loss of compactness at finite cost
  versus divergence of cost.
- **An obstruction theorem** for why the time-dependent optical metric cannot be
  stated on a quotient the way Perlick's can. The optical projection is
  horizontally conformal with dilation equal to the Perlick factor, and
  conformality is *forced* by Fuglede–Ishihara; but descent needs fibre-constant
  dilation (Meena–Zawadzki), equivalently umbilical fibres with gradient mean
  curvature (Bishop), and that fails away from a Killing selector. This converts
  a presentational choice into a theorem, and it was the weakest part of Paper I.
- **A remark on two structures deliberately not assumed**: the constraint
  submanifold of Caponio et al. (its construction uses conservation) and the
  Palais–Smale property of GMP 1998 (its constraint is flow-preserved).

## 7. Bibliography

Added: Haws–Kiser 1995; Taş (arXiv:2512.08776, **Dec 2025**, not 2026);
Giannoni–Masiello–Piccione ×3 (timelike Fermat 1998; Morse for massive particles
2000; **conformally stationary** 2001 — closest to our frozen-$A$ regime);
Caponio–Corona–Giambò–Piccione 2024 (fixed energy, affine Noether charge —
closest published variational structure to the controlled rail, with the
difference stated); Bishop 1972; Fuglede 1978; Ishihara 1979; Baird–Wood 2003;
Meena–Zawadzki 2024; Myers 1941; Hamilton 1986; Böhm–Wilking 2008.

**On R2's misattribution.** The paper at JMP 38(12), 6367–6381,
`doi:10.1063/1.532217` is by Giannoni, Piccione and **Verderesi**, not Masiello
(CrossRef, Semantic Scholar, and the AIP page R2 themselves link). We do **not**
say so in the response. We cite it correctly and add the genuine
Giannoni–Masiello–Piccione works; the correction makes itself.

**On the "trivial application of T2" claim.** Both preprints were obtained and
read: `10.13140/RG.2.2.31647.21921` (contains T2) and
`10.13140/RG.2.2.35935.37281` (the soliton/geodesic-paths one, whose Ricci-soliton
definition produced the identity of §3.4). T2 is the Fuglede–Ishihara dichotomy. The response concedes the real
relation (§3.6 above), then declines the stated one, narrowly: T2 gives the
existence of a conformal factor, not its value; and its hypotheses — by Bishop's
theorem, **reference [14] of that preprint itself** — are equivalent to totally
umbilical fibres, which is exactly what the non-stationary problem gives up.

**Not incorporated, deliberately**: the curvature-operator strand as a *method*.
The reason is computed, not asserted — see §3.6.

---

## 8. Deliverables

| File | State |
|---|---|
| `paper2.pdf` | 63 pp, 0 unresolved refs, 0 floats too large |
| `paper2_highlighted.pdf` | via `latexdiff --type=CFONT` |
| `response_to_referees_CQG.pdf` | 9 pp |
| `paper1_JMP.pdf` / `paper1.pdf` | 52 / 35 pp |
| `.zenodo.json` | 1.3.0; tag `v1.3.0` pushed |

**`latexdiff` needed three repairs** to compile, and will again: `ulem`
underlining breaks across `\emph`, so use `--type=CFONT`; float-context markers
land before `\hline` inside tabulars; one `\multicolumn` row carries a stray
marker.

## 8bis. Traps that cost real time

- **`grep -c Overfull` does not catch `Float too large for page`.** Different
  warning class. Check both.
- **Wolfram Language terminates a statement at a newline** if what precedes is
  syntactically complete. Writing
  `f[x_] := a + b x^2` *newline* `- c x^3;` **silently drops the cubic term** —
  no error, no warning, wrong answer. It cost a full debugging cycle here, and
  the symptom was a scan that came out exactly inverted. Always end the
  continued line with the operator.
- **Mathematica comments nest**, so `(*)` inside a `(* ... *)` block opens a
  comment that is never closed.
- A generated `\newcommand` that already contains `$…$` must not be used inside
  `$…$`.
- `NSolve[... && 2 < r < 10^6, r, Reals]` can return unevaluated; for a
  polynomial, take exact roots and filter.
- Declaring a symbol positive in SymPy can make a degeneration locus vanish.

## 9. Still open

1. **Tag every numbered claim** with its evidence level — the three main results
   carry theirs; the rest do not yet.
2. **Source bundle** for CQG (`.tex` + figures + `.cls`) not yet assembled.
3. **The bibliography has never been audited** for whether cited works say what
   we attribute to them — a caveat inherited from Paper I, still true.
