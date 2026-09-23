# Sources for the bibliographic audit

The PDFs in this directory are **not committed** — they are published papers
under publisher copyright (APS, IOP, AIP, Springer, Elsevier). `.gitignore`
excludes `Fonti/*.pdf`. This file is tracked and records what was collected and
what was checked, so the audit is reproducible from the DOIs alone.

**Naming convention.** Every file is `Author-Year-topic.pdf`: first author's
surname (or concatenated surnames for two or three authors, `EtAl` beyond that),
four-digit year, then a short hyphenated topic. No download identifiers, no arXiv
numbers. Preprints with no year on the title page use `undated` rather than a
guessed one. The directory holds **54 files, one per work, no duplicates** — a
`sha256` pass removed 44 duplicate downloads on 2026-08-23.

Verification scripts referenced below are in `../paper2/verification/`; the
machine check of the bibliography against the DOI registry is
`../paper2/provenance/audit_refs.py`.

---

## Why this exists

Two kinds of error hide from anyone reading only the compiled PDF:

1. **A citation that resolves to the wrong paper.** `Perlick1991` carried
   `10.1063/1.529473`, which is Bengtsson, JMP 32, 3158 — not the paper the whole
   programme extends. `MeenaZawadzki2024` pointed at Liu–Zhao.
2. **A source that does not say what we attribute to it.** Only reading finds
   these. Three so far.

---

## Audit status

Legend — **✓** confirmed · **✗** correction made · **+** produced new material

| Source | Read | Outcome |
|---|---|---|
| Perlick 2000, *Ray Optics* (book) | full ToC | **✗** the monograph is entirely about **light rays** — no timelike or brachistochrone chapter — yet Paper II cited it for "the brachistochrone was formulated by Perlick". Split from the 1991 paper and now cited by section for the ray-optical/Fermat machinery only. Index confirms independently: no *brachistochrone* entry, but *arrival time functional* (159, 164) and *conjugate point* (105, 165, 175, 178, 203, 207) — our topics, done there for light |
| Perlick 1991 | core | **✓ +** two optical factors are exactly his Eqs. (42), (46); his Rindler example already shows the branch asymmetry; he uses "brachistochrone" for *stationary* travel time |
| Bishop 1972 (via Singh–Gupta) | thm + proof | **✗** we had added a completeness hypothesis the theorem does not require |
| Neishtadt 2014 | intro + ToC | **✗** it is about resonance passage, not separatrix crossing; citation narrowed |
| Neishtadt–Vasiliev 1999 | intro + refs | **✗ +** *not* the paper we cite (3D volume-preserving, and they state the two problems "cannot be reduced to one another"). Its ref. [3] confirms our target exactly, and dates it **1986** — our key said 1987. Supplied the two companion jump-formula papers |
| Gibbons–Werner 2008 | core | **✓ +** their Eq. (8) reproduces our null limit identically; Gauss–Bonnet gives a second topological reading of the Maxwell set |
| Meena–Zawadzki 2024 | core | **✗ +** Thm 5 has *three* conditions, not one; Prop. 1 and Cor. 7 supplied two results we had been asserting |
| Caponio et al. 2024 | core | **+** the semi-holonomic constraint; and a risk closed — their submanifold argument uses conservation, which we do not have |
| GMP 1998 | core | **+** third-party reading of Perlick; Palais–Smale is for a flow-preserved constraint, not ours |
| GMP 2000 | core | **+** Morse index = geometric index, the step our local-minimality claim needed |
| GMP 2001 | core | **✗ +** we had called it "closest to our regime"; it is about *null* geodesics and rests on the Killing property. Light-convexity frames the Maxwell count |
| Lecian, curvature eigenvalues | full | **+** the weighted-seed correspondence; T8 draws the boundary |
| Lecian, T2 preprint | full | analysis of the "trivial application" claim |
| Lecian, solitons/geodesics | half | the λ + Ric(u,u) identity |
| Lecian, Yamabe (Axioms) | intro | **+** metrisation is standard; Birkhoff explains why conformal Kerr is not Einstein |
| Lecian, Cauchy (Axioms) | intro | **+** independent route to our control-domain boundary; rail flow incomplete where geodesics are not |
| Lecian, 2-planes | intro | **+** energy-dependence localises our threshold as new |
| Fuglede 1978 | §§5–8 + intro | **✗✗ +** what we called "the Fuglede–Ishihara dichotomy" is his **definition** of semiconformal (§5), not a theorem — the argument built on it was circular. The real theorem (§7) is *harmonic morphism ⟺ semiconformal + harmonic*, which says what semiconformality does **not** give. Also: his machinery is Riemannian, our total space Lorentzian |
| Ishihara 1979 | intro + §§1–3, Thm 5.1 | **✓ +** his Thm 5.1 *is* a genuine dichotomy, but hypothesised on *preserving harmonic functions*, not on semiconformality — our garbled version had swapped the two. His §1 criterion (submersion harmonic ⟺ minimal fibres) then **settles** the harmonic-morphism question for us: fibres are 1d, not geodesic, obstruction = static-observer thrust |
| Doran 2000 | full | **✓ +** our azimuthal shift is exactly his Eq. (21) composed with BL→EF, verified symbolically; his chart is built on *free-falling* observers, which is why it cannot rescue the control domain |
| Kodama 1980 | core | **✓ +** claims confirmed for our normalisation; his own K differs in *norm* (his Eq. 3.6), now noted; **he states K is timelike outside the apparent horizon** — third independent route to our control-domain boundary |
| Böhm–Wilking 2008 | thm 1 + intro | **✗ +** their Thm 1 needs a *compact* manifold and *2-positive* operator; and their intro dates the results — Hamilton **1982** is the 3d one, 1986 the 4d one, so we were citing the wrong Hamilton for a 3d base |
| Cary–Escande–Tennyson 1986 | core §§I–II | **✗ +** corrected my own wording of a day earlier: the jump is *deterministic at lowest order* (area of the final lobe), phase dependence only at $O(\varepsilon\ln\varepsilon)$. Supplies the mechanism — the frozen period diverges on the separatrix, so $\delta=\varepsilon T_0$ is large however small $\varepsilon$ |
| Giannoni–Piccione 2002 | §1 | **✓ +++** the best novelty statement in the audit. Their trial paths are *timelike* like ours; they impose $g(a,Y)=0$ (frictionless slide), which by $\dd_\tau(-g(u,W))=-g(a,W)-\tfrac12(\mathcal L_Wg)(u,u)$ leaves the charge drifting at the Killing defect. They then state that non-stationary brachistochrones may not solve a *second-order* equation at all, multipliers being eliminable "unless in the stationary case". Our rail is the other branch: hold the charge, pay the thrust |
| Hamilton 1982 | Thm 1.1 | **✓** confirmed at the source, not just via Böhm–Wilking: *compact* 3-manifold, strictly positive Ricci ⟹ admits constant positive curvature |
| Hamilton 1986 | Thm 1.1 | **✓** *compact* four-manifold, positive curvature operator ⟹ diffeomorphic to $S^4$ or $\mathbb{RP}^4$. Both Hamilton attributions now rest on the primary sources |
| Chen 1977 | §1.1 + refs | **✗** he attributes the **shuffle** to **Ree** (his ref. [57] = Ann. of Math. 68 (1958) 210), not to himself; we had credited Chen. Ree added to the bibliography |
| Ashtekar–Krishnan 2003 | §§I–II, Def. 1 | **✗ +** *space-like* is part of their **definition** of a dynamical horizon, and they exclude the timelike case only under the DEC. Our unconditional claim for ingoing Vaidya needed $m'(v)>0$: on the tube $g(T,T)=4m'(v)$, so accretion gives a dynamical horizon and evaporation their separately-treated timelike tube |
| Misner–Sharp 1964 | §I | **✓** their Eq. (1.11) defines the mass function as "an appropriate *total energy* of each fluid sphere" — the attribution is right |
| Hayward 1996 | abstract + §II | **✓ +** states our claim verbatim: "the conserved Kodama current has charge $E$", with $E=\tfrac12r[1-g^{-1}(\dd r,\dd r)]$ (his Eq. 4); his trapping criterion $E\gtrless r/2$ puts the marginal sphere at $r=2m(v)$ |
| Lindquist–Schwartz–Misner 1965 | §§I–III | **✗** they analyse the **outgoing** (retarded) Vaidya form, where positivity of the radiated energy forces $\dd m/\dd u\le0$ — the *opposite* mass-rate sign to the ingoing accreting case we cite them for |
| Booth 2005 | §§1–2 | **✓ +** a review, appropriate for the bundle; and he flags (with Ashtekar–Krishnan, independently) that "apparent horizon" in numerical relativity is *not* the Hawking–Ellis term. Paper I now says which sense it uses |
| Abreu–Visser 2010 | §§I–III | **✗ +** crediting Kodama himself: the Kodama vector does *not* in general reduce to the Killing vector in a static spacetime, only to something **parallel** to it. We wrote "reduces to $\xi$". In fact $K=\partial_t/(EL)$, equal iff $g_{tt}g_{rr}=-1$ — true for Schwarzschild and Vaidya, false for a general static interior |
| Ashtekar–Krishnan 2003 | — | pending |
| Hayward–Mukohyama–Ashworth 1999 | §§1–2 | **✓** their $k=*\dd r$ is our Kodama vector, "a dynamic analogue of a stationary Killing vector"; and they state plainly that "Misner & Sharp originally defined $E$", confirming that attribution |
| McVittie 1933 | §§I–II | **✓** confirms Kaloper at the source: his Eq. (1) is Schwarzschild in *isotropic* coordinates, and he builds the solution by requiring the matter to be "at rest" in his system — zero coordinate velocity and momentum, i.e. the "no accretion" property |
| Nolan 1998 | — | pending |
| Kaloper–Kleban–Martin 2010 | §§I–III | **✗ +** McVittie is **not** a conformal rescaling — its lapse carries no scale factor and its mass is constant ("there is no accretion"). Thakurta belongs to the *conformal* class instead. They also show Nolan's central assertion is **incorrect**, and call that literature "riddled with basic errors" |
| Mello–Maciel–Zanchin 2017 | abstract + §§I–II | **✓✓ ++** they analyse **our metric**. Published result: for increasing unbounded $a(t)$ Thakurta "does **not** describe a cosmological black hole" — independent confirmation of our reading. Also: Thakurta is a *generalised* McVittie at the boundary case $\dot m/m=\dot a/a$, and is **not** Sultana–Dyer (factor in $\eta$ vs EF advanced time), a confusion they flag |
| Vaidya 1951 | §§1–3 | **✗** his paper is the **radiating** star — outgoing, emitting. The ingoing accreting form we use is its time-reverse; standard, and standardly called his, but not what he computes |
| Kerr 1963 | — | attribution of the seed metric; standard, no claim beyond it |
| Carter 1968 | abstract + §1 | **✓ +** exactly our claim: the fourth constant comes "from the unexpected separability of the Hamilton–Jacobi equation". Also corroborates the twist computation — the principal null congruences "have nonzero rotation (except when $a$ vanishes)" |
| Taş 2025 | abstract + §1 | **✓ +** description accurate; but a real overlap to declare — his Schwarzschild reduction *is* our arrival-time optical metric, and he reaches transverse-variation/stability too. They agree because in a static spacetime the controlled rail *is* the fixed-energy geodesic |
| Giannoni–Piccione–Tausk 2002 (= the 1999 preprint) | §1 | **✗ +** it is the **travel**-time theory, ours the arrival branch — they say the two problems are "essentially different", so the branch must be matched where we appeal to their Morse theory. Also carries the non-stationary obstruction three years before GP 2002, there as the stated reason for restricting to stationary |
| Piccione–Tausk 2000 | — | collected as a substitute for GPV 1997, not cited |
| Ichikawa 2023 | — | collected as adjacent to the genus-2 cluster, not cited |
| D'Hoker–Enriquez–Schlotterer–Zerbini 2026 | — | pending |
| Caponio–Javaloyes–Sánchez 2024 | abstract + ToC | **+** their *wind* Finslerian structures exist precisely for $|W|\ge1$ — the regime past our control-domain boundary. §4 causal $K$, §5 arbitrary $K$, §8.3 $K$-horizons. So the geometry continues past $r_e$; our compact-control problem does not |
| Baird–Wood 2003 | abstract | **✓** independently confirms the corrected Fuglede statement verbatim: harmonic morphisms *are* harmonic maps that are semiconformal; and it lists Killing-field submersions as examples |
| Sultana–Dyer 2005 | §1 | **✓ ++** supplies the *published* necessary-and-sufficient criterion for our own terminology: a CSLS is a conformal Killing horizon iff the congruence twist vanishes there. Applied to TK it **fails**, twist $=-a/2M$ equatorially, conformally invariantly — so $r_e$ is not a horizon, independently of $A(\eta)$ |
| Born–Fock 1928 | §§1–2 | **✗** category error: their theorem is *quantum* — energy operator, Schrödinger equation, transition probabilities — and we cited it for a classical asymptotic expansion. Their own p. 165 attributes the classical action-variable statement to **Ehrenfest** |
| Haws–Kiser 1995 | full | **✗** it is a *pedagogy* article — a Mathematica package for racing curves, pre-calculus upward — not a review of "the classical problem and its geometry". Its one substantive part, §3, is the brachistochrone **with kinetic friction**, i.e. the elementary case where the conserved energy fails |
| D'Hoker–Hidding–Schlotterer 2025 | §§4.5–4.7 + abstract | **✓ ++** vindicates our hedge with the authors' own words: from genus two the **meromorphic** sector (where our third-kind letter lives) admits no two-point meromorphic analogue "without additional marked points" (§4.5), and closure under primitives is offered as *evidence*, not proof (§4.6) |
| D'Hoker–Schlotterer 2024 | via [69] of DHHS | **✓** resolved without the PDF: DHHS cite it as where the higher-genus Fay identities "and their proof" are given (§5). Corrected my own sentence calling them "still being developed" |
| Baune et al. 2024 | abstract + §§5.5–7 | **✓ +** author names confirmed against the title page (the earlier Konstantin/Egor fix holds). Their test curve is genus-2 real hyperelliptic with real roots — our class. They state plainly that convergence is *argued*, no error estimates, code unreleased, prototype: strong support for our hedge, and a clean contrast with our archived theta-nome evaluation |
| Giannoni–Piccione–Verderesi 1997 | §I | **✓ ++** their $U_k=\{\beta<k^2\}$, i.e. $|Y|<k$ (their Eqs. 2–3), **is our control domain** — Lemma 2's $\hat E>|W|$, with their strict inequality matching where our slice degenerates. A thirty-year precedent for the boundary referee 1's Major 1 is about. Also a third independent source for "the two variational problems are essentially different" |
| Dyer–Honig 1979 · Sultana–Dyer 2004 | via SD2005 refs [7],[8] | **✗** the CSLS/CKH *terminology* is Dyer–Honig 1979 (JMP **20** 409) and the twist *criterion* is Sultana–Dyer **2004** (JMP **45** 4764) — SD 2005 restates both. We had credited 2005 for both. Found by a NotebookLM sweep, verified on the PDF's reference list and CrossRef |
| Randers 1941 | full | **✓ +** attribution correct; and he *names* the gauge freedom $b\mapsto b+d\phi$ (his "$k$ transformation", distinguished from Weyl's) that our Doran remark relies on, and states the exactness criterion our non-rotating restriction turns on |
| Kovner 1990 | core §§1–2 | **✗✗ +** *not* a context citation — his abstract covers "arbitrary, stationary, and **nonstationary** metrics, for **massive** and massless particles". We had bundled him under the stationary case, a live priority risk. Now distinguished on four grounds: his extremals are *geodesics*, his only constraint is the mass shell, he asks for extremality not minimality, and his one genuinely non-stationary application is *perturbative lensing* — a static lens rippled by a weak gravitational wave (§§IV–V), asking where images form, not what the time-optimal constrained worldline is |
| Gibbons et al. 2009 | core §§1–2, 3.4 | **✗ +** their triality is for *null* geodesics, so our "reproduces Randers/Zermelo" needed narrowing to the form; **+** their Finsler condition $|b|_a<1$ fails exactly on the ergosurface, and the Randers data are conformal invariants — a fourth independent route to our control domain |
| Bao–Robles–Shen 2004 | §0 | **✓ +** attribution correct; their "maximal domain" $|W|_h<1$ is the Zermelo form of the same number as Gibbons' $|b|_a<1$, giving our control domain the reading *the wind outruns the ship*; and their Thm 3.1 places constant flag curvature on the homothety rung, strictly below our selector |
| Brown–Levin 2011 | abstract + §1 | **✗** theirs are the **multiple** elliptic polylogarithms — several variables, configuration space of $n{+}1$ points — which "generalize the classical elliptic polylogarithms". Our $W_{jk}$ is single-variable, length and weight two, so it belongs to the **classical** class (Beilinson–Levin, Zagier). They also build on **Chen's** reduced bar construction, confirming that citation |
| Myers 1941 | full | **✗ +** our argument was a non sequitur — he *proves* sectional decay is not enough ($S^2\times S^2$, §4). Retested on Ric: $R_{\rm opt}<0$ pointwise at every energy, so the hypothesis fails by sign. His Lemma needs no completeness |
| Filippov 1962 | full | **✗ +** bib entry was broken (DOI orphaned outside it); every Thm. 1 hypothesis re-checked and met with margin; his §III *proves* our no-relaxation claim, §V is the counterexample behind our remark, §IV explains why the domain is compact |

## Referee-proposed reading

Taken from the reports themselves (`Response/PaperII/CQG response/`), not
reconstructed. The two referees differ sharply on this point.

**Referee 1** (`report_CQG_116884.pdf`) proposes *nothing*. Major comment 12 states
that the reference coverage is "broadly current and appropriate" and that the
remaining issue is "a matter of positioning rather than a request for additional
citations". Nothing in that report asks for a source we do not already cite.

**Referee 2** (`CQGReviewerstask2026I.pdf`) gives an explicit numbered list under
*Items of bibliography to be presented in the introductory text*, plus references
[1]–[11] in the body.

| Proposed by referee 2 | Here? |
|---|---|
| 1) Perlick 1991, JMP **32**(11) 3148–3157 | ✓ read |
| 2) *Giannoni, Masiello & Piccione* 1997, JMP **38**(12) 6367–6381 | ✓ **read** — the referee's author list is wrong: the title page gives Giannoni, Piccione, **Verderesi**, which is what our bibliography already had |
| 3) *Giannoni, Piccione & Tausk* 2002, "The arrival time brachistochrones", *Class. Quantum Grav.* | ✓ read — but the referee has conflated two papers. "The arrival time brachistochrones in general relativity" is **Giannoni & Piccione**, *J. Geom. Anal.* **12** (2002) 375 (`10.1007/BF02922047`); the Tausk paper is the *travel*-time Morse theory, DCDS **8** 697. We hold and cite both |
| 4) Taş 2026, `arXiv:2512.08776` | ✓ read (v3) |
| 5) Giannoni–Piccione–Tausk 1999, `arXiv:math-ph/9905007` | ✓ read (= the DCDS paper above) |
| [5] Haws & Kiser 1995, Amer. Math. Monthly **102** 328–336 | ✓ read |
| [9] Hamilton 1986, JDG **24** 153–179 | ✓ read |
| [10] Böhm & Wilking 2008, Annals **167** 1079–1097 | ✓ read |
| [11] Lecian, *Eigenvalues of the curvature of Einsteinian weighted solitons* | ✓ read |
| Myers 1941, Duke Math. J. **8**(2) 401–404 | ✓ read |
| Lecian, *New Geodesics 2-planes solitons connected at infinity* | ✓ read |
| Lecian, *Geodesics completeness and Cauchy hypersurfaces*, Axioms **14**(12) 896 | ✓ read |
| Lecian, Axioms **15**(4) 267 (Yamabe) | ✓ read |
| Lecian, *Pseudo-Riemannian solitons after umbilicity conditions* (T2) | ✓ read |
| Lecian, *Generalized Schwarzschild solitons, spherically-symmetric weights* | ✓ read |

So: **all sixteen items are held and read.**

### Round 2 (decision of 16 September 2026)

The referee numbering is **inverted** with respect to Round 1: the technical
report is Referee 2 here and asks for nothing; the bibliographic referee is
Referee 1. Reports in `Response/PaperII/CQG response/round2/`.

Referee 1 lists four citations, but they are **three works**: the arXiv preprint
and the EPJ Plus article carry the same title and are the same paper.

| Proposed by referee 1, round 2 | Here? |
|---|---|
| Chanda, *More on Jacobi metric: Randers–Finsler metrics, frame dragging and geometrisation techniques*, `arXiv:1911.06321` | ✓ held — `Chanda-2024-Randers-Finsler-frame-dragging-geometrisation.pdf` |
| Chanda, same title, *Eur. Phys. J. Plus* **139**(11) 2024, `10.1140/epjp/s13360-024-05775-y` | ✓ **same work as the line above** — the file held is v13 (5 Nov 2024), which is the published version; §5.2.1 treats the Kerr metric |
| Chanda, Gibbons, Guha, Maraner & Werner, *Jacobi-Maupertuis Randers-Finsler metric for curved spaces and the gravitational magnetoelectric effect*, *JMP* **60**(12) 122501 (2019), `10.1063/1.5098869`, `arXiv:1903.11805` | ✓ held — `ChandaEtAl-2019-Jacobi-Maupertuis-Randers-Finsler-magnetoelectric.pdf` |
| Lecian, *Experimental validations of the geodesics 2-plane from GR pseudospherical cylinders*, `10.13140/RG.2.2.24140.58246` | ✓ read — `Lecian-2026-geodesics-2-plane-experimental-validations.pdf`, 13 pp, Sapienza, 7 Sep 2026. Distinct from *New Geodesics 2-planes solitons connected at infinity* also held here |

Two DOIs were confirmed against Crossref rather than copied from the report; the
ResearchGate identifier resolves to no Crossref record at all, and the file was
obtained from a logged-in session.

**What the Chanda papers actually contain.** Chanda–Gibbons–Guha–Maraner–Werner 2019
builds the Jacobi–Maupertuis Randers–Finsler metric for a **stationary** spacetime at
**fixed conserved** $p_0$, eq. (16):
$ds_J=\sqrt{\frac{p_0^2-(mc)^2g_{00}}{g_{00}}\gamma_{ij}dx^idx^j}+p_0\frac{g_{0i}}{g_{00}}dx^i$
with $\gamma_{ij}=-g_{ij}+g_{0i}g_{0j}/g_{00}$, valid where the Randers convexity
condition $\sqrt{a^{ij}b_ib_j}<1$ holds; eq. (24) is the explicit Kerr case. Two
facts in it bear directly on the referee's demand. First, their §1 states that the
section-4 result "agree[s] with a remark in" Perlick, *The brachistochrone problem
in a stationary space-time*, JMP **32** (1991) 3148 — the paper Paper II already
builds on (reference [9] in the arXiv version, [16] in the published one). The
agreement is with Perlick's **Proposition 3.3**, the *free-fall* statement, not
with his Proposition 3.2, which is the $t$-brachistochrone: see the correction
below. Second, the construction presupposes a conserved $p_0$ and a *geodesic*;
our rail is neither free nor at conserved charge once $A$ runs, since
$\dd\Ehat/\dd\tau=\varepsilon$.

**Correction (16 Sep 2026), and the paper that settles the referee's demand.**
Chanda's eq. (16) reduces the *free geodesics* at fixed energy. The arrival-time
brachistochrone is a different functional: eliminating proper time from the
constraint gives
$F_T=\sqrt{\frac{E^2}{f(E^2-f)}h_{ij}dx^idx^j}-\omega_idx^i$, whose quadratic part
is exactly $\Lambda^2h$ of the manuscript's eq. (12), whereas Chanda's is
$\frac{E^2-f}{f}h$. The two differ by the position-dependent factor
$E^2/(E^2-f)^2$, so they do **not** share geodesics — in Schwarzschild at $E^2=2$
the circular extremal sits at $r=2\sqrt2M$ for the arrival cost and at
$r=(1+\sqrt5)M$ for Jacobi–Maupertuis.

This is already in **Perlick 1991, p. 3153**. His **Proposition 3.2** states that a
curve is a $t$-brachistochrone of specific energy $e^C$ iff it minimises
$\tilde S_C=\int[\sqrt{\bar h_C(\xi',\xi')}-\psi(\xi')]\,ds$ — a Randers functional,
which he himself compares to "a charged particle moving in a magnetostatic field".
His **Proposition 3.3** is the free-fall counterpart with a different metric
$\hat h_C$, and his comparison of Props 3.1–3.3 shows the brachistochrone problem
is equivalent to free fall only *in some other stationary spacetime*, with the
modified potential of his eq. (50). So the Randers form of the brachistochrone was
published in 1991, in the paper this manuscript already rests on.

Chanda 2024 (the EPJ Plus paper) adds, at eqs. (3.2.4)–(3.2.5): "optical metrics are
not Jacobi metrics for null curves, and their similarity for Riemannian static
metrics … is merely coincidence." So the identification of our optical/Perlick
structure with a JMRF metric is one the referee's own preferred source denies in the
stationary rotating case, and grants only in the static Riemannian one.

**What the Lecian preprint actually contains, on the point the referee raises.**
§7, *Analytical details of geodesics deviations … the geodesics separations*, is the
textbook Jacobi equation: the covariant separation velocity (15), its second
covariant derivative (16), and the geodesic-deviation equation
$D^2\xi^\mu/d\tau^2=R^\mu{}_{\nu\alpha\beta}u^\nu u^\alpha\xi^\beta$ (17). The
background is the pseudospherical-cylinder soliton
$ds^2=-A(r)dt^2+B(r)dr^2+r^2(d\theta^2+\sinh^2\theta\,d\varphi^2)$ of eq. (1) —
hyperbolic angular sector, static weights. §8 lists measurement techniques
(gyroscopic precession, VLBI photon-ring profiling, gravitational-wave strain).

No equation in it involves Kerr, a Randers–Finsler structure, a brachistochrone or
a control constraint, so it supplies a *method* to apply — write the Jacobi equation
on our geometry and evaluate the separation — not a result that bears on our
spacetime. Two internal notes for whoever cites it: §5 repeats §4 verbatim
(eqs. (9)–(11) are eqs. (6)–(8)), and reference [16], "The LISA Cosmic Carrier Team,
*Phys. Rev. D* **111**(6) 064012 (2025)", should be verified before being relied on.

Two items previously recorded here were my own reconstruction and are wrong:
Cesari 1983 was never proposed by either referee, and neither was Filippov 1962 —
both entered our bibliography on our own initiative. They are listed under
*Still missing* only because we could not obtain them, not because a referee asked.

Also collected, not proposed but adjacent: Piccione–Tausk, *Variational aspects of
the geodesic problem in sub-Riemannian geometry* (`math/9911215`) — the closest
free substitute for Giannoni–Piccione–Verderesi.

## Still missing

- **Cesari 1983**, *Optimization — Theory and Applications*, Springer. Book, no
  DOI. Load-bearing for Paper I's existence theorem.
- **Thakurta 1981**, Indian J. Phys. **55B**, 304. No DOI, journal not indexed.
  It is the metric the whole of Paper II is built on.
- **Neishtadt 1986**, Sov. J. Plasma Phys. **12** 568–73. Still missing; the copy
  obtained was Neishtadt–Vasiliev, *Nonlinearity* **12** (1999) 303, a different
  problem class. Bibliographic data confirmed from that paper's reference list.
- **Bishop 1972** in the original. Currently carried by Singh–Gupta's statement
  *with proof* and by Meena–Zawadzki, which agree.

## Coverage against the two manuscripts

Counting every work cited by `paper1_JMP.tex` or `paper2.tex`, excluding software
and our own Zenodo deposits: **79 cited works, 50 held here, 29 not held.**

The 29 not held are books and classical background — Cesari, Pontryagin, Liberzon,
Chandrasekhar, Whittaker–Watson, Baker, Mumford, Griffiths–Harris, Faraoni,
Bender–Orszag, Baird–Wood — together with short classical notes cited for standard
facts (Misner–Sharp 1964, Lindquist–Schwartz–Misner 1965, Martel–Poisson 2001,
Natário 2009, Hackmann 2008/2010, Booth 2005), the genus-two background
(Beilinson–Levin, Bloch, Chen, Zagier, Buchstaber–Enolskii), the two Hamilton
papers, Timofeev 1978, and the three listed under *Still missing* below. The
machine audit confirms these exist and are correctly identified; it cannot confirm
they say what we attribute to them, which is why the reading has concentrated on
the claim-bearing tail instead.

Two files here are deliberately **not** cited: Piccione–Tausk 2000, collected as
the closest free substitute for the paywalled Giannoni–Piccione–Verderesi 1997,
and Ichikawa 2023, adjacent to the genus-two cluster.
