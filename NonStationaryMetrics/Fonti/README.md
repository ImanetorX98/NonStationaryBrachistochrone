# Sources for the bibliographic audit

The PDFs in this directory are **not committed** — they are published papers
under publisher copyright (APS, IOP, AIP, Springer, Elsevier). `.gitignore`
excludes `Fonti/*.pdf`. This file is tracked and records what was collected and
what was checked, so the audit is reproducible from the DOIs alone.

**Naming convention.** Every file is `Author-Year-topic.pdf`: first author's
surname (or concatenated surnames for two or three authors, `EtAl` beyond that),
four-digit year, then a short hyphenated topic. No download identifiers, no arXiv
numbers. Preprints with no year on the title page use `undated` rather than a
guessed one. Published and preprint versions of the same work may both be held,
with suffixes such as `-arXiv`; byte-identical duplicates are removed by a
`sha256` pass (the first, on 2026-08-23, removed 44 duplicate downloads).

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
| Meena–Zawadzki 2024 | core | **✗ +** Thm 5 has *three* conditions, not one; Prop. 1 and Cor. 7 supplied two results we had been asserting. Riemannian statements: not imported automatically into our Lorentzian setting, and their navigation wind is not identified with the rail selector |
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
| Ishihara 1979 | intro + §§1–3, Thm 5.1 | **✓ +** his Thm 5.1 *is* a genuine dichotomy, but hypothesised on *preserving harmonic functions*, not on semiconformality — our garbled version had swapped the two. His §1 criterion (submersion harmonic ⟺ minimal fibres) is Riemannian: non-geodesic fibres (1d; obstruction = static-observer thrust) indicate the answer, but in our Lorentzian setting with non-constant dilation the tension is computed in full, as in the manuscript |
| Doran 2000 | full | **✓ +** our azimuthal shift is exactly his Eq. (21) composed with BL→EF, verified symbolically; his chart is built on *free-falling* observers, which is why it cannot rescue the control domain |
| Kodama 1980 | core | **✓ +** claims confirmed for our normalisation; his own K differs in *norm* (his Eq. 3.6), now noted; **he states K is timelike outside the apparent horizon** — third independent route to our control-domain boundary |
| Böhm–Wilking 2008 | thm 1 + intro | **✗ +** their Thm 1 needs a *compact* manifold and *2-positive* operator; and their intro dates the results — Hamilton **1982** is the 3d one, 1986 the 4d one, so we were citing the wrong Hamilton for a 3d base |
| Cary–Escande–Tennyson 1986 | core §§I–II | **✗ +** the jump is *deterministic at lowest order* (area of the final lobe), phase dependence only at $O(\varepsilon\ln\varepsilon)$. Supplies the mechanism — the frozen period diverges on the separatrix, so $\delta=\varepsilon T_0$ is large however small $\varepsilon$ |
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
| Hayward–Mukohyama–Ashworth 1999 | §§1–2 | **✓** their $k=*\dd r$ is our Kodama vector, "a dynamic analogue of a stationary Killing vector"; and they state plainly that "Misner & Sharp originally defined $E$", confirming that attribution |
| McVittie 1933 | §§I–II | **✓** confirms Kaloper at the source: his Eq. (1) is Schwarzschild in *isotropic* coordinates, and he builds the solution by requiring the matter to be "at rest" in his system — zero coordinate velocity and momentum, i.e. the "no accretion" property |
| Nolan 1998 | — | pending |
| Kaloper–Kleban–Martin 2010 | §§I–III | **✗ +** McVittie is **not** a conformal rescaling — its lapse carries no scale factor and its mass is constant ("there is no accretion"). Thakurta belongs to the *conformal* class instead. They also show Nolan's central assertion is **incorrect**, and call that literature "riddled with basic errors" |
| Mello–Maciel–Zanchin 2017 | abstract + §§I–II | **✓✓ ++** they analyse the **non-rotating** Thakurta metric (angular momentum set to zero, §II). Published result: for increasing unbounded $a(t)$ it "does **not** describe a cosmological black hole" — a useful spherical comparison for our reading, not a classification of the global horizons of the rotating case. Also: Thakurta is a *generalised* McVittie at the boundary case $\dot m/m=\dot a/a$, and is **not** Sultana–Dyer (factor in $\eta$ vs EF advanced time), a confusion they flag |
| Vaidya 1951 | §§1–3 | **✗** his paper is the **radiating** star — outgoing, emitting. The ingoing accreting form we use is its time-reverse; standard, and standardly called his, but not what he computes |
| Kerr 1963 | — | attribution of the seed metric; standard, no claim beyond it |
| Carter 1968 | abstract + §1 | **✓ +** exactly our claim: the fourth constant comes "from the unexpected separability of the Hamilton–Jacobi equation". Also corroborates the twist computation — the principal null congruences "have nonzero rotation (except when $a$ vanishes)" |
| Taş 2025 | abstract + §1 | **✓ +** description accurate; but a real overlap to declare — his Schwarzschild reduction *is* our arrival-time optical metric, and he reaches transverse-variation/stability too. They agree because in a static spacetime both problems reduce to the same arrival-time optical metric; the controlled rail is not the free fixed-energy geodesic |
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
| D'Hoker–Schlotterer 2024 | via [69] of DHHS | **✓** resolved without the PDF: DHHS cite it as where the higher-genus Fay identities "and their proof" are given (§5). An earlier description of them as "still being developed" was corrected |
| Baune et al. 2024 | abstract + §§5.5–7 | **✓ +** author names confirmed against the title page (the earlier Konstantin/Egor fix holds). Their test curve is genus-2 real hyperelliptic with real roots — our class. They state plainly that convergence is *argued*, no error estimates, code unreleased, prototype: strong support for our hedge, and a clean contrast with our archived theta-nome evaluation |
| Giannoni–Piccione–Verderesi 1997 | §I | **✓ ++** their $U_k=\{\beta<k^2\}$, i.e. $|Y|<k$ (their Eqs. 2–3), **is our control domain** — Lemma 2's $\hat E>|W|$, with their strict inequality matching where our slice degenerates. A thirty-year precedent for the control-domain boundary. Also a third independent source for "the two variational problems are essentially different" |
| Dyer–Honig 1979 · Sultana–Dyer 2004 | via SD2005 refs [7],[8] | **✗** the CSLS/CKH *terminology* is Dyer–Honig 1979 (JMP **20** 409) and the twist *criterion* is Sultana–Dyer **2004** (JMP **45** 4764) — SD 2005 restates both. We had credited 2005 for both. Verified on the PDF's reference list and on Crossref |
| Randers 1941 | full | **✓ +** attribution correct; and he *names* the gauge freedom $b\mapsto b+d\phi$ (his "$k$ transformation", distinguished from Weyl's) that our Doran remark relies on, and states the exactness criterion our non-rotating restriction turns on |
| Kovner 1990 | core §§1–2 | **✗✗ +** *not* a context citation — his abstract covers "arbitrary, stationary, and **nonstationary** metrics, for **massive** and massless particles". We had bundled him under the stationary case, a live priority risk. Now distinguished on four grounds: his extremals are *geodesics*, his only constraint is the mass shell, he asks for extremality not minimality, and his one genuinely non-stationary application is *perturbative lensing* — a static lens rippled by a weak gravitational wave (§§IV–V), asking where images form, not what the time-optimal constrained worldline is |
| Gibbons et al. 2009 | core §§1–2, 3.4 | **✗ +** their triality is for *null* geodesics, so our "reproduces Randers/Zermelo" needed narrowing to the form; **+** their Finsler condition $|b|_a<1$ fails exactly on the ergosurface, and the Randers data are conformal invariants — a fourth independent route to our control domain |
| Bao–Robles–Shen 2004 | §0 | **✓ +** attribution correct; their "maximal domain" $|W|_h<1$ is the Zermelo form of the same number as Gibbons' $|b|_a<1$, giving our control domain the reading *the wind outruns the ship*; and their Thm 3.1 places constant flag curvature on the homothety rung, strictly below our selector |
| Brown–Levin 2011 | abstract + §1 | **✗** theirs are the **multiple** elliptic polylogarithms — several variables, configuration space of $n{+}1$ points — which "generalize the classical elliptic polylogarithms". Our $W_{jk}$ is single-variable, length and weight two, so it belongs to the **classical** class (Beilinson–Levin, Zagier). They also build on **Chen's** reduced bar construction, confirming that citation |
| Myers 1941 | full | **✗ +** our argument was a non sequitur — he *proves* sectional decay is not enough ($S^2\times S^2$, §4). Retested on Ric: $R_{\rm opt}<0$ pointwise at every energy, so the hypothesis fails by sign. His Lemma needs no completeness |
| Filippov 1962 | full | **✗ +** bib entry was broken (DOI orphaned outside it); every Thm. 1 hypothesis re-checked and met with margin; his §III *proves* our no-relaxation claim, §V is the counterexample behind our remark, §IV explains why the domain is compact |

## Sources suggested during review

Works suggested for the introduction and the related-work discussion during the review
of Paper II. Every item is held here and has been read; bibliographic data were checked
against the title pages and Crossref rather than copied from the suggestion.

| Suggested work | Here? |
|---|---|
| Perlick 1991, JMP **32**(11) 3148–3157 | ✓ read |
| Giannoni, Piccione & Verderesi 1997, JMP **38**(12) 6367–6381 | ✓ read — suggested under the author list *Giannoni, Masiello & Piccione*; the title page gives Giannoni, Piccione, **Verderesi**, as in our bibliography |
| "The arrival time brachistochrones" | ✓ read — two papers share this theme. "The arrival time brachistochrones in general relativity" is **Giannoni & Piccione**, *J. Geom. Anal.* **12** (2002) 375 (`10.1007/BF02922047`); Giannoni–Piccione–Tausk is the *travel*-time Morse theory, DCDS **8** 697. Both are held and cited |
| Taş 2026, `arXiv:2512.08776` | ✓ read (v3) — `Tas-2025-brachistochrone-ruled-timelike-surfaces.pdf`, 37 pp. Def. 3 (ruled surface Σ(s,u)=γ_s(u)), §4 (flat: totally geodesic planes), Prop. 4 (Jacobi equation of the rulings; needs strong Legendre, so L=F²/2 for Randers), conclusions (Kerr and non-stationary open) |
| Giannoni–Piccione–Tausk 1999, `arXiv:math-ph/9905007` | ✓ read (= the DCDS paper above) |
| Haws & Kiser 1995, Amer. Math. Monthly **102** 328–336 | ✓ read |
| Hamilton 1986, JDG **24** 153–179 | ✓ read |
| Böhm & Wilking 2008, Annals **167** 1079–1097 | ✓ read |
| Myers 1941, Duke Math. J. **8**(2) 401–404 | ✓ read |
| Lecian, *Eigenvalues of the curvature of Einsteinian weighted solitons* | ✓ read |
| Lecian, *New Geodesics 2-planes solitons connected at infinity* | ✓ read |
| Lecian, *Geodesics completeness and Cauchy hypersurfaces*, Axioms **14**(12) 896 | ✓ read |
| Lecian, Axioms **15**(4) 267 (Yamabe) | ✓ read |
| Lecian, *Pseudo-Riemannian solitons after umbilicity conditions* (T2) | ✓ read |
| Lecian, *Generalized Schwarzschild solitons, spherically-symmetric weights* | ✓ read |
| Lecian, *Experimental validations of the geodesics 2-plane from GR pseudospherical cylinders*, `10.13140/RG.2.2.24140.58246` | ✓ read — `Lecian-2026-geodesics-2-plane-experimental-validations.pdf`, 13 pp, 7 Sep 2026. Distinct from *New Geodesics 2-planes solitons connected at infinity* |
| Chanda, *More on Jacobi metric: Randers–Finsler metrics, frame dragging and geometrisation techniques*, `arXiv:1911.06321`; *Eur. Phys. J. Plus* **139**(11) 2024, `10.1140/epjp/s13360-024-05775-y` | ✓ held — one work: `Chanda-2024-Randers-Finsler-frame-dragging-geometrisation.pdf` is v13 (5 Nov 2024), the published version; §5.2.1 treats the Kerr metric |
| Chanda, Gibbons, Guha, Maraner & Werner, *Jacobi-Maupertuis Randers-Finsler metric for curved spaces and the gravitational magnetoelectric effect*, *JMP* **60**(12) 122501 (2019), `10.1063/1.5098869`, `arXiv:1903.11805` | ✓ held — `ChandaEtAl-2019-Jacobi-Maupertuis-Randers-Finsler-magnetoelectric.pdf` |

The ResearchGate identifier above resolves to no Crossref record.

**What the Chanda papers contain.** Chanda–Gibbons–Guha–Maraner–Werner 2019
builds the Jacobi–Maupertuis Randers–Finsler metric for a **stationary** spacetime at
**fixed conserved** $p_0$, eq. (16):
$ds_J=\sqrt{\frac{p_0^2-(mc)^2g_{00}}{g_{00}}\gamma_{ij}dx^idx^j}+p_0\frac{g_{0i}}{g_{00}}dx^i$
with $\gamma_{ij}=-g_{ij}+g_{0i}g_{0j}/g_{00}$, valid where the Randers convexity
condition $\sqrt{a^{ij}b_ib_j}<1$ holds; eq. (24) is the explicit Kerr case. Their §1
states that the section-4 result "agree[s] with a remark in" Perlick, *The
brachistochrone problem in a stationary space-time*, JMP **32** (1991) 3148 — the
paper Paper II builds on. The agreement is with Perlick's **Proposition 3.3**, the
*free-fall* statement, not with his Proposition 3.2, which is the
$t$-brachistochrone (see below). The construction presupposes a conserved $p_0$ and a
*geodesic*; the rail of Paper II is neither free nor at conserved charge once $A$
runs, since $\dd\Ehat/\dd\tau=\varepsilon$.

**Arrival-time functional versus Jacobi–Maupertuis.** Chanda's eq. (16) reduces the
*free geodesics* at fixed energy. The arrival-time brachistochrone is a different
functional: eliminating proper time from the constraint gives
$F_T=\sqrt{\frac{E^2}{f(E^2-f)}h_{ij}dx^idx^j}-\omega_idx^i$, whose quadratic part
is exactly $\Lambda^2h$ of eq. (12) of Paper II, whereas Chanda's is
$\frac{E^2-f}{f}h$. The two differ by the position-dependent factor
$E^2/(E^2-f)^2$, so they do **not** share geodesics — in Schwarzschild at $E^2=2$
the circular extremal sits at $r=2\sqrt2M$ for the arrival cost and at
$r=(1+\sqrt5)M$ for Jacobi–Maupertuis.

This is already in **Perlick 1991, p. 3153**. His **Proposition 3.2** states that a
curve is a $t$-brachistochrone of specific energy $e^C$ iff it minimises
$\tilde S_C=\int[\sqrt{\bar h_C(\xi',\xi')}-\psi(\xi')]\,ds$ — a Randers functional,
which he compares to "a charged particle moving in a magnetostatic field".
His **Proposition 3.3** is the free-fall counterpart with a different metric
$\hat h_C$, and his comparison of Props 3.1–3.3 shows the brachistochrone problem
is equivalent to free fall only *in some other stationary spacetime*, with the
modified potential of his eq. (50). The Randers form of the brachistochrone was
therefore published in 1991.

Chanda 2024 (the EPJ Plus paper) adds, at eqs. (3.2.4)–(3.2.5): "optical metrics are
not Jacobi metrics for null curves, and their similarity for Riemannian static
metrics … is merely coincidence." The identification of the optical/Perlick structure
with a JMRF metric is thus denied by the same source in the stationary rotating case,
and granted only in the static Riemannian one.

**What the Lecian experimental-validations preprint contains.** §7, *Analytical
details of geodesics deviations … the geodesics separations*, is the textbook Jacobi
equation: the covariant separation velocity (15), its second covariant derivative
(16), and the geodesic-deviation equation
$D^2\xi^\mu/d\tau^2=R^\mu{}_{\nu\alpha\beta}u^\nu u^\alpha\xi^\beta$ (17). The
background is the pseudospherical-cylinder soliton
$ds^2=-A(r)dt^2+B(r)dr^2+r^2(d\theta^2+\sinh^2\theta\,d\varphi^2)$ of eq. (1) —
hyperbolic angular sector, static weights. §8 lists measurement techniques
(gyroscopic precession, VLBI photon-ring profiling, gravitational-wave strain).
No equation in it involves Kerr, a Randers–Finsler structure, a brachistochrone or
a control constraint: it supplies a *method* (write the Jacobi equation on the
geometry and evaluate the separation), not a result on our spacetime. §5 repeats §4
(eqs. (9)–(11) are eqs. (6)–(8)).

Cesari 1983 and Filippov 1962 entered the bibliography independently of these
suggestions; Cesari is listed under *Still missing* because it could not be obtained.

Also collected, adjacent: Piccione–Tausk, *Variational aspects of the geodesic
problem in sub-Riemannian geometry* (`math/9911215`) — the closest free substitute
for Giannoni–Piccione–Verderesi.

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

*Census as of the second revision (September 2026).* Counting every work cited by
`paper1_JMP.tex` or `paper2.tex` at that time, excluding software and our own Zenodo
deposits: **79 cited works, 50 held here, 29 not held.** The sources added in the
third revision are listed in the sections below.

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

## Further sources for Paper II

Collected after the second revision of Paper II. The PDFs are not redistributed (they stay
out of Git); each row gives the bibliographic data and what was checked in the file. Only the
passages relevant to the manuscripts were read; no observational data were re-analysed.

### Geometry and optimal control

| File | Reference | Checked |
|---|---|---|
| `BakryEmery-1985-diffusions-hypercontractives.pdf` | Séminaire de Probabilités XIX, LNM **1123** (1985) 177–206; `10.1007/BFb0075847` | Prop. 3, p. 187: for L = Δ + b on a Riemannian **or pseudo-Riemannian** manifold, Γ₂ = (Hess f\|Hess g) + (Ric − ∇^sym b)(∇f,∇g); for b = grad h this is Ric − Hess h (the m = ∞ tensor). The finite-m term −(1/m) dφ⊗dφ is not in this paper (it is in Case 2012). The scanned text layer is poor: read the page images. Cited in "A Bakry–Émery identity" |
| `Case-2012-smooth-metric-measure-spaces-quasi-Einstein.pdf` | Internat. J. Math. **23** (2012) 1250110; `10.1142/S0129167X12501108`; arXiv:1011.2723 | The m-Bakry–Émery tensor Ric + ∇²φ − (1/m) dφ⊗dφ (p. 1); quasi-Einstein metrics include conformally Einstein metrics for negative m (p. 2) |
| `Case-2013-sharp-metric-obstructions-quasi-Einstein.pdf` | J. Geom. Phys. **64** (2013) 12–30; `10.1016/j.geomphys.2012.10.006`; arXiv:1110.3010 | p. 2: m = 2 − n corresponds to an Einstein metric in a conformal class. Cited for this value of m |
| `GiannoniMasielloPiccione-1998-Morse-theory-light-rays-stably-causal.pdf` | Ann. IHP Phys. Théor. **69** (1998) 359–412; Numdam `AIHPA_1998__69_4_359_0` | §1–3: Morse relations for light rays from a point to a timelike curve in light-convex subsets of stably causal manifolds (Thm 1.7); Kerr outside the stationary limit (§3.4) |
| `GraichenPetitKugi-2008-state-constraint-transformation.pdf` | Proc. 47th IEEE CDC (2008); `10.1109/CDC.2008.4739035` | Saturation of the state constraint and regularisation of a fixed-final-time problem; interior conditions are junctions between arcs. No normality theorem for a free-final-time problem such as ours |
| `GraichenPetitKugi-2008-HAL-hal-00483077v1.pdf` | HAL hal-00483077v1 (same paper) | Text identical to the proceedings copy (cover and footers differ). Contains no statement on normality ("normality", "λ₀", "abnormal" and "inactive" do not occur) |

### Observational context

| File | Reference | Checked |
|---|---|---|
| `SahuEtAl-2022-isolated-black-hole-astrometric-microlensing.pdf` | ApJ **933**, 83 (2022); `10.3847/1538-4357/ac739e` | Astrometry, photometry and parallax separate mass, distance and proper motion of a dark lens; 7.1 ± 1.3 M☉ is the paper's result |
| `CiurloEtAl-2020-dust-enshrouded-objects-Galactic-black-hole.pdf` | Nature **577**, 337–340 (2020); `10.1038/s41586-019-1883-y` | Six G objects; gas and dust; Keplerian orbits reconstructed from astrometry and radial velocities |
| `MastersonEtAl-2025-JWST-tidal-disruption-events.pdf` | ApJL **988**, L48 (2025); `10.3847/2041-8213/ade153` | Four TDEs; high-ionisation lines from accretion and infrared dust emission. Not evidence of absent accretion during the flare |
| `GaiaCollab-2024-Gaia-BH3-33-solar-mass-arXiv.pdf` | A&A **686**, L2 (2024); `10.1051/0004-6361/202449763` | Gaia BH3: 32.70 ± 0.82 M☉; very metal-poor halo giant companion ([Fe/H] ≈ −2.5); isochrones 12–14 Gyr |
| `ElBadryEtAl-2023-Gaia-BH1-arXiv.pdf` | MNRAS **518**, 1057 (2023); `10.1093/mnras/stac3140` | Gaia BH1: 9.62 ± 0.18 M☉; companion age ≳ 1 Gyr (gyrochronology) |
| `ElBadryEtAl-2023-Gaia-BH2-arXiv.pdf` | MNRAS **521**, 4323 (2023); `10.1093/mnras/stad799` | Gaia BH2: 8.94 ± 0.34 M☉; age 5–13 Gyr (8–12 Gyr if thick disk) |
| `GRAVITY-2022-mass-distribution-Galactic-Centre-arXiv.pdf` | A&A **657**, L12 (2022); `10.1051/0004-6361/202142465`; arXiv:2112.07478 | Tab. 4: mass and distance of Sgr A*, statistical errors |
| `EHT-2019-M87-shadow-mass-VI-arXiv.pdf` | ApJL **875**, L6 (2019); `10.3847/2041-8213/ab1141`; arXiv:1906.11243 | Mass of M87* from imaging and GRMHD calibration; statistical and systematic errors; no model-independent spin |
| `EHT-2022-SgrA-astrophysical-models-V-arXiv.pdf` | ApJL **930**, L16 (2022); `10.3847/2041-8213/ac6672`; arXiv:2311.09478 | Kerr-plus-plasma model comparison; spin constraints depend on the astrophysical assumptions |
| `MillerJonesEtAl-2021-CygnusX1-mass-distance-arXiv.pdf` | Science **371**, 1046–1049 (2021); `10.1126/science.abb3363`; arXiv:2102.09091 | Mass and distance of Cygnus X-1, an accreting system |
| `ZhaoEtAl-2021-CygnusX1-spin-arXiv.pdf` | ApJ **908**, 117 (2021); `10.3847/1538-4357/abbcd6`; arXiv:2102.09093 | Spin bound from continuum fitting, under disk and alignment assumptions |
| `Planck-2020-2018-cosmological-parameters-VI-arXiv.pdf` | A&A **641**, A6 (2020); `10.1051/0004-6361/201833910`; arXiv:1807.06209 | Fiducial H₀ and Ωm in flat ΛCDM; does not show that the paper's A is the local cosmological factor |
| `PrsaEtAl-2016-IAU-nominal-solar-quantities-arXiv.pdf` | AJ **152**, 41 (2016); `10.3847/0004-6256/152/2/41`; arXiv:1605.09788 | Tab. 1: nominal solar GM, the conversion constant for geometric scales |
| `PitjevaEtAl-2021-solar-GM-change-rate-EPM2019-AA.pdf` | A&A **647**, A141 (2021); `10.1051/0004-6361/202039893`; arXiv:2201.09804 | §3, eq. (1): d(GM☉)/dt / GM☉ = (−10.2 ± 1.4)×10⁻¹⁴ yr⁻¹ (3σ), from EPM2019 (over 850 000 observations, 1913–2017); eq. (12): −13.4 < Ṁ☉/M☉ < −8.7 (×10⁻¹⁴ yr⁻¹); eq. (14): −2.9 < Ġ/G < +4.6 (×10⁻¹⁴ yr⁻¹, 3σ). Used to exclude a solar gravitational parameter growing at H₀ ≈ 6.89×10⁻¹¹ yr⁻¹ |

### The Thakurta metric and cosmological coupling

| File | Reference | Checked |
|---|---|---|
| `HutsiEtAl-2021-Thakurta-not-cosmological-black-holes-arXiv.pdf` | EPJC **81**, 999 (2021); `10.1140/epjc/s10052-021-09803-4` | The conditions producing Thakurta describe a single black hole **accreting** in a dust universe; the required energy flux is not supplied by radiation, PBHs or dark matter |
| `BoehmEtAl-2021-comment-Thakurta-arXiv.pdf` | arXiv:2105.14908 (comment) | Defence of Thakurta: Misner–Sharp mass, no assumption of matter accretion |
| `HutsiEtAl-2021-reply-comment-Thakurta-arXiv.pdf` | arXiv:2106.02007 (reply) | The mass growth is an ingoing energy flux, that is accretion |
| `BoehmEtAl-2021-eliminating-LIGO-bounds-PBH-arXiv.pdf` | JCAP **03** (2021) 078; `10.1088/1475-7516/2021/03/078` (arXiv v3) | Thakurta proposal for PBHs, mass ∝ a(t) |
| `BoehmEtAl-2021-eliminating-LIGO-bounds-PBH-JCAP.pdf` | as above, published version | §2.1–2.3: m is today's physical mass, Misner–Sharp mass m a(t) = m/(1+z); "its increase is not due to any accretion of matter"; for R ≪ 1/H the metric tends to static Schwarzschild. Eq. (2.3) also contains a cosmological term in the total Misner–Sharp mass; m a is the central term |
| `HaradaMaedaSato-2022-Thakurta-not-cosmological-black-hole-arXiv.pdf` | PLB **833**, 137332 (2022); `10.1016/j.physletb.2022.137332` | Spherical Thakurta with power-law conformal factor: no event or trapping horizon of the cosmological-black-hole type; violation of all standard energy conditions in specific regions. Not a global classification for arbitrary A or non-zero spin |
| `KobakhidzePicker-2022-apparent-horizons-Thakurta-arXiv.pdf` | EPJC **82**, 347 (2022); `10.1140/epjc/s10052-022-10312-1` | Opposite view: with the Kodama foliation the horizon is future outer trapping |
| `MacielZanchin-2024-comment-apparent-horizons-Thakurta-arXiv.pdf` | EPJC **84**, 1109 (2024); `10.1140/epjc/s10052-024-13489-9` | Kobakhidze–Picker use coordinates covering a region beyond the singularity |
| `PatelKoppar-1988-Kerr-like-radiating-metric-expanding-universe-ActaPhysHung.pdf` | Acta Phys. Hung. **64**, 353–360 (1988); `10.1007/BF03053840` | Radiating Kerr–NUT (Vaidya–Patel–Bhatt) times e^{2F(t)}, with t paired with the null coordinate u (not BL time): perfect fluid + null radiation + Λ; constant m; for a = 1, N = 0 "Kerr metric in the background of an expanding universe". §1 records that Guha Thakurta (Indian J. Phys. 55B, 304, 1981) uses a "perfect fluid with heat flow". Energy conditions not discussed |
| `FarrahEtAl-2023-cosmological-coupling-black-holes-arXiv.pdf` | ApJL **944**, L31 (2023); `10.3847/2041-8213/acb704` | Cosmological coupling M ∝ aᵏ, k ≈ 3 from SMBHs in elliptical galaxies |
| `MisteleEtAl-2023-comment-Farrah-cosmological-coupling-arXiv.pdf` | RNAAS **7**, 101 (2023); `10.3847/2515-5172/acd767` | Critique of Farrah et al. |
| `AndraeElBadry-2023-Gaia-cosmological-coupling-arXiv.pdf` | A&A **673**, L10 (2023); `10.1051/0004-6361/202346350` | Gaia BH1/BH2: with k = 3 the birth mass of BH2 would be < 2.2 M☉ (77%) |
| `LacyEtAl-2024-coupling-accretion-history-SMBH-arXiv.pdf` | ApJL **961**, L33 (2024); `10.3847/2041-8213/ad1b5f` | Constraints from the accretion history of SMBHs |
| `Rodriguez-2023-coupling-NGC3201-arXiv.pdf` | ApJL **947**, L12 (2023); `10.3847/2041-8213/acc9b6` | Constraints from black holes in the globular cluster NGC 3201 |
| `CrokerWeiner-2019-symmetry-pressure-Friedmann-I-arXiv.pdf` | ApJ **882**, 19 (2019); `10.3847/1538-4357/ab32da` | Formalism of cosmological coupling |
| `CrokerEtAl-2020-symmetry-pressure-Friedmann-II-arXiv.pdf` | ApJ **889**, 115 (2020); `10.3847/1538-4357/ab5aff` | Stellar remnants and coupling |
| `CrokerEtAl-2021-coupled-compact-objects-LIGO-arXiv.pdf` | ApJL **921**, L22 (2021); `10.3847/2041-8213/ac2fad` | Coupled compact objects, one-parameter model of LIGO–Virgo masses |
| `CalzaEtAl-2024-coupled-black-holes-PTA-arXiv.pdf` | Sci. Rep. **14**, 31296 (2024); `10.1038/s41598-024-82661-8` | Consequences for pulsar timing arrays |
| `CadoniEtAl-2023-coupling-nonsingular-black-holes-arXiv.pdf` | JCAP **11** (2023) 007; `10.1088/1475-7516/2023/11/007` | Metric a²(η)[−e^α dt² + e^β dr² + r²dΩ²]: k = 1 for spherically symmetric objects; JWST compatible with k = 1 at 1σ, ellipticals at z ≈ 0.8–0.9 prefer k > 1. Eq. (2.5) keeps an internal time dependence, so it is not the same stationary conformal seed as Paper II |
| `CadoniEtAl-2026-coupled-black-holes-regular-horizons-arXiv.pdf` | arXiv:2601.03296v2 | Coupled black holes with regular horizons |
| `WuEtAl-2026-coupled-black-holes-dark-sector-arXiv.pdf` | EPJC **86**, 288 (2026); `10.1140/epjc/s10052-026-15533-2` | Exact solution with M ∝ aᵏ in an anisotropic dark-sector background |
| `FaraoniJacques-2007-cosmological-expansion-local-physics-arXiv.pdf` | PRD **76**, 063510 (2007); `10.1103/PhysRevD.76.063510` | Generalised McVittie with arbitrary m(t): imperfect fluid with heat flux; cosmic expansion and local physics |
| `BenAchourCisternaHassaine-2025-dynamical-axisymmetric-compact-objects-arXiv.pdf` | JCAP **06** (2026) 061; `10.1088/1475-7516/2026/06/061`; arXiv:2512.19542 | First exact dynamical axisymmetric compact object in FLRW (Einstein–scalar); cites Thakurta as a conformal rescaling with energy-condition violations |
| `BenAchourCisternaHassaine-2026-dynamical-axisymmetric-compact-objects-JCAP.pdf` | as above, published version | As the arXiv version |
| `BabichevCharmousisLecoeur-2023-rotating-BH-cosmological-scalar-tensor-arXiv.pdf` | arXiv:2305.17129 | Rotating black holes in a cosmological background in scalar–tensor theories |

### Compact binaries used as consistency references

| File | Reference | Checked |
|---|---|---|
| `AntoniadisEtAl-2013-massive-pulsar-compact-binary-Science.pdf` | Science **340**, 1233232 (2013); `10.1126/science.1233232` | PSR J0348+0432: M = 2.01 ± 0.04 M☉, white dwarf 0.172 M☉, P_b = 2.46 h, white-dwarf cooling age ~2 Gyr; Ṗ_b = −8.6 ± 1.4 µs/yr, consistent with GR |
| `CromartieEtAl-2020-relativistic-Shapiro-massive-pulsar-NatAstron.pdf` | Nature Astronomy **4**, 72 (2020); `10.1038/s41550-019-0880-2` | PSR J0740+6620: 2.14 (+0.10/−0.09) M☉ from Shapiro delay; 4.77-day orbit |
| `KramerEtAl-2021-double-pulsar-strong-field-arXiv.pdf` | PRX **11**, 041050 (2021); `10.1103/PhysRevX.11.041050` | Double pulsar J0737−3039: orbital decay consistent with GR to 0.013% |
| `JacobyEtAl-2006-B2127+11C-M15-double-NS-arXiv.pdf` | ApJ **644**, L113 (2006); `10.1086/505742` | Double neutron star in M15: 1.358 and 1.354 M☉; Ṗ_b = (−3.95 ± 0.13)×10⁻¹², consistent with GR to 3% |

### Quasi-local energy and energy conditions

| File | Reference | Checked |
|---|---|---|
| `PrainEtAl-2015-Hawking-Hayward-conformal-arXiv.pdf` | CQG **33**, 145008 (2016); `10.1088/0264-9381/33/14/145008` | Hawking–Hayward quasi-local energy under conformal rescaling, without spherical symmetry; Misner–Sharp as a special case |
| `Hammad-2016-conformal-mapping-quasilocal-masses-HH-arXiv.pdf` | CQG **33**, 235016 (2016); `10.1088/0264-9381/33/23/235016` | The geometric HH mass transforms with the "wrong" conformal factor; correction analogous to Misner–Sharp |
| `Hammad-2016-conformal-mapping-quasilocal-masses-HH-CQG.pdf` | as above, published version | Distinction between geometric mass and material interpretation |
| `JarabaCamposBaldiotti-2026-quasilocal-Hamiltonian-Kerr-Schild-arXiv.pdf` | arXiv:2608.16534 | Misner–Sharp generalised to rotating generalised Kerr–Schild geometries (Iyer–Wald formalism) |
| `MartinMorunoVisser-2018-essential-core-Hawking-Ellis-types-arXiv.pdf` | CQG **35**, 125003 (2018); `10.1088/1361-6382/aac147` | Hawking–Ellis classification of stress–energy tensors; type IV (no causal eigenvector) violates all standard energy conditions |

### Related work on relativistic brachistochrones

| File | Reference | Checked |
|---|---|---|
| `WangJia-2026-brachistochrone-constant-velocity-static-spherical-arXiv.pdf` | L. Wang, J. Jia, arXiv:2607.29265v1 [gr-qc] (31 Jul 2026), 9 pp. | Brachistochrone for a traveller at **constant local speed** (not at maintained energy) in static spherically symmetric spacetimes: coordinate time ∫√(B dr² + C dθ²)/√A, i.e. null geodesics of the optical metric, time scaled by 1/v (§V, generalised Fermat principle). Closed form inside a truncated singular isothermal sphere (eq. 22), numerical in Schwarzschild (eqs. 28–30) and in a relativistic Plummer model. No rotation, no conformal factor, no thrust or propellant analysis; axisymmetric and non-stationary cases named as future work (§VI). Cites Kamath, JMP 29, 2268 (1988) as the first curved-spacetime brachistochrone (weak field) |
| `Kamath-1988-brachistochrone-almost-flat-space-JMP.pdf` | S. G. Kamath, J. Math. Phys. **29**, 2268–2272 (1988); `10.1063/1.528158`; publisher version | First general-relativistic extension of the Goldstein–Bender relativistic brachistochrone: linearised (isotropic) weak field of a static spherical source, particle **falling from rest at infinity** (the marginal energy), coordinate time T = ∫ dl/v. The curve θ(u), u = 1/r, is a hyperelliptic integral (eq. 10′–11); exact in elliptic integrals for one value of the Snell constant (§IV); for the Sun reduced numerically to r^{1/2}θ = const (eq. 27). Full Schwarzschild case only commented (§VI). No thrust, propellant or experiment |
| `AgmonMann-2026-general-features-brachistochrones-AJP.pdf` | D. Agmon, A. Mann, Am. J. Phys. **94**, 439–444 (2026); `10.1119/5.0303237`; CC BY-NC | Newtonian brachistochrone for any conservative force: necessary and sufficient condition F_n = mv²/R_c; **"Rule of 2"**, the guiding force is N = −2F_n (eq. 7), and the "Mirror Rule" (eq. 8), both already in Routh (1898). For a central force T₀ = r sin α / v is conserved (eqs. 21–26), and the guide can be replaced by a magnetic field B ∝ F(r)/r (eq. 29). Relativistic laboratory-time version N/F_n = 2 − β² announced (§II), not derived here |
| `GoldsteinBender-1986-relativistic-brachistochrone-JMP.pdf` | H. F. Goldstein, C. M. Bender, J. Math. Phys. **27**, 507–511 (1986); `10.1063/1.527199`; publisher version | Special-relativistic brachistochrone of a particle starting at rest in a **uniform** electric or gravitational field (flat spacetime, velocity-dependent "gravitational" force mg/√(1−v²/c²)). General solution y = ∫ k v/√(1−k²v²) dx (eq. 6); three regimes k² < 1, = 1, > 1 for laboratory time (§III A, C); separate **proper-time** brachistochrones (§III B, D), which lose the unbounded solutions. The earliest explicit distinction between laboratory-time and proper-time brachistochrones |

### Orbital mechanics and mission references for the Solar-System comparison

| File | Reference | Checked |
|---|---|---|
| `Battin-1999-mathematics-methods-astrodynamics-AIAA.pdf` | R. H. Battin, *An Introduction to the Mathematics and Methods of Astrodynamics*, Revised Edition, AIAA (1999); `10.2514/4.861543` | §9.1, eqs. (9.1)–(9.3): turn angle of a hyperbolic approach, sin ν = 1/(1 + r_m v∞²/μ), total deflection 2ν; §9.2, eqs. (9.9)–(9.10): flyby turn angle from inbound and outbound relative velocities and the corresponding passing distance; §8.5: Laplace sphere of influence. The pericentre formula of `probe_saturn_comparison.py` is eq. (9.10) |
| `JacobsonEtAl-2006-Saturn-system-gravity-field-AJ.pdf` | R. A. Jacobson et al., AJ **132**, 2520–2526 (2006); `10.1086/508812` | Final solution: GM of Saturn 37 931 207.7 ± 1.1 km³ s⁻² (value used in the scripts) |
| `ArchinalEtAl-2018-IAU-WGCCRE-2015-rotational-elements-CMDA.pdf` | B. A. Archinal et al., Celest. Mech. Dyn. Astron. **130**, 22 (2018); `10.1007/s10569-017-9805-5` | Saturn pole α₀ = 40.589 − 0.036T, δ₀ = 83.537 − 0.004T (Table 1), equatorial radius 60 268 ± 4 km (1-bar level); both used in the scripts |
| `GuoFarquhar-2008-New-Horizons-mission-design-SSR.pdf` | Y. Guo, R. W. Farquhar, Space Sci. Rev. **140**, 49–74 (2008); `10.1007/s11214-007-9242-y` | Design maximum launch energy C3 = 164 km² s⁻²; achieved C3 = 157.75 km² s⁻² (Table 7), Atlas V 551 with Star 48B upper stage, Jupiter gravity assist |
| `Spilker-2019-Cassini-Huygens-Saturn-system-13-years-Science.pdf` | L. Spilker, Science **364**, 1046–1051 (2019); `10.1126/science.aat3760` | Thirteen years of Cassini in the Saturn system; Grand Finale: 22 orbits through the gap between the rings and the planet (2017), a flight precedent for close passages and ring-plane crossings |
| `PattersonBenson-2008-NEXT-ion-propulsion-status-performance-NASA-TM.pdf` | M. J. Patterson, S. W. Benson, NASA/TM-2008-214986, AIAA 2007-5199; `10.2514/6.2007-5199` | NEXT xenon ion thruster: 0.54–6.9 kW input power, maximum thrust > 236 mN, maximum specific impulse > 4170 s (exhaust speed ≈ 41 km/s), peak efficiency > 70%, xenon throughput requirement > 300 kg |
| `BensonRiehlOleson-2007-NEXT-ion-propulsion-Saturn-system-exploration-AIAA.pdf` | S. W. Benson, J. P. Riehl, S. R. Oleson, AIAA 2007-5230 (NTRS 20080006604); `10.2514/6.2007-5230` | NEXT systems of two to four thrusters for Saturn missions, used in the inner Solar System with a Venus or Earth gravity assist and a solar-electric module (17 kW array); 2400–4000 kg delivered to Saturn |
| `GibsonEtAl-2018-KRUSTY-Kilopower-ground-test-results-AIAA.pdf` | M. A. Gibson et al., AIAA 2018-4973 (NTRS 20180005435); `10.2514/6.2018-4973` | Kilopower fission power, 1–10 kWe class; KRUSTY ground test of the 1 kWe flight-concept core at about 4 kWt and 800 °C |
| `ThorntonBorder-2003-radiometric-tracking-deep-space-navigation-JPL.pdf` | C. L. Thornton, J. S. Border, *Radiometric Tracking Techniques for Deep-Space Navigation*, JPL DESCANSO Monograph 1 (Wiley, 2003); `10.1002/0471728454` | X-band DSN performance: Doppler random error 0.03 mm/s and range random error 60 cm for 60-s averages, range bias 2 m; ΔDOR angles at the 5–50 nrad level |

### Space clocks and time transfer (solar-sail clock comparison, Paper II `par:obs-sail`)

Collected on 2026-10-22 to anchor the clock requirement of the solar-sail proposal (fractional stability of about
4×10⁻¹³ over 72 h, or 1.5×10⁻¹² over 25 h, for a 10% measurement of the response of X = t − τ; about 3×10⁻¹⁵ over a
day for the ~1 ns inversions without daughters). Delva2018, Burt2021 and Samain2014 are cited in `par:obs-sail`; ACES is not cited.

| File | Reference | Checked |
|---|---|---|
| `BurtEtAl-2021-trapped-ion-atomic-clock-space-Nature.pdf` (published version, 19 pp.); also `BurtEtAl-2021-trapped-ion-atomic-clock-space-preprint.pdf` (Research Square `10.21203/rs.3.rs-122117/v1`, CC BY 4.0) | E. A. Burt et al., Nature **595**, 43–47 (2021); `10.1038/s41586-021-03571-7` | NASA Deep Space Atomic Clock (Hg⁺ trapped ion), launched June 2019, powered on August 2019: in space 3–5×10⁻¹⁵ at one day (depending on temperature variations) and 3×10⁻¹⁵ at 23 days, time deviation < 4 ns at 23 days, uncorrected drift +3.0(0.7)×10⁻¹⁶/day; short-term 1–2×10⁻¹³ at 1 s on the ground; < 50 W. States that other space clocks have short-term stability 1–10×10⁻¹²/√τ and drift 1–10×10⁻¹⁵/day. The published version confirms the one-day range, 3×10⁻¹⁵ to 5×10⁻¹⁵ «depending on the temperature variations for a given run», and the 23-day value with time deviation < 4 ns. No mass figure in either version |
| `DelvaEtAl-2018-gravitational-redshift-eccentric-Galileo-arXiv.pdf` | P. Delva et al., Phys. Rev. Lett. **121**, 231101 (2018); `10.1103/PhysRevLett.121.231101`; arXiv:1812.03711 | Galileo 5 and 6 carry two passive hydrogen masers (PHM) and two rubidium clocks (RAFS) each. Table I: standard deviation of clock-bias residuals after a daily linear fit, 0.11–0.20 ns for the PHMs used and 0.69 ns for one PHM, 2.21 ns for a RAFS. Redshift deviation (0.19 ± 2.48)×10⁻⁵ |
| `CacciapuotiEtAl-2024-Atomic-Clock-Ensemble-in-Space-arXiv.pdf` | L. Cacciapuoti et al. (28 authors), «Atomic Clock Ensemble in Space», arXiv:2411.02912 [gr-qc] (5 Nov 2024) | Pre-launch status of ACES. PHARAO (laser-cooled Cs) specified to 1×10⁻¹³/√τ for 1 s < τ < 10⁶ s and accuracy 1–2×10⁻¹⁶; SHM active space hydrogen maser of 46 kg, 70 l, 93 W, stability close to ground active masers; clock signal at 1×10⁻¹⁶. «Will be launched … on a Space X Dragon». Launch on 21 April 2025 (Falcon 9) and installation on the ISS on 28 April 2025 from the DLR blog «Déjà vu in space» (dlr.de, 2025), with commissioning then in progress; no in-orbit performance read |
| `SamainEtAl-2014-T2L2-characterization-calibration-flight-instrument-Metrologia.pdf` | E. Samain et al., Metrologia **51**, 503–515 (2014); `10.1088/0026-1394/51/5/503`; file from the DLR institutional repository (elib 91979) | T2L2 on Jason-2 (launched June 2008): laser time transfer specified to time stability better than 1 ps over 1000 s and accuracy better than 100 ps; geometric correction uncertainty about 3 ps |
