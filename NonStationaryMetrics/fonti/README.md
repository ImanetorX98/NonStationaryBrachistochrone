# Sources for the bibliographic audit

The PDFs in this directory are **not committed** — they are published papers
under publisher copyright (APS, IOP, AIP, Springer, Elsevier). `.gitignore`
excludes `fonti/*.pdf`. This file is tracked and records what was collected and
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
| Misner–Sharp 1964 | — | pending |
| Lindquist–Schwartz–Misner 1965 | — | pending |
| Booth 2005 | — | pending |
| Hayward 1996 | — | pending |
| Abreu–Visser 2010 | — | pending |
| Ashtekar–Krishnan 2003 | — | pending |
| Hayward–Mukohyama–Ashworth 1999 | — | pending |
| McVittie 1933 | — | pending |
| Nolan 1998 | — | pending |
| Kaloper–Kleban–Martin 2010 | — | pending |
| Mello–Maciel–Zanchin 2017 | — | pending |
| Vaidya 1951 | — | pending |
| Kerr 1963 | — | pending |
| Carter 1968 | — | pending |
| Taş 2025 | — | pending |
| Giannoni–Piccione–Tausk 1999 | — | pending |
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
| Randers 1941 | full | **✓ +** attribution correct; and he *names* the gauge freedom $b\mapsto b+d\phi$ (his "$k$ transformation", distinguished from Weyl's) that our Doran remark relies on, and states the exactness criterion our non-rotating restriction turns on |
| Kovner 1990 | core §§1–2 | **✗✗ +** *not* a context citation — his abstract covers "arbitrary, stationary, and **nonstationary** metrics, for **massive** and massless particles". We had bundled him under the stationary case, a live priority risk. Now distinguished on three grounds: his extremals are *geodesics*, his only constraint is the mass shell, and he asks for extremality not minimality |
| Gibbons et al. 2009 | core §§1–2, 3.4 | **✗ +** their triality is for *null* geodesics, so our "reproduces Randers/Zermelo" needed narrowing to the form; **+** their Finsler condition $|b|_a<1$ fails exactly on the ergosurface, and the Randers data are conformal invariants — a fourth independent route to our control domain |
| Bao–Robles–Shen 2004 | §0 | **✓ +** attribution correct; their "maximal domain" $|W|_h<1$ is the Zermelo form of the same number as Gibbons' $|b|_a<1$, giving our control domain the reading *the wind outruns the ship*; and their Thm 3.1 places constant flag curvature on the homothety rung, strictly below our selector |
| Brown–Levin 2011 | — | pending (content claim: elliptic dilogarithm) |
| Myers 1941 | full | **✗ +** our argument was a non sequitur — he *proves* sectional decay is not enough ($S^2\times S^2$, §4). Retested on Ric: $R_{\rm opt}<0$ pointwise at every energy, so the hypothesis fails by sign. His Lemma needs no completeness |
| Filippov 1962 | full | **✗ +** bib entry was broken (DOI orphaned outside it); every Thm. 1 hypothesis re-checked and met with margin; his §III *proves* our no-relaxation claim, §V is the counterexample behind our remark, §IV explains why the domain is compact |

## Referee-proposed reading

Referee 2's list, and the two books referee 1 leans on, tracked separately since
the reports are answered against them.

| Proposed | Here? |
|---|---|
| Perlick 1991 | ✓ read |
| Giannoni–Piccione–Tausk 1999, Morse for travel-time brachistochrones (`math-ph/9905007`) | ✓ arXiv |
| Taş, brachistochrone-ruled timelike surfaces (`2512.08776`) | ✓ arXiv — note the copy is **v3, Apr 2026**; our bib cites the Dec 2025 v1 |
| Böhm–Wilking 2008, positive curvature operators (`math/0606187`) | ✓ **read** |
| Myers 1941 | ✓ **read** |
| Filippov 1962 | ✓ **read** |
| Lecian ×5 | ✓ |
| **Giannoni–Piccione–Verderesi 1997**, sub-Riemannian (JMP 38, 6367) | ✗ paywalled |
| **Giannoni–Piccione–Tausk 2002**, arrival-time (DCDS 8, 697) | ✗ paywalled |
| **Haws–Kiser 1995**, Amer. Math. Monthly 102, 328 | ✓ **read** |
| **Hamilton 1986**, four-manifolds with positive curvature operator | ✗ paywalled — but the statement we needed is quoted in Böhm–Wilking's introduction, which we read |
| **Cesari 1983** (referee 1's existence machinery) | ✗ book |

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
