# Sources for the bibliographic audit

The PDFs in this directory are **not committed** — they are published papers
under publisher copyright (APS, IOP, AIP, Springer, Elsevier). `.gitignore`
excludes `fonti/*.pdf`. This file is tracked and records what was collected and
what was checked, so the audit is reproducible from the DOIs alone.

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
| Doran 2000 | full | **✓ +** our azimuthal shift is exactly his Eq. (21) composed with BL→EF, verified symbolically; his chart is built on *free-falling* observers, which is why it cannot rescue the control domain |
| Kodama 1980 | core | **✓ +** claims confirmed for our normalisation; his own K differs in *norm* (his Eq. 3.6), now noted; **he states K is timelike outside the apparent horizon** — third independent route to our control-domain boundary |
| Böhm–Wilking 2008 | thm 1 + intro | **✗ +** their Thm 1 needs a *compact* manifold and *2-positive* operator; and their intro dates the results — Hamilton **1982** is the 3d one, 1986 the 4d one, so we were citing the wrong Hamilton for a 3d base |
| Hayward 1996 | — | pending |
| Abreu–Visser 2010 | — | pending |
| Ashtekar–Krishnan 2003 | — | pending |
| Hayward–Mukohyama–Ashworth 1999 | — | pending |
| Randers 1941 | — | pending (context citation) |
| Kovner 1990 | — | pending (context citation) |
| Gibbons et al. 2009 | core §§1–2, 3.4 | **✗ +** their triality is for *null* geodesics, so our "reproduces Randers/Zermelo" needed narrowing to the form; **+** their Finsler condition $|b|_a<1$ fails exactly on the ergosurface, and the Randers data are conformal invariants — a fourth independent route to our control domain |
| Bao–Robles–Shen 2004 | — | pending (context citation) |
| Brown–Levin 2011 | — | pending (content claim: elliptic dilogarithm) |
| Myers 1941 | full | **✗ +** our argument was a non sequitur — he *proves* sectional decay is not enough ($S^2\times S^2$, §4). Retested on Ric: $R_{\rm opt}<0$ pointwise at every energy, so the hypothesis fails by sign. His Lemma needs no completeness |
| Filippov 1962 | full | **✗ +** bib entry was broken (DOI orphaned outside it); every Thm. 1 hypothesis re-checked and met with margin; his §III *proves* our no-relaxation claim, §V is the counterexample behind our remark, §IV explains why the domain is compact |
| unidentified, 31 pp | — | pending (arrived with group 2) |

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
| **Haws–Kiser 1995**, Amer. Math. Monthly 102, 328 | ✗ paywalled |
| **Hamilton 1986**, four-manifolds with positive curvature operator | ✗ paywalled |
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

## The ~60 references not yet examined

Mostly classical background (Kerr, Carter, Chandrasekhar, Baker, Mumford, Fay,
DLMF, Whittaker–Watson) cited for standard facts. The machine audit confirms
they exist and are correctly identified; it cannot confirm they say what we
attribute to them.
