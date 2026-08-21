# Paper I — JMP variant: state of play

Handoff note. Everything below is committed on `main` and pushed to
`nsb` (`github.com/ImanetorX98/NonStationaryBrachistochrone`).

---

## 1. What this is

`paper1/paper1_JMP.tex` — 52 pp, a variant of Paper I targeted at the
*Journal of Mathematical Physics* (AIP). Physics content unchanged from
`paper1.tex`; what differs is that the formal core is **proved** rather than
asserted, and six defects found during the audit are fixed.

Compiles clean: 0 overfull boxes, 0 unresolved references, 0 "float too large".

Companion: `paper2` (Thakurta–Kerr) is under review at **CQG**. It is a
distinct manuscript — shares authorship, not results or figures.

## 2. Numbering convention — do not break this

- **I.1–I.6** are the six shared results the companion paper cites by number.
  They are pinned with `\Ithm{k}` before each environment, because they do not
  appear in document order (I.4 and I.6 live in the appendices, I.5 in §4).
  **Never let these auto-number.**
- **I.A–I.J** are auxiliary results added by this variant. They share a
  separate counter (`auxlemma`/`auxtheorem`, letters). They renumber freely
  when new ones are inserted, which is fine — everything is `\ref`'d.

Current letters, in document order: I.A rail indicatrix · I.B time-dependent
optical metric · I.C no conjugate point on monotone arcs · I.D caustic of the
apsidal family · I.E the sweep never reaches π · I.F determinant = resultant ·
I.G scaling of the reduction · I.H third-kind letter and Hermite reduction ·
I.I residues of the source antiderivative · I.J block form.

## 3. What is proved

Formal core (§2), FLRW (§3), Vaidya rail (§4.1), turning-point laws (§4.3),
closed forms (§4.4), obstructions (§4.5), appendices A–C. Highlights:

- **Perlick extension** (I.B). The rail arrival-time functional is
  `F = ω_a ν^a + sqrt(a_ab ν^a ν^b)` with `a_ab = Ê² h_ab / (f(Ê²−f))`, the
  data free to depend on the adapted time. The computation is *pointwise in
  the event*, so stationarity enters only by removing the time argument.
  Perlick's fixed-energy metric is the stationary member; `h_ab/f` the
  `Ê→∞` limit. Hypothesis (H2) is exactly the positive-definiteness of that
  metric, and its failure is the freezing surface.
- **Determinant identity** (I.F): `det M = 32 Res_r(S,S') = −32 s₆ disc_r S`
  for an **arbitrary** sextic — not a Vaidya fact.
- **Closed forms**: all `c_k^m` and `A^m` as rational functions of
  `(Ê², J, m)` over the reduced discriminant `Δ`; substituting them back
  annihilates the reduction identity identically.
- **§4.5**: no conjugate point on monotone arcs (I.C); through periapsis a
  caustic *does* exist, proved by compactness (I.D); and the swept azimuth
  never reaches π (I.E), so the caustic is the **only** single-excursion
  obstruction and the Maxwell mechanism needs more than one excursion.
- **§3.1**: FLRW extended to arbitrary spatial curvature; `|W| = a` for every
  k, so the brachistochrone is the spatial geodesic and curvature enters only
  through the geodesic distance. Closed slices give an explicit Maxwell set at
  the antipodal cut locus.

## 4. Defects found and fixed

1. Prop I.3 mixed the cost-to-go HJB convention with the eikonal calibration;
   restated in forward form with an explicit calibration inequality.
2. Appendix C claimed the source coefficients are polynomial in J and finite
   at `J_deg`. They are not — they carry `Δ` in the denominator and have a
   simple pole there. The passage now distinguishes numerator from
   coefficients, which is *why* appendix B builds on the degenerate curve
   directly instead of taking a limit.
3. `eq:vaidya-full` mixed blocks normalised at `r_0` with one that was not;
   evaluating it as written was off by 158% at the sample point. Fixed with an
   explicit bracket plus a stated convention.
4. §4.1 gave a wrong relation between the two branch Hamiltonians (off by
   `−w/Ê`). Both now come from one support function.
5. `B` was used once in appendix B and never defined. It is `−a₃/(4a₄)`.
6. Three appendix-B coefficients silently set `m = 1`; `m` restored.

## 5. Verification

`paper1/verification/` — every algebraic step re-derived in **two** CAS
(Mathematica + SymPy), each check an exact symbolic zero, not a tolerance.

```
wolframscript -file verify_paper1_core.wls       # core
python3       verify_paper1_core.py              # cross-check, exits 1 on failure
wolframscript -file verify_appB_residues.wls
wolframscript -file verify_appB_blocks.wls
wolframscript -file verify_perlick_recovery.wls
wolframscript -file verify_optical_metric.wls
wolframscript -file verify_jacobi_conjugate.wls
wolframscript -file verify_caustic.wls
wolframscript -file verify_chord_bound.wls
wolframscript -file verify_apsidal_finiteness.wls
python3       verify_section44_closedform.py     # eq:vaidya-full end-to-end, 1e-39
```

`REPRODUCTION_section4.2.md` — §4.2 is phenomenology of integrated
trajectories, so it was reproduced by re-running the archived scripts rather
than proved. All three claims hold (penetration threshold, bounce
`r_min = 3.010492` at `v = 7.2857`, opposite winners in the timing).

**Caveat worth remembering**: those are reproductions of the author's own
integrations. They confirm the text reports what the scripts produce, not that
the scripts model what they intend to.

`build_paper1.sh` makes a clean checkout turnkey: runs the ten figure
generators, collects their output into `paper/Immagini`, compiles. Verified on
a fresh clone of `main`.

## 6. What is NOT proved, and is labelled as such

- **Three conjectures**, all of irreducibility type: the weight-two letters
  `W_jk` irreducible to weight one; the minimal dilogarithm-basis dimension;
  a canonical single-valued completion of the genus-two dilogarithm. If any
  were false, the closed forms would *simplify* — nothing becomes wrong. No
  result depends on them.
- `Φ''(J_c) ≠ 0` (the caustic being a simple fold) — numerical.
- **Declared out of scope**: the outgoing-Vaidya BVP, a third
  (Kodama-observer) clock, the ergosphere interior.
- **Never audited**: the bibliography (whether cited works say what is
  attributed to them), and the content of the 12 pre-existing figures.

## 7. Submission state (JMP)

Files: `paper1_JMP.pdf`, `cover_letter_paperI_JMP.pdf`, `figure_alt_text.txt`.

Portal answers already settled: Subject Area = *General Relativity and
Gravitation*; Special Topic = No; English editing = No; Related Manuscript =
`Unknown`; Dual submission = **No** (companion is a distinct manuscript);
Figure permissions = none previously published; Data availability = already
prepared; Funding = none, table left empty.

Suggested reviewers proposed: Perlick (Bremen), Caponio (Bari), Piccione (São
Paulo), Werner (Duke Kunshan), Hackmann (Bremen) — **emails were deliberately
not fabricated**; take them from institutional pages or recent arXiv PDFs, and
re-verify affiliations.

Zenodo: **v1.2.0, DOI `10.5281/zenodo.22035415`**, cited as
`Rosignoli2026CodeJMP`. Added as a *new* bib entry rather than editing
`Rosignoli2026Code`, which is shared with paper2 (under CQG review) and the
master; those still point at the older release deliberately.

## 8. Open items

1. **The shared `abelfunctions` bib note** still cites the old DOI
   `21782443`. The statement is true (the patched tree is archived there), but
   the JMP PDF then shows two DOIs for "the reproducibility package". Fixing
   it touches paper2 and the master. **Decision pending.**
2. **Alt text for tables** — the portal warns that on acceptance alt text will
   be required for figures *and tables*. Only the 12 figures are covered.
3. **`run_regression.py`** does not include the new verification suite; its
   regex extractor expects a different stdout format. `verify_paper1_core.py`
   already exits non-zero, so half the work is done.
4. **REVTeX port** — the manuscript is typeset in `iopart`. Fine for initial
   submission (AIP wants one compiled PDF); needed only at revision.
5. **`Φ_max < π` asymptotics** — now proved outright (I.E), so this is closed.
   What remains open is only `Φ'' ≠ 0`.

## 9. Traps that cost time — do not repeat

- Declaring `r` positive in SymPy makes the genus-degeneration locus *vanish*:
  the double root is at `r_d < 0`. Use plain symbols.
- SymPy will not combine `sqrt(A) sqrt(B)` into `sqrt(AB)` without positivity;
  use `powsimp(..., force=True)` and say why it is valid.
- Nested radicals: to check that an expression solves a quadratic, **substitute
  it into the quadratic** rather than comparing radicals. That is what finally
  closed the Perlick recovery symbolically.
- The chord comparison for `Φ_max < π` fails pointwise in `r` and holds in the
  optical radius `ρ = sqrt(J(r))`. Wrong variable ⇒ wrong conclusion.
- `grep -c Overfull` does **not** catch `Float too large for page`. Both must
  be checked; the latter is what pushed a caption off page 18.
