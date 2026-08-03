# Paper I — Codex handoff, round 2 (post 34-page build)

**This supersedes `HANDOFF_CODEX_production.md`.** That production batch is done: the
authoritative build is now `paper1/paper1.pdf`, **34 pages, compiled 4 Aug 2026** (Fig 1
reduced to two spherical panels, Fig 9 caption with fit data, PDF metadata, References
heading, no overfull, no near-empty pages). The updated CQG editorial review
(`output/pdf/CQG_editorial_review_latest_codex.pdf`) verdict is:

> *Send for external peer review as a CQG Research Paper — after a short mandatory pass.
> "Almost, but not today." Correct the ten mandatory items, synchronise the Zenodo
> records, and submit. Further broad rewriting is not justified.*

Apply the edits below **to the current 34-page build**, then rebuild once
(`pdflatex → bibtex → pdflatex → pdflatex`).

---

## GOLDEN RULES — read before editing

- **Do not reopen the mathematics.** None of these items requires changing a derivation,
  a sign, a residue, or the off-shell construction. They are notation, wording,
  self-containedness, bibliography and archive-metadata repairs.
- **Never infer a DOI** from a title. Where a DOI/title/version decision is needed, it is
  an **author decision** — leave a clearly marked `% TODO(author)` and do not fabricate.
- **Do not touch the PRESERVE list** (§ end of this file). In particular keep verbatim:
  Theorem I.5 in the main text; equation (22) as the Figure-9 target; Figure 9's exact
  regression values and script hash; the "Proved / Conditional and numerical / Outlook"
  conclusion; the conjectural/open labels; the Vaidya boundary source
  `S_D = [r p_r] − λ`; the conditional-HJB caveat; the turning-point non-uniformity
  statement.
- After each notation change, **grep the whole source** (tex, captions, figure scripts
  that emit in-plot text) so the symbol is consistent everywhere, not just at first use.

---

## MANDATORY (M1–M10) — must fix before submission

### M1 — Undefined `T` at the start of §4.1
On PDF p.11, §4.1 opens with `−g(T,T) = (f u^v − 1)^2 / Ê^2` then "hence `−u·K = Ê`".
`T` is never defined (global search finds no other occurrence).
**Fix:** if `T` is an auxiliary tangent, define its components and role *before* the
identity. If the controlled constraint is already `−u_v = Ê` with `K = ∂_v`, the auxiliary
identity is unnecessary — **delete it**. Clarity repair only; **do not change eq (12)**.

### M2 — Unannounced switch from `Ê` to `E`
Control charge is `Ê` through eqs (1), (7)–(18) and Table 1. §4.4 then defines the sextic
with `DE = (E^2−1)r + 2m`, and Appendices A–B continue with `E`. No sentence states
`E ≡ Ê`.
**Fix:** state `E ≡ Ê` at the first switch, or use `Ê` throughout. See also S/notation:
the product-like symbol `DE` should become a single named quantity `D_E` (or `Δ_E`), or a
sentence must state it is one defined polynomial factor, not a product of two variables.

### M3 — Two incompatible definitions of `ε`
Table 1 says `ε = ṁ`. Appendix C and Figure 9 use `ε = M ṁ / m` (dimensionless).
**Fix:** reserve `ṁ = dm/dv` for the rate and `ε = M ṁ / m` for the dimensionless
bookkeeping parameter. If units `G = c = M = 1` are used in a numerical subsection, say so.
If the slow coordinate is `s = v/M`, print that definition **before eq (C.1)**. Essential:
the central convergence plot (Fig 9) is plotted against `ε`, so the two must not conflict.

### M4 — Collision between the sextic `S` and the source `S_D`
Eq (22) defines `S_D(λ) = ∫₀^λ 𝒟H dλ'`. A later paragraph nonetheless writes
`S = ∫ ΘH_v dλ'`, and the Fig 9 caption says the on-shell term drops "the costate
integral S". But `S` is the genus-two spectral sextic under every radical.
**Fix:** use `S_D` for the source in the paragraph, the Figure 9 caption, and any
script-facing prose. Keep `S` for the sextic only. (Small visual distinction, entirely
different objects — this is exactly the collision that breeds sign/factor confusion.)

### M5 — Path parameter, cost clock and augmented state
§2.1 says the state is only `x = (r, φ)` and the evolution parameter is the selected clock
`s ∈ {t, τ, η, v}`. But for a Vaidya proper-time cost the metric depends on advanced time
`v`; if `v` is not the path parameter it must be in the state, and if `v` is the path
parameter then `τ` is a cost and must not be listed as the evolution parameter of that
Hamiltonian. Also the HJB equation currently "minimises over `u`", but `u` already denotes
the four-velocity (and later the retarded coordinate).
**Fix:** add one short protocol paragraph or table distinguishing (1) the monotone path
parameter, (2) the clock one-form being minimised, (3) the full state including a time
coordinate when the metric is non-autonomous, (4) the control — **rename the control to
`θ`** (or another unused symbol) everywhere in the HJB statement.

### M6 — False claim that a nome series is printed
Appendix C says the genus-one nome series `g^{(1)}` is "printed in Section 4.4". It is not
present there or anywhere. §4.4 also says the full Kronecker–Eisenstein evaluation is in
Appendix B, while Appendix B gives the Weierstrass block form and residue data but **not** a
standalone q-series recipe.
**Fix (choose the shorter, consistent option):** state that the PDF fixes the algebraic
reduction and the special-function representation, while period construction and numerical
q-series evaluation require the archived routine. Do **not** simultaneously claim full
reconstruction from the PDF *and* reliance on an unprinted routine. Remove the "printed
`g^{(1)}`" claim.

### M7 — Rename the Appendix-C "Separatrix limit"
Main text and Appendix B correctly say `J_deg` gives a double root `r_d < 0` off the
physical arc — not a physical external separatrix. Appendix C's last paragraph is still
titled "Separatrix limit (J → J_deg)" and repeatedly uses "separatrix correction / form".
**Fix:** rename to **"Algebraic genus-degeneration limit"** and use "degeneration
correction / form" for Vaidya. Reserve `J_sep` and the word "separatrix" for an accessible
dynamical boundary (the rotating companion, Paper II).

### M8 — "Spatial-curvature" wording for a spatially flat FLRW
The abstract says the correction "separates the spatial-curvature effects of the FLRW base
from the radial mass-flow effects"; the introduction repeats the spatial-curvature vs
radial-flow contrast. The FLRW metric used is **spatially flat**, so this is a
contradiction. The real logic: homogeneous time dependence produces freezing *without* a
spatial gradient, whereas Vaidya introduces radial gradients and mass-flow memory.
**Fix:** use "homogeneous expansion" for FLRW and "radial spatial-gradient and mass-flow
effects" for Vaidya, in **both the abstract and the introduction** (and the cover letter if
the phrase appears there).

### M9 — Reference [12] / cover letter title vs public DOI  *(author decision)*
Manuscript and cover letter cite DOI `10.5281/zenodo.21740000` under the title
"Brachistochrones in conformal Kerr (Thakurta) spacetimes: ergosphere trichotomy,
separatrices, and adiabatic response". The **registered Zenodo/DataCite title** is
"Relativistic Brachistochrones in Conformal Kerr Spacetimes: Adiabatic Theory, Ergosphere
Transitions and Hyperelliptic Closed Forms".
**Fix (author must choose one):** (a) update the citation and cover letter to the
*deposited* title; or (b) publish a new companion Zenodo version with the intended title
and cite that version. **Bibliographic title must match the retrievable object.** If
"Thakurta" stays in a BibTeX title, protect its capitalization with braces `{Thakurta}`
(the current PDF lower-cases it). Mark with `% TODO(author)` — do not invent metadata.

### M10 — Public snapshots older than the final build  *(author decision)*
The cover letter calls DOI `10.5281/zenodo.21739998` "the present paper", but that public
record contains `paper1-10.pdf`, not this 34-page build. The software DOI
`10.5281/zenodo.21707378` is active/v1.0/MIT but was issued before the final Figure 1 and
Figure 9 production wrappers.
**Fix (author action):** publish immutable final versions (new Zenodo version for Paper I,
and if needed Paper II and the code package), then update the cover letter / data
statement / Table A2 hash from the *verified* new metadata; **or** soften the wording so no
older record is called "the present paper". Do not mutate a DOI claim by inference. Mark
`% TODO(author)`.

---

## STRONGLY RECOMMENDED (S1–S10) — do where practical

- **S1** Reduce companion-Paper-II forward references to ≈5 strategically placed mentions
  (there are currently ~two dozen; prune later repetitions, keep the formal companion
  citation in the introduction).
- **S2** Remove the "retained separately" sentence from Figure 1's caption (it describes
  file management, not reader information); the caption should end after explaining FLRW
  and Vaidya.
- **S3** Define **PMP**, **HJB** and **BVP** at first use (IOP asks acronyms be defined in
  the abstract and main text).
- **S4** Rename the generic tangent argument in `F(x, v)` to avoid collision with Vaidya
  advanced time `v` (use a neutral symbol).
- **S5** Move essential parameters from Figures 4–8 plot titles into the captions (e.g.
  `(E, J, m', v_1, r_1)`); use line styles as well as colour.
- **S6** Add line styles/markers and vector export to remaining multi-curve figures (2–8,
  A1–A2) for grayscale robustness — *without* changing numerical content.
- **S7** Make references **[15]**, **[31]**, **[32]**, **[39]** independently locatable in
  the *rendered* list:
  - [15] Caponio–Javaloyes–Sánchez, show the title + issue 1501 + DOI `10.1090/memo/1501`
    (source BibTeX already has them; the `iopart-num` style is suppressing them).
  - [31] Brown–Levin, must visibly read "Multiple elliptic polylogarithms", arXiv:1110.6917.
  - [32] match the official Zenodo software-record title (keep v1.0 + exact commit only if
    that is the archive supporting the submitted paper).
  - [39] abelfunctions — the archived snapshot matters more than the live GitHub URL;
    typeset `is_*` unambiguously and cite the precise archive.
  Note: the visible-rendering problem is the `iopart-num` style dropping DOI/title fields
  that ARE in the source; fix at the rendered-entry level, not by inventing data.
- **S8** Run a verified DOI pass (no guessing). Prioritise the direct-brachistochrone,
  Kodama/Vaidya, dynamical-horizon and software references.
- **S9** Compress repeated special-function and Paper-II implementation prose in §4.4 and
  Appendices B–C (keep eqs (19)–(24), Theorem I.5, the boundary source, the representation
  hierarchy and Figure 9; move repeated degeneration/companion discussion to appendices).
- **S10** Verify a consistent author **affiliation** across manuscript, ORCID, cover letter
  and Zenodo. The code-release DataCite metadata currently carries a *different* affiliation
  string from the preprint records — align them (`% TODO(author)` if it needs a Zenodo edit).

---

## NOTATION & TERMINOLOGY AUDIT (review §9) — consolidated

| symbol | issue | required action |
|---|---|---|
| `Ê` vs `E` | mixed | define `E ≡ Ê` once or use one symbol (**M2**) |
| `ε` | both `ṁ` and `Mṁ/m` | reserve `ṁ` = rate, `ε` = dimensionless; scale slow clock (**M3**) |
| `S` vs `S_D` | sextic vs source | use `S_D` for source incl. Fig 9 (**M4**) |
| `T` | only in `g(T,T)` | define or delete (**M1**) |
| `u` | 4-velocity + HJB control + retarded time | keep `u^a` for 4-velocity, `u` for retarded time only, `θ` for control (**M5**) |
| `v` | advanced time + tangent arg of `F(x,v)` | rename the generic tangent (**S4**) |
| `s` | evolution clock + slow parameter | distinguish path param / dimensionless slow coord / optimised clock (**M5**) |
| `λ` | elapsed inward clock | define orientation at first use in the theorem; repeat compactly in App C |
| `J_pen` | physical penetration threshold | **preserve** — distinction is correct |
| `J_deg` | algebraic double root at `r_d<0` | call it genus degeneration, not separatrix (**M7**) |
| `J_sep` | possible spectral separatrix (Paper II) | do **not** introduce in Paper I unless needed |
| `DE` | polynomial `(E²−1)r+2m` | prefer `D_E`/`Δ_E`, or state it is one named factor |
| `F` | optical cost / integrand / rational fn | a subscript would reduce load in App B–C |
| PMP/HJB/BVP | undefined acronyms | define at first use (**S3**) |
| "closed form" | sometimes needs periods/code | say "closed special-function representation" where numerical layer is real |
| "separatrix" | physical vs algebraic | reserve for a dynamical boundary; "genus-degeneration locus" for Vaidya `r_d<0` |

---

## FINAL PRODUCTION CHECK (review §12.3 / §14) — after edits

- **P1** Full `pdflatex/bibtex` build: **0 undefined refs, 0 overfull boxes**.
- **P2** Grep PDF + build files for `??`, stale equation numbers, duplicate labels.
- **P3** Confirm Theorem I.5 still in main text and eq (22) is still the Figure-9 target.
- **P4** Confirm Table A2 and eq (B.2) remain within margins (do **not** widen Table A2 —
  keep monospaced paths; it currently wraps correctly).
- **P5** Re-render and inspect every page (figures, bibliography, page transitions).
- **P6** Confirm PDF title/author/subject/keywords and embedded fonts.
- **P7** Confirm final DOI links resolve to the exact titles/versions stated (**M9/M10**).
- **P8** Prepare the chosen single- or double-anonymous package consistently. (If
  double-anonymous review is selected, a separate anonymised PDF is required — the present
  file is not anonymised.)

---

## RECOMMENDED SEQUENCE (review §14)

1. **Definitions first:** `T`, `E/Ê`, `ε`, `S/S_D`, path/cost/state conventions, control
   symbol. Recompile; check eqs (12), (22), (C.1)–(C.6).
2. **Special-function claims:** remove the false printed-`g^{(1)}` claim, pick the honest
   reconstructibility wording, rename the algebraic degeneration.
3. **Public records (author):** decide final titles, create new Zenodo versions as needed,
   then update ref [12], ref [32] and the cover letter from *verified* metadata.
4. **Bibliography pass:** fix [15], [31], [39] rendering; add only verified permanent
   identifiers; keep all 43 entries cited.
5. **Trim Paper-II leakage** (S1) without touching Paper-I derivations.
6. **Polish captions/plots** (S5, S6) — parameters self-contained; vectorise only if it
   does not risk numerical content.
7. **Build & render**; repeat the clean 34–35-page checks; inspect every page.
8. **Submit:** upload the complete PDF, supplementary/code info, and the corrected cover
   letter; select the peer-review anonymity model.

---

## PRESERVE — do NOT reopen (review §13)

Controlled invariant `−u·W = Ê` as the central definition; the selector hierarchy with its
spherical-symmetry limitation; the operational rail-force / control-power discussion; the
compact strictly-convex indicatrix assumptions; the fixed-launch / fixed-target /
free-arrival protocol; endpoint-transversality vs interior-first-integral distinction; the
costate vs mechanical-angular-momentum distinction; the conditional HJB and cut/Maxwell
caveats; FLRW as the homogeneous branch-degenerate base; the unnormalised Kodama selector
for ingoing Vaidya; the formal status of negative-rate ingoing continuation; **Theorem I.5**
in the main text and its Appendix-C derivation; the turning-point non-uniformity statement;
the Vaidya boundary source `S_D = [r p_r] − λ` (after notation cleanup); **Figure 9 and its
exact regression values**; explicit conjectural/open labels; the "Proved / Conditional and
numerical / Outlook" conclusion; the data DOI, script map and transparent AI disclosure
(update only for exact versioning).

---

## Paper I title (unchanged — for `pdftitle` and citation)

> **Controlled-rail brachistochrones in non-stationary spacetimes: conformal symmetry,
> Kodama energy, and Vaidya dynamics**

**Note on Paper II (see M9):** the companion is currently cited as *"Brachistochrones in
conformal Kerr spacetimes: ergosphere trichotomy, separatrices, and adiabatic response"*,
but the **deposited Zenodo title** is *"Relativistic Brachistochrones in Conformal Kerr
Spacetimes: Adiabatic Theory, Ergosphere Transitions and Hyperelliptic Closed Forms"*.
These must be reconciled before submission (author decision — do not guess).
