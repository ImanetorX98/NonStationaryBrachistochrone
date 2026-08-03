# Paper I — production handoff for Codex

**Status at handoff:** all editorial MUST-FIX and the cheap LaTeX SHOULD-FIX are done
and committed. Paper compiles: **37 pp, 0 undefined refs, 1 pre-existing overfull**
(line ~1188, a deep-appendix formula display — see task 5). What remains below is the
**figure + PDF production batch**, which needs Python script re-runs and cannot be done
at the LaTeX level alone.

- Main file: `NonStationaryMetrics/paper1/paper1.tex` (documentclass `iopart`, style `iopart-num`).
- Figure scripts: `NonStationaryMetrics/VaidyaMetric/`, `NonStationaryMetrics/FLRWmetric/`,
  `NonStationaryMetrics/ThakurtaMetric/`.
- Build: `cd NonStationaryMetrics/paper1 && pdflatex paper1 && bibtex paper1 && pdflatex paper1 && pdflatex paper1`
- **Do NOT alter** (frozen by the editorial review): Theorem I.5 in the main text; the
  conditional-minimization caveat; the turning-point non-uniformity statement; the
  special-function hierarchy and its conjectural admissions; the Figure 9 slope-2
  validation; references [8]–[11], Zenodo [32], the AI disclosure.
- **Never infer a DOI.** If a DOI is unknown, leave it out.

---

## Task 1 — Figure 1 (`fig:indicatrici`): enlarge, remove third (TK) panel

**Why:** Figure 1 renders too small, and its third panel shows the rotating
Thakurta–Kerr indicatrix, which is Paper-II material (leakage). Removing the TK panel
frees width so the remaining panels can be enlarged.

- Script: `NonStationaryMetrics/ThakurtaMetric/genera_figure_thakurta.py`
  (this is the generator for `fig_indicatrici`).
- Drop the TK/rotating panel; keep only the spherical (FLRW + Vaidya) indicatrices that
  Paper I actually uses.
- Enlarge remaining panels to fill the text width; target legible axis labels/ticks at
  final print size (labels ≥ 9 pt, ticks ≥ 8 pt, legend ≥ 8 pt in the *rendered* figure).
- Re-export the vector PDF the `\includegraphics` in `paper1.tex` points to (search the
  tex for `fig_indicatrici`). Keep the same output filename so the tex needs no change.
- In `paper1.tex`, check the Figure 1 caption still matches (no dangling reference to the
  removed rotating panel).

## Task 2 — Figure 9 caption: add reproducibility numbers

**Why:** the review wants the slope-2 validation figure to be self-documenting.

- Find the Figure 9 block in `paper1.tex` (it is the adiabatic slope-2 / first-order
  validation figure; grep `slope` and `2.00`).
- Add to the **caption** (not the plot): the ε (rate) fit interval, the number of sampled
  points, the residual norm of the fit, and the fitted slope(s). Pull the exact numbers
  from the generating script's output (`vaidya_first_order_offshell.py`, whose SHA-256
  prefix `5dc9fe177f2f` is already cited in Table A2). **Do not change the slope-2 result
  itself** — only report the numbers already produced.

## Task 3 — Page reflow: remove the near-empty pages 21 / 31 / 35

**Why:** three pages are almost blank (float-placement / `\clearpage` artefacts),
padding the paper. Goal: drop 2–3 pages without moving any float out of its section.

- Recompile, then `pdftotext paper1.pdf - | ...` or open the PDF to confirm which pages are
  near-empty in the *current* build (numbering may have shifted after the last edits).
- Fix by relaxing float placement: prefer `[htbp]` over `[t]`/`[p]`, move or soften the
  `\clearpage` before `\section{Conclusions}` if it is the cause, and let long float
  runs interleave with text. Do **not** shrink figures or delete content.
- Re-verify: 0 undefined refs, no float appears before it is first referenced, page count
  drops by 2–3.

## Task 4 — PDF metadata

**Why:** submission PDF currently has empty document metadata.

- `paper1.tex` loads `hyperref`. Add a `\hypersetup{...}` in the preamble (after the
  `\usepackage{hyperref}` line) with:
  - `pdftitle` = **the paper title** (see bottom of this file),
  - `pdfauthor` = `Iman Rosignoli`,
  - `pdfsubject` = `Classical and Quantum Gravity submission`,
  - `pdfkeywords` = e.g. `relativistic brachistochrone; controlled rail; Kodama energy;
    Vaidya spacetime; FLRW; optimal control; Randers–Zermelo geometry`.
- Confirm with `pdfinfo paper1.pdf` after rebuild.

## Task 5 — Final consistency pass

- Kill the remaining overfull hbox (~line 1188): it is a wide math display in the deep
  appendix. Break the line, use `\!`/`\,` spacing, or wrap in a `\resizebox`/`split` —
  cosmetic only, do not change the mathematics.
- Full pass: `grep '??' paper1.log`, check for stale/duplicated labels, renumber sanity
  (all `\ref`/`\eqref` resolve), Table A2 last row still not colliding with the SHA-256
  line, References heading present exactly once.
- Final build must report: **0 undefined, 0 overfull, N pages** (N = 34–35 after reflow).

## Task 6 (optional, greyscale-safe figures)

If time allows: in the figure scripts remove verbose in-plot titles (e.g. the Kodama
figure title), add distinct line styles + markers so the figures survive greyscale
printing, and export vector PDF. Regenerate and re-include. Lower priority than 1–5.

---

## Paper title (give this to Codex, and use it for `pdftitle`)

> **Controlled-rail brachistochrones in non-stationary spacetimes: conformal symmetry,
> Kodama energy, and Vaidya dynamics**

(Companion Paper II: *Brachistochrones in conformal Kerr spacetimes: ergosphere
trichotomy, separatrices, and adiabatic response*.)
