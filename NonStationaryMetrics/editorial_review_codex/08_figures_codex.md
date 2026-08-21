# Editorial notes — Figures and tables

## Global assessment

All 38 pages were rendered and inspected. No clipping, missing glyphs or corrupted plots were found. The main issue is not technical failure but hierarchy and legibility. The manuscript has the look of an internal validation dossier: large Matplotlib panels, verbose titles inside plots, long captions and substantial page area devoted to diagnostics.

## Specific observations

- **Figure 1 (page 7):** the three indicatrices are too small, with tiny legends and axis labels. The Thakurta–Kerr panel also imports Paper-II content. Redraw as a clean vector schematic or restrict it to FLRW/Vaidya.
- **Figures 2–3 (pages 11–12):** individually clear but too large for the modest FLRW message. Combine them as panels.
- **Figures 4–6:** useful for Vaidya phenomenology, though typography and in-plot titles should be standardised.
- **Figure 7 (page 18):** dominated by a disclaimer. It shows a limitation more than a main result and interrupts the derivation; move it to supplement/appendix or replace with a short textual statement.
- **Figures 8–9 (pages 23–24):** appear between pieces of the Conclusions. This damages closure and gives the formal negative-rate continuation undue prominence. Place all result figures before the Conclusions.
- **Figure 10 (page 25):** scientifically important because it distinguishes `O(epsilon)` from `O(epsilon^2)` residual scaling. Keep it, enlarge labels, state fit range and uncertainty, and avoid referring to “Eq. (40)” if the displayed manuscript equation is C.1.
- **Figure A1:** useful but can be smaller.
- **Figure A2:** captions claim minimisation from a restricted perturbation test; qualify the language.

## Tables

- **Table 1:** remove Paper-II-only notation and group the remaining symbols by geometry/control/clock.
- **Table A1:** acceptable but the explanatory caption is too long for a five-number table.
- **Table A2:** visibly cramped. The final row collides with the script/residual columns and the table mixes prose separators with data rows. Rebuild with `tabularx`/`longtable`, one script per row, and a separate column for scope or equation.

## Style checklist

- Export line plots as vector PDF/EPS.
- Use 8–12 pt text at final size, consistent with IOP guidance.
- Use line styles/markers as well as colour.
- Remove narrative headlines from inside plots; put the conclusion in the caption.
- Use `(a)`, `(b)` panel labels and concise, self-contained captions.
- Keep notation identical between figure, caption and text (`v`, `tau`, `mu`, `m'`).
