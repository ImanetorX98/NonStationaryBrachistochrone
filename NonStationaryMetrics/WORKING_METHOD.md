# Working method

Written down because it emerged rather than being planned, and because it is what
made the CQG revision productive rather than defensive. It is not a diary — the
diary is `POST_REVIEW_CHANGES.md` and the git log.

---

## 1. The rule that produced everything else

> **Every claim gets verified before it goes in, including the ones that feel obvious.**

The verification suite is nine scripts and ~130 checks. Roughly a third of them
found something. The ones that felt safest were not safer:

- "the cusp amplitude K" — my closed form was wrong, Mathematica's was right
- "Bishop's theorem says umbilical fibres" — two hypotheses missing
- "GMP 2001 is closest to our regime" — it is about null geodesics
- "$N(2M) = -4E^4$" — right by hand, but the script silently dropped the cubic term

**Corollary:** a check that passes on the first attempt deserves suspicion. Several
of the "passes" above were passing because the test was malformed.

## 2. Never assert what you have not read

Three times in this revision I wrote that something had been examined when it had
not. Each time it was caught — twice by the author, once by re-reading my own
output. In a document going to a referee who checks citations individually, that
class of error is expensive.

**Practice now:** if a claim about a source goes into the manuscript or the
response, the source is open in front of me. If it cannot be obtained, the text
says so explicitly rather than implying familiarity.

**The audit was worth its cost.** Of six sources read properly: two errors of ours
found, four results we could then use, one risk closed. That rate justifies
reading the rest.

## 3. Generate, never transcribe

Any number that appears in two places will diverge. It happened here to
checksums (three different values for the same file) and to fitted slopes (figure
and table disagreeing, and both stale).

**Practice:** one archived command produces the number, writes a LaTeX fragment,
and the manuscript `\input`s it. Cross-references from the response into the
manuscript are read out of `paper2.aux` for the same reason — a first draft cited
section 2.2 where the paper says 2.1.

Scripts: `provenance/make_provenance.py`, `provenance/make_crossrefs.py`.

## 4. Concede the real point, then disagree narrowly

The strongest replies in the response are the ones that start by granting
something. The referee who claimed our results were "a trivial application" was
pointing at a real structure — we found three genuine connections by taking the
claim seriously, and the disagreement that remains is one narrow sentence instead
of a confrontation.

**Practice:** before declining, spend the effort to find what is right about the
objection. If nothing is, say so in one sentence and move on. Never argue about
tone; never correct the referee's slips (a misattributed author, a misspelled
name) — cite correctly and let the record speak.

## 5. Separate what is claimed from what is shown

Four protocols in the manuscript now do this: type (E)/(C)/(A) for every curve,
the endpoint convention, the nomenclature table, and the four levels of evidence.
The last is the load-bearing one — exact symbolic identity, independent numerical
agreement, recovery of established limits, computer-assisted certificate — because
it forces every claim to name its status.

**Practice:** when a result cannot be proved at the level you want, say which level
it *is* at rather than softening the language.

## 6. Freeze what has been submitted

`paper1/submitted_JMP_2026-08/` is the text exactly as sent to JMP, with a README
naming its one known defect. The working copy moves; the submitted text stays
recoverable. This also lets us tell CQG that the conceptual corrections are
public on Zenodo without disturbing the JMP process.

## 7. Traps that cost real time

| Trap | Symptom |
|---|---|
| Wolfram terminates a statement at a newline after a complete expression | a cubic term written on a continuation line beginning with `-` is **silently dropped**; the parameter scan came out exactly inverted |
| Mathematica comments nest | `(*)` inside `(* … *)` opens a comment never closed |
| `grep -c Overfull` | does **not** catch `Float too large for page` — different warning class |
| a generated `\newcommand` containing `$…$` | must not be used inside `$…$` |
| `NSolve[… && 2 < r < 10^6, r, Reals]` | can return unevaluated; for a polynomial take exact roots and filter |
| SymPy `positive=True` | can make a degeneration locus vanish |
| `latexdiff` | `ulem` breaks across `\emph` → use `--type=CFONT`; markers land before `\hline` in tabulars; one `\multicolumn` row needs hand repair |
| off-by-one in a bracketing helper | `r_turn()` returned the grid point *inside* the forbidden region — a curve vanished from a published panel with no error |

---

## Progress

### Deliverables

| File | State |
|---|---|
| `paper2.pdf` | 66 pp — regenerate before submission |
| `paper2_highlighted.pdf` | regenerate (`latexdiff`, three known repairs) |
| `response_to_referees_CQG.pdf` | 10 pp — regenerate after final edits |
| `submission_CQG_R1.zip` | rebuild from final source |
| `paper1_JMP.pdf` | 54 pp, revised working copy |
| `paper1/submitted_JMP_2026-08/` | frozen |
| `.zenodo.json` | 1.5.0; tags `v1.3.0`, `v1.4.0` pushed |

### Review response

All 12 major + 10 minor comments of referee 1, and referee 2's bibliography,
answered. Five errors caught by the review, four found while chasing them, five
new results. See `POST_REVIEW_CHANGES.md`.

### Bibliographic audit — 7 of 11 sources read

| Source | Read | Produced |
|---|---|---|
| Lecian, T2 preprint | full | analysis of the "trivial application" claim |
| Lecian, solitons/geodesic paths | full | the $\lambda+\mathrm{Ric}(u,u)$ identity |
| Meena–Zawadzki | core | Bishop confirmed; Prop. 1 and Cor. 7 |
| GMP 2000 (Morse, massive) | core | index = geometric index; Maxwell topologically |
| Caponio et al. | core | semi-holonomic constraint; **risk closed** |
| GMP 1998 (timelike Fermat) | core | third-party reading of Perlick; PS caveat |
| Lecian, curvature eigenvalues | pp. 1–11 | the weighted-seed correspondence |
| GMP 2001 | pp. 1–3 | characterisation corrected |
| Axioms 15(4) 267 (Yamabe flow) | intro | confirms the weight construction; Birkhoff framing |
| Axioms 14(12) 896 (Cauchy) | intro | independent boundary convergence; incompleteness contrast |
| Lecian, 2-planes solitons | intro | energy-dependence sharpens our threshold |

### Open

1. Finish the audit — the four unopened sources, one of which we cite.
2. Regenerate all four CQG deliverables from the final source.
3. Zenodo release from `v1.5.0`; update the DOI in manuscript and cover letter.
4. The ~90 pre-existing references have never been audited.
