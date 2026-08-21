#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_crossrefs.py -- keep the author response's cross-references to the
manuscript honest.

The response document cites the revised manuscript by number ("Proposition 3",
"section 2.1", ...).  Those numbers are produced by LaTeX when paper2.tex is
compiled, and they move whenever a section or environment is inserted.  Typing
them by hand into a separate document is the same failure mode that produced the
drifting hashes and slopes of referee major comments 7 and 10 -- and it did
happen: a first draft of the response cited section 2.2, Protocol 2.2 and
Proposition 3.1 where the manuscript actually says 2.1, Protocol 1 and
Proposition 2.

So they are read out of paper2.aux instead.  Run this after compiling paper2.tex
and before compiling the response.

Writes: crossrefs.tex  (\\input by response_to_referees_CQG.tex)
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER2 = HERE.parent
AUX = PAPER2 / "paper2.aux"

# label -> macro name used in the response
WANTED = {
    "sec:domain":           "SecDomain",
    "sec:bvp":              "SecBvp",
    "sec:trichotomy":       "SecClassification",
    "sec:inversion":        "SecInversion",
    "sec:thakurta":         "SecThakurta",
    "sec:submersion":       "SecSubmersion",
    "app:validation":       "AppValidation",
    "lem:compact":          "LemCompact",
    "prot:status":          "ProtStatus",
    "prot:endpoint":        "ProtEndpoint",
    "prot:names":           "ProtNames",
    "prot:evidence":        "ProtEvidence",
    "prop:classification":  "PropClassification",
    "prop:exterior-sep":    "PropExteriorSep",
}

RE_LABEL = re.compile(r"\\newlabel\{([^}]+)\}\{\{([^}]*)\}")


def main() -> int:
    if not AUX.exists():
        sys.stderr.write(f"{AUX} not found -- compile paper2.tex first\n")
        return 1
    found: dict[str, str] = {}
    for m in RE_LABEL.finditer(AUX.read_text(errors="replace")):
        label, number = m.group(1), m.group(2)
        if label in WANTED:
            found[label] = number

    missing = [k for k in WANTED if k not in found]
    out = [
        "% ---------------------------------------------------------------",
        "% GENERATED FILE -- do not edit by hand.",
        "% Cross-reference numbers read out of paper2.aux by",
        "% provenance/make_crossrefs.py, so that the response cannot cite a",
        "% number the manuscript does not have.",
        "% ---------------------------------------------------------------",
    ]
    for label, macro in WANTED.items():
        if label in found:
            out.append(r"\newcommand{\ms" + macro + r"}{" + found[label] + "}")
        else:
            out.append(r"\newcommand{\ms" + macro + r"}{\textbf{??}}")
    (HERE / "crossrefs.tex").write_text("\n".join(out) + "\n")

    for label, macro in WANTED.items():
        mark = " " if label in found else "!"
        print(f"  {mark} {label:24s} -> {found.get(label, 'MISSING')}")
    if missing:
        sys.stderr.write(f"\n{len(missing)} label(s) missing from paper2.aux\n")
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
