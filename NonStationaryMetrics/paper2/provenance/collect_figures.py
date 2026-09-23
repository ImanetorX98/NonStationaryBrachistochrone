#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collect_figures.py -- prove that every compiled figure is the current output of
its generator, and collect the ones that are not.

WHY THIS EXISTS.  The manuscript reads its figures from ``paper/Immagini`` via
``\\graphicspath``, while most generators write into their OWN directory:
``savefig(fig, HERE, ...)`` in KerrScripts, and so on.  Nothing connected the two.
An external audit put it exactly right: an inventory of candidate generators does
not prove the identity of the compiled figure with the generator's latest output.

It was not hypothetical.  A correction to the master penetration bar -- excluding
J = 0 from the cusp interval, to agree with the classification proposition -- was
made, the generator was re-run, and the PDF did not change, because the generator
had written to ``KerrScripts/`` and the paper was still compiling the older file
in ``paper/Immagini``.  The fix was invisible in the manuscript for as long as it
took to check the hashes.

WHAT THIS DOES.

  --check     for every figure the manuscript includes, compare the compiled
              file with the generator's output.  Report DIFFERS / MISSING /
              STALE and exit nonzero.  Writes nothing.
  --collect   copy the generator outputs over the compiled ones, reporting each
              copy.  Then --check passes.

Historical copies under submission_* and older trees are deliberately ignored:
those are frozen snapshots of what was filed and must not be touched.

Run:  python3 collect_figures.py --check
      python3 collect_figures.py --collect
"""
from __future__ import annotations

import hashlib
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER2 = HERE.parent
NSM = PAPER2.parent
ROOT = NSM.parent

TEX = PAPER2 / "paper2.tex"
COMPILED = NSM / "paper"          # what \graphicspath points at

# Where each generator actually writes its output.  Read off the savefig calls;
# keep this table next to the generators it describes, not in prose.
GENERATOR_DIR = {
    "fig_master_penetration_taut.png": ROOT / "KerrScripts",
    "fig_atlas_tau.png": ROOT / "KerrScripts",
    "fig_atlas_t.png": ROOT / "KerrScripts",
    "fig_separatrix_3traiettorie.png": ROOT / "KerrScripts",
    "fig_jpt_genus2.png": ROOT / "KerrScripts",
    "fig_jcm_capture.png": ROOT / "KerrScripts",
    "fig_thakurta_cuspide_ergosfera.png": COMPILED / "Immagini",
    "fig_phi_validation_true_dynamic.pdf": COMPILED / "Immagini",
    # Added after the September 2026 cross-audit.  The eight entries above were
    # the only ones under this check, so fourteen of the paper's twenty-two
    # figures were unwatched -- which is how a corrected generator once failed
    # to reach the compiled PDF without anything reporting it.
    "fig_bvp_estremi_fissi.png": ROOT / "KerrScripts",
    "fig_colormap_conforme_asimmetrico.png": ROOT / "KerrScripts",
    "fig_colormap_conforme_estremi_fissi.png": ROOT / "KerrScripts",
    "fig_colormap_spin_asimmetrico.png": ROOT / "KerrScripts",
    "fig_colormap_spin_estremi_fissi.png": ROOT / "KerrScripts",
    "fig_penetrate_bounce_orbit.pdf": ROOT / "KerrScripts",
    "fig_penetration_threshold.pdf": ROOT / "KerrScripts",
    "fig_tbranch_separatrix.png": ROOT / "KerrScripts",
    "fig_phi_closed_validation.pdf": NSM / "ThakurtaMetric",
    "fig_phi_penetrating.pdf": NSM / "ThakurtaMetric",
    "fig_thakurta_kerr_inversione_AJ.png": NSM / "ThakurtaMetric" / "Thakurtafigures",
    "fig_thakurta_kerr_plunge_t_tau.png": NSM / "ThakurtaMetric" / "Thakurtafigures",
    "fig_thakurta_kerr_rami.png": NSM / "ThakurtaMetric" / "Thakurtafigures",
    # fig_indicatrici_thakurta_kerr.pdf is produced by
    # ThakurtaMetric/genera_indicatrici_separate_codex.py, which writes through
    # its own path logic rather than paper_style.savefig; left out deliberately
    # rather than guessed at.
}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def included_figures() -> list[str]:
    tex = TEX.read_text()
    return sorted(set(re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", tex)))


def main(argv: list[str]) -> int:
    collect = "--collect" in argv
    figs = included_figures()
    print(f"{len(figs)} figures included by {TEX.name}; "
          f"{'collecting' if collect else 'checking'} against the generators\n")

    n_ok = n_diff = n_missing = n_untracked = 0
    for rel in figs:
        name = Path(rel).name
        target = COMPILED / rel
        if not target.exists():
            print(f"  MISSING   {rel}  (not in {COMPILED.name})")
            n_missing += 1
            continue
        src_dir = GENERATOR_DIR.get(name)
        if src_dir is None:
            # figure with no generator in the table: report, do not guess
            n_untracked += 1
            continue
        src = src_dir / name
        if not src.exists():
            print(f"  NO OUTPUT {name}  (generator has not been run in "
                  f"{src_dir.relative_to(ROOT)})")
            n_missing += 1
            continue
        if src.resolve() == target.resolve():
            n_ok += 1
            continue
        if sha(src) == sha(target):
            n_ok += 1
            continue
        if collect:
            shutil.copy2(src, target)
            print(f"  COLLECTED {name}\n"
                  f"              from {src.relative_to(ROOT)}\n"
                  f"              hash {sha(target)[:12]}")
            n_ok += 1
        else:
            print(f"  DIFFERS   {name}")
            print(f"              compiled  {sha(target)[:12]}  "
                  f"{target.relative_to(ROOT)}")
            print(f"              generator {sha(src)[:12]}  "
                  f"{src.relative_to(ROOT)}")
            n_diff += 1

    print()
    print(f"  {n_ok} up to date, {n_diff} differing, {n_missing} missing, "
          f"{n_untracked} with no generator in the table")
    if n_untracked:
        print("  (figures with no generator listed are not checked here; adding "
              "one to GENERATOR_DIR brings it under this check)")
    if n_diff and not collect:
        print("\n  rerun with --collect to copy the generator outputs over the "
              "compiled figures, then recompile the manuscript.")
    return 1 if (n_diff or n_missing) else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
