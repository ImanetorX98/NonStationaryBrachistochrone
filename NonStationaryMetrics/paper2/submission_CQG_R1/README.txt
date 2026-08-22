CQG-116884 -- revised source bundle
===================================

"Brachistochrones in conformal Kerr spacetimes: control domain, exterior
separatrices, and adiabatic response"

Build:
    latexmk -pdf paper2.tex
paper2.bbl is included, so BibTeX is not required.

Contents
  paper2.tex, paper2.bbl, paper2.pdf
  provenance/adiabatic_slopes.tex   GENERATED -- see below
  Immagini/                         all figures
  iopart.cls, iopart10.clo, iopart12.clo, iopams.sty

On provenance/adiabatic_slopes.tex
  This file is generated, not written by hand. The fitted adiabatic slopes, the
  epsilon-window and the SHA-256 checksums printed in the manuscript are emitted
  by one archived command in the reproducibility package,
      python3 paper2/provenance/make_provenance.py
  which runs each validation script, parses the slopes and their covariance
  sigmas, reads the epsilon-window out of the source, and hashes every script
  involved. The manuscript \input's the result, so the printed values cannot
  drift from the archive. This replaces the hand-transcribed values of the
  previous version, which had -- see referee 1, major comments 7 and 10.

Reproducibility package
  https://github.com/ImanetorX98/NonStationaryBrachistochrone
  Symbolic verification suite: paper2/verification/ -- 11 scripts, each check an
  exact symbolic zero or a high-precision convergence test.
  Artefact -> command -> checksum map: paper2/provenance/MANIFEST.tsv
  Environment lock (Python stack and CAS version):
  NonStationaryMetrics/requirements-lock.txt

Companion paper
  The conceptual corrections to the minimum principle discussed in the Author
  Response are in the revised companion Paper I, archived on Zenodo. The version
  under consideration elsewhere is preserved unchanged in the repository under
  paper1/submitted_JMP_2026-08/.
