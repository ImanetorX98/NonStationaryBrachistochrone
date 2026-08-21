CQG-116884 -- revised source bundle
===================================

Manuscript: "Brachistochrones in conformal Kerr spacetimes: control domain,
exterior separatrices, and adiabatic response"

Build:
    latexmk -pdf paper2.tex
(or pdflatex x3; paper2.bbl is included, so BibTeX is not required)

Contents
--------
  paper2.tex                 main source
  paper2.bbl                 resolved bibliography
  paper2.pdf                 compiled from exactly these files
  provenance/                generated macros (see below)
  Immagini/                  all 22 figures
  iopart.cls, iopart10.clo, iopart12.clo, iopams.sty
                             IOP class files, included for self-containment

On provenance/adiabatic_slopes.tex
----------------------------------
This file is GENERATED, not written by hand. The fitted adiabatic slopes, the
epsilon-window and the SHA-256 checksums printed in the manuscript are emitted
by a single archived command in the reproducibility package,

    python3 paper2/provenance/make_provenance.py

which runs each validation script, parses the slopes and their covariance
sigmas, reads the epsilon-window out of the source, and hashes every script
involved. The manuscript \input's the result, so the printed values cannot drift
from the archive. This replaces the hand-transcribed values of the previous
version, which had drifted -- see referee 1, major comments 7 and 10.

Reproducibility package
-----------------------
  https://github.com/ImanetorX98/NonStationaryBrachistochrone

  Symbolic verification suite: paper2/verification/ -- 9 scripts, 123 checks,
  each an exact symbolic zero or a high-precision convergence test.
  Manifest mapping each artefact to one command and one checksum:
  paper2/provenance/MANIFEST.tsv
  Environment lock (Python stack and CAS version):
  NonStationaryMetrics/requirements-lock.txt
