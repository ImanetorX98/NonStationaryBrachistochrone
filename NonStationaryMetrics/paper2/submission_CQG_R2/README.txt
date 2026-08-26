CQG-116884 -- revised submission (R2)

Self-contained LaTeX source for the revised manuscript.  Compiles from a clean
directory with:

    latexmk -pdf paper2.tex

Contents
    paper2.tex                 revised main file, clean (no highlighting)
    paper2.pdf                 clean PDF built from this source
    paper2.bbl                 resolved bibliography
    refs.bib, companionI.bib   bibliography sources
    provenance/                generated fragment carrying the fitted slopes,
                               epsilon-window and provenance digest, which the
                               manuscript \input's so printed numbers cannot
                               drift from the archive
    Immagini/                  the 22 figures included by the manuscript
    iopart.cls, iopams.sty,    IOP class files
    iopart10.clo, iopart12.clo

The full reproducibility archive -- verification scripts, figure generators and
the laboratory design study -- is at doi:10.5281/zenodo.22110181 (v1.6.0).
