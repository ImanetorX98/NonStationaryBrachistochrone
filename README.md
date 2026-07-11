# NonStationaryBrachistochrone

Numerical and symbolic codebase for the paper **"Constrained-worldline
brachistochrones in non-stationary spacetimes: Kodama rails, closed forms, and
plunge inversion"** (I. Rosignoli, 2026).

The relativistic brachistochrone — the fastest constrained "rail" worldline
between two events — is extended from stationary (Perlick 1991) to
**non-stationary** spacetimes. The rail invariant `−u·W = Ê` is an actively
controlled quantity, legitimate via Pontryagin's maximum principle, with the
selector `W` fixed by the hierarchy **Killing → conformal Killing → Kodama
vector**.

Three spacetimes form a ladder of increasing structure:

| Folder | Spacetime | Role |
|---|---|---|
| `FLRWmetric/` | FLRW | degenerate base (freezing, no branch splitting) |
| `VaidyaMetric/` | Vaidya | dynamical black hole, `W = ∂_v` the Kodama vector |
| `ThakurtaMetric/` | Thakurta–Kerr | conformal rotating black hole `g = A(η)² g_Kerr` |
| `KerrMetric/` | Kerr (equatorial) | Weierstrass separatrix, Doran continuation, fixed-endpoint BVP/colormaps |
| `PaperFigures/` | — | shared conceptual figures |
| `paper/` | — | LaTeX source (`main.tex`, `refs.bib`), figures (`Immagini/`), compiled `main.pdf` |

Each metric folder contains a `*Results.md` logbook documenting the analytic
results, validations, and — deliberately retained — the exploratory scripts and
negative results (which ansatz classes failed the held-out trajectory and
phase-space residual tests). `sumUp.md` is a session index mapping each figure to
its result and generating script.

## Main results

- **Kodama energy** conserved along the rail (`−u·K = Ê` identically), a
  *controlled* conservation law (drift `−m'(u^v)²/r` for free fall).
- Closed-form branch **Hamiltonians** for each metric (t / τ / conformal-time
  branches) from the elliptic indicatrix.
- **Equatorial separatrix** `φ(r)` in Weierstrass functions `℘, σ, ζ`, continued
  through the ergosphere in Doran coordinates.
- Analytic **cusp theorem** `φ − φ_e ∝ (r − r_e)^{3/2}` and the ergosphere
  trichotomy.
- **Plunge inversion** resolved as protocol-dependent: rotation and the conformal
  factor invert under a same-launch comparison; physical evaporation never does;
  no inversion survives fixed endpoints (`n_t/n_τ = E/f > 1`).

## Requirements

```
pip install -r requirements.txt
```
(`numpy`, `scipy`, `sympy`, `matplotlib`, `mpmath`, `tqdm`)

Scripts import the shared `paper_style.py` at the repository root, so run them
from within the repository (the metric-folder scripts add the root to
`sys.path`). Figures are written next to their generating scripts (e.g.
`ThakurtaMetric/Thakurtafigures/`).

## License

Released under the MIT License (see `LICENSE`).
