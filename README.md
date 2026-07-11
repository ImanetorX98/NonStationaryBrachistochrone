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
| `ThakurtaMetric/` | Thakurta–Kerr | conformal rotating compact object `g = A(η)² g_Kerr` |
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

## Main scripts and the figures they produce

Each figure is written to a `*figures/` subfolder next to its script (e.g.
`ThakurtaMetric/Thakurtafigures/`); the versions used in the paper are collected
in `paper/Immagini/`. The symbolic scripts (`*_sympy.py`) derive and check the
Hamiltonians and identities and print residuals rather than plotting. The
`colormap_*` scripts integrate the Hamiltonian flow over a grid and are heavy
(tens of minutes to a couple of hours; they show a `tqdm` progress bar).

### `FLRWmetric/` — degenerate base
| Script | Produces / does |
|---|---|
| `genera_figure_flrw.py` | `fig_flrw_cinematica` (rail speed, freezing at `a=Ê`), `fig_flrw_worldlines`, `fig_flrw_variazionale` (branch degeneracy: `t` and `τ` minimized by the same comoving line) |
| `flrw_hamiltonians_sympy.py`, `flrw_optical_metric_sympy.py`, `flrw_brachistochrone_sympy.py` | symbolic rail Hamiltonians and optical index |

### `VaidyaMetric/` — Kodama rail
| Script | Produces / does |
|---|---|
| `genera_figure_vaidya.py` | `fig_vaidya_orizzonti` (event vs apparent horizon), `fig_vaidya_traiettorie` (τ-brachistochrones + teleological `p_v`), `fig_vaidya_variazionale`, `fig_vaidya_kerr_a0` (`a→0` residuals) |
| `vaidya_brachistochrone_vparam.py` | `fig_vaidya_bounce` (periapsis bounce, `v`-parametrization), `fig_vaidya_timing` (`r_min` vs arrival epoch) |
| `vaidya_penetration_map.py` | `fig_vaidya_penetration_map` (capture threshold `J_c(v₀)`, accretion/evaporation asymmetry) |
| `kodama_conservation.py` | `fig_kodama_conservazione` (`−u·K = Ê` to 1e−16 along the rail) |
| `plunge_vaidya_t_tau.py` | `fig_vaidya_plunge_t_tau` (plunge-radius law, `t` deeper) |
| `inversione_fisica.py` | `fig_vaidya_no_inversione_evaporazione` (no physical evaporative inversion) |
| `verifica_minimo_brachi.py` | `fig_verifica_minimo_brachi` (direct minimum-principle check) |
| `vaidya_hamiltonians_sympy.py` | symbolic `H_v`, `H_τ` from the ingoing-wind indicatrix |

### `ThakurtaMetric/` — conformal rotating compact object
| Script | Produces / does |
|---|---|
| `genera_figure_thakurta.py` | `fig_thakurta_cattura` (capture by expansion), `fig_thakurta_kerr_superfici` (rigid horizon vs breathing freezing surface), `fig_thakurta_kerr_residui`, `fig_indicatrici` |
| `genera_figure_thakurta_rami.py` | `fig_thakurta_kerr_rami` (the three branches from one ellipse) |
| `thakurta_kerr_plunge_t_tau.py` | `fig_thakurta_kerr_plunge_t_tau` (rotational plunge-inversion locus in `(a,J)`) |
| `inversione_conformale_AJ.py` | `fig_thakurta_kerr_inversione_AJ` (conformal inversion; `cond(r,A)=0` existence theorem) |
| `thakurta_kerr_plunge_inversion.py` | `fig_thakurta_kerr_plunge_map` (capture/escape map) |
| `cuspide_ergosfera.py` | `fig_thakurta_cuspide_ergosfera` (ergosphere cusp, power law 2/3 vs 2) |
| `tricotomia_figura.py` | `fig_thakurta_tricotomia_Jneg` (retrograde `J<0` capture band) |
| `thakurta_kerr_quasicostanti.py`, `genera_drift_map.py`, `genera_drift_map_A.py` | quasi-constant transfer and drift maps `(θ,r)` and `(A,r)` |
| `colormaps_inversione_aJ_AJ.py` | `fig_colormaps_inversione` (same-launch inversion in `(J,a)` and `(J,A)`) |
| `thakurta_kerr_sympy.py`, `thakurta_kerr_optical_sympy.py`, `thakurta_kerr_dynamic_K.py`, `thakurta_kerr_h_counterterm.py` | symbolic Hamiltonians, non-autonomous Randers reduction, 3D dynamic `K` and counterterm `S₁` |

### `KerrMetric/` — equatorial closed forms and fixed-endpoint studies
| Script | Produces / does |
|---|---|
| `kerr_separatrix_trajectories.py` | `fig_separatrix_3traiettorie` (separatrix: ODE vs quadrature vs evaluated Weierstrass, through the ergosphere) |
| `kerr_separatrix_gallery.py` | `fig_separatrix_gallery` (three-method superposition, four parameter sets) |
| `kerr_separatrix_weierstrass.py`, `kerr_separatrix_validation.py` | Weierstrass closed form `φ(r)` and its symbolic/numeric validation |
| `kerr_doran_shift_weierstrass.py` | Doran-continuation shift (second elliptic curve) |
| `bvp_estremi_fissi.py` | `fig_bvp_estremi_fissi` (fixed-endpoint `t/τ`, symmetric and asymmetric) |
| `bvp_kerr_simmetrico.py`, `bvp_conforme_inversione.py` | `fig_bvp_kerr_inversione`, `fig_bvp_conforme_inversione` (no fixed-endpoint inversion vs spin / conformal factor) |
| `colormap_spin_estremi_fissi.py`, `colormap_conforme_estremi_fissi.py` | `fig_colormap_spin_estremi_fissi`, `fig_colormap_conforme_estremi_fissi` (`Δr>0` maps, symmetric endpoints) |
| `colormap_spin_asimmetrico.py`, `colormap_conforme_asimmetrico.py` | `fig_colormap_spin_asimmetrico`, `fig_colormap_conforme_asimmetrico` (asymmetric-endpoint robustness maps) |

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
