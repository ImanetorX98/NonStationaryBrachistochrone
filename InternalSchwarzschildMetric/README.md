# CoordinateTimeGravityTrain

Numerical scripts for gravity-train and brachistochrone trajectories inside a
uniform-density sphere, including Newtonian paths and relativistic internal
Schwarzschild models for coordinate time and proper time.

The repository contains the source code used to generate the figures and
animations. Generated images, videos, cached files, and local working notes are
excluded from version control.

## Requirements

Python 3.10+ with:

- numpy
- matplotlib
- Pillow

Install the Python dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

For MP4 animations, `ffmpeg` must also be available on the system path.

## Main Scripts

- `simula_caduta_tunnel.py`: Newtonian gravity-train trajectories and optional plots or animations.
- `simula_brachistocrona_schwarzschild_interna.py`: coordinate-time brachistochrone in the internal Schwarzschild metric.
- `simula_brachistocrona_tempo_proprio_schwarzschild_interna.py`: proper-time comparison under fixed-energy release from the surface.
- `genera_figure_paper.py`: batch generation for the main series, validation, error, and trajectory figures.
- `fig_brachistocrone_t_antipodi_numerico.py`: antipodal coordinate-time brachistochrone family.
- `fig_famiglia_brachistocrone_delta180_matplotlib.py`: publication-style antipodal family figure with Matplotlib.
- `fig_nlo_delta90_semicerchio.py`: NLO trajectory family for a 90 degree separation.
- `fig_nlo_delta60_primo_quadrante.py`: NLO trajectory family for a 60 degree separation.
- `fig_indice_rifrazione_nr_vs_mu.py`: effective refractive-index curve versus compactness.
- `animazione_confronto_brachistocrone.py`: comparison animation between Newtonian and relativistic paths.

## Examples

```bash
python3 simula_caduta_tunnel.py --delta-deg 90 --plot
python3 simula_brachistocrona_schwarzschild_interna.py --delta-deg 120 --u 0.2 --plot
python3 genera_figure_paper.py
python3 fig_brachistocrone_t_antipodi_numerico.py
```

Outputs are written either next to the script or under `figures/`, depending on
the command-line options.
