# TOVInternalMetric

Numerical programs for coordinate-time and proper-time brachistochrones inside
Tolman-Oppenheimer-Volkoff stellar interiors.

The programs use geometrized units `G = c = 1`.  The preferred workflow for
realistic runs is CompOSE cold neutron-star data:

```text
eos.nb.ns
eos.thermo.ns
```

The importer reads the baryon-density grid, pressure, and energy density, then
converts `MeV fm^-3` to `km^-2` for the TOV integration.

Generic CSV tables are also supported when pressure and energy density have
already been converted consistently:

```text
pressure = pressure(density)
epsilon = epsilon(density)
```

The analytic relativistic polytrope is still available as a toy model,

```text
p = K rho^Gamma
eps = rho + p / (Gamma - 1)
```

The theory reference is `Versioni/TOVBrachistocrone.pdf`.

## Requirements

```bash
python3 -m pip install -r requirements.txt
```

The scripts use `numpy` for numerics and `reportlab` for APS-friendly vector
PDF export.  SVG files and CSV source data are also produced.

## Generate All Figures

```bash
python3 genera_grafici_tov.py
```

Output is written to `figures/`:

- `fig_01_tov_profiles`: TOV pressure, energy density, mass, metric functions and release kinematics.
- `fig_02_effective_indices`: coordinate-time and proper-time effective indices.
- `fig_03_brachistochrone_curves`: spatial curves for selected apertures.
- `fig_04_turning_radius_vs_delta`: turning radius as a function of aperture.
- `fig_05_travel_times_vs_delta`: coordinate and proper travel times.
- `fig_06_compactness_scan_*`: dependence on compactness along a polytropic sequence.

Each figure has a matching CSV file under `figures/data/`.

## Useful Examples

Run with a tabulated EOS:

```bash
python3 genera_grafici_tov.py \
  --eos compose \
  --compose-dir eos_tables/compose/qmc_rmf_1 \
  --eos-name qmc_rmf_1 \
  --central-density 0.8 \
  --sequence-density 0.25,0.35,0.45,0.6,0.8,1.0,1.2 \
  --outdir figures_compose/qmc_rmf_1
```

Generate only the curve comparison:

```bash
python3 genera_grafici_tov.py --figure curves --curve-deltas 45,90,135
```

Use a different central pressure:

```bash
python3 genera_grafici_tov.py --central-pressure 8e-4
```

Scan a different polytropic sequence:

```bash
python3 genera_grafici_tov.py --figure compactness --sequence-pc 5e-5,1e-4,2e-4,4e-4,8e-4,1.6e-3
```

See `eos_tables/README.md` for generic CSV tables and `eos_tables/compose/qmc_rmf_1/SOURCE.md`
for the included CompOSE dataset.
