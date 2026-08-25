# EOS Tables

Put tabulated equations of state here when running the TOV brachistochrone
programs with `--eos table`.

## Required Format

CSV with a header row.  The default column names are:

```text
density,pressure,epsilon
```

- `density`: density variable used by the table.
- `pressure`: pressure as a function of density.
- `epsilon`: total energy density used in the TOV equations.

If `epsilon` is omitted, the code interprets `density` as total energy density.
For rest-mass density tables, include an explicit `epsilon` column.

Alternative column names are accepted:

- density: `density`, `rho`, `rho0`, `rest_mass_density`, `energy_density`, `epsilon`, `eps`
- pressure: `pressure`, `p`, `press`
- energy density: `epsilon`, `eps`, `energy_density`, `total_energy_density`, `e`

Explicit names can also be passed on the command line:

```bash
python3 ../genera_grafici_tov.py \
  --eos table \
  --eos-table path/to/table.csv \
  --density-column rho \
  --pressure-column p \
  --epsilon-column eps \
  --central-density 0.002
```

## Included Sample

`sample_polytrope_table.csv` is a small tabulated version of the analytic
polytrope used for code validation.  It is not a realistic neutron-star EOS.

