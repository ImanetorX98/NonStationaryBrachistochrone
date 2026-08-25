# CompOSE QMC-RMF 1

Source: CompOSE EoS entry `297`, QMC-RMF tables by Brodie et al.

- EoS page: <https://compose.obspm.fr/eos/297>
- Density grid: <https://compose.obspm.fr/download//3D/Brodie/Qmc_rmf_1/eos.nb.ns>
- Thermodynamic table: <https://compose.obspm.fr/download//3D/Brodie/Qmc_rmf_1/eos.thermo.ns>

The imported cold neutron-star files are:

- `eos.nb.ns`: baryon-density grid `n_b` in `fm^-3`.
- `eos.thermo.ns`: cold beta-equilibrated thermodynamic quantities.

The TOV importer computes:

```text
p = n_b Q1
e = n_b m_n (1 + Q7)
```

and converts both from `MeV fm^-3` to `km^-2` for geometrized TOV integration.

