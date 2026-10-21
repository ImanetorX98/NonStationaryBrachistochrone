# NonStationaryBrachistochrone

Reproducibility package for two companion papers on **controlled-rail brachistochrones**:
the fastest timelike worldline from a fixed departure event to a fixed spatial target,
with free arrival time, when a controlled thrust holds the traveller at a maintained energy, in stationary and **non-stationary** spacetimes (FLRW,
Vaidya, and the conformal Kerr geometry of Thakurta).

| | Title | Concept DOI |
|---|---|---|
| **Paper I** | *Controlled-rail brachistochrones in non-stationary spacetimes: conformal symmetry, Kodama energy, and Vaidya dynamics* | [10.5281/zenodo.21739997](https://doi.org/10.5281/zenodo.21739997) |
| **Paper II** | *Brachistochrones in conformal Kerr spacetimes: control domain, exterior separatrices, and adiabatic response* | [10.5281/zenodo.21739999](https://doi.org/10.5281/zenodo.21739999) |
| **Code** | this repository | [10.5281/zenodo.21707377](https://doi.org/10.5281/zenodo.21707377) |

These are Zenodo concept DOIs: each always resolves to the latest version, and the
Zenodo record lists the DOI of every individual version.

## Where to start

- [`NonStationaryMetrics/README.md`](NonStationaryMetrics/README.md) — the project, its
  layout, installation, and how to compile both papers.
- [`NonStationaryMetrics/paper1/verification/README.md`](NonStationaryMetrics/paper1/verification/README.md)
  and [`NonStationaryMetrics/paper2/verification/README.md`](NonStationaryMetrics/paper2/verification/README.md)
  — claim-by-claim verification maps: where each result appears in the manuscript and
  which script reproduces it.
- [`NonStationaryMetrics/Fonti/README.md`](NonStationaryMetrics/Fonti/README.md) — the
  bibliographic audit.

## Layout

| Path | Contents |
|---|---|
| `NonStationaryMetrics/` | The non-stationary study: both manuscripts, all verification and figure scripts. |
| `KerrScripts/` | Stationary-Kerr base used by Paper II: separatrix closed forms (Weierstrass), genus-2 periods and theta functions (SageMath), fixed-endpoint boundary-value problems, colour maps and figures. |

Generated binaries (figures, PDFs, raw data) are not tracked; the scripts regenerate them.

## Licence and citation

The code is released under the MIT licence ([`LICENSE`](LICENSE)). Citation metadata
are in [`CITATION.cff`](CITATION.cff); please cite the code together with the companion
papers, through the DOIs listed in [`NonStationaryMetrics/README.md`](NonStationaryMetrics/README.md).
