# SageMath environment for the genus-two scripts

Only four scripts cited in Paper II need SageMath; everything else runs on the Python
environment of `requirements-lock.txt`. This folder makes that Sage environment
reproducible.

| File | Purpose |
|---|---|
| `setup_sage_env.sh` | Creates a conda environment with SageMath 10.9 (conda-forge) and installs `abelfunctions` at a fixed revision with the patch below. |
| `abelfunctions_sage10.patch` | Compatibility patch for `abelfunctions` under Sage 10.x (40 lines, two files). |
| `run_sage_checks.py` | Runs the four scripts with the right interpreter and checks the values quoted in Paper II against stated tolerances; exits non-zero on any failure. |
| `test_run_sage_checks.py` | Adversarial test of the runner on simulated outputs (no Sage needed): malformed, non-finite, missing or out-of-range values, a wrong number of eigenvalues, a failing or hanging script must all fail. |

## Reference environment

- SageMath **10.9**, conda-forge build `sage 10.9 py313hc25858d_0`, Python 3.13 (conda-forge), on macOS
  arm64 (native: no Rosetta, which is no longer available on macOS 27). The conda-forge packages also
  exist for Linux x86_64/aarch64; the recipe has been run on macOS arm64 only.
- `abelfunctions` from <https://github.com/abelfunctions/abelfunctions>, revision
  `d98843e4d25cebe9dc2993c86f469776142c5820` (MIT licence), plus `abelfunctions_sage10.patch`:
  Sage 10.x removed `is_LaurentSeries`, `is_LaurentSeriesRing`, `is_PowerSeries` and
  `is_PowerSeriesRing`; the patch replaces them with `isinstance` checks and drops an unused import.
- Build with **Cython 3.0** (`>=3.0,<3.1`): with Cython 3.1 or later the generated C code of this
  package does not compile (undeclared `__pyx_tp_new_vectorcall_*`). The setup script takes Cython 3.0
  from a temporary directory, so the environment keeps its own Cython.

## Build and check

```bash
cd NonStationaryMetrics/sage_env
./setup_sage_env.sh sage-nsb                # about 4-5 GB; a few minutes with a warm conda cache
python3 run_sage_checks.py --env sage-nsb   # about one minute
python3 test_run_sage_checks.py             # seconds; no Sage needed
```

## The four scripts

| Script | Interpreter | Needs `abelfunctions` | Quantity checked (Paper II) |
|---|---|---|---|
| `NonStationaryMetrics/KerrSessionScripts/kerr_tau_general_genus2_periods.sage.py` | environment Python | no | Riemann symmetry of the period matrix, about 7e-12; Im τ positive definite |
| `KerrScripts/kerr_tbranch_genus2_periods.sage.py` | environment Python | no | Riemann symmetry, about 3e-12; Im τ positive definite |
| `KerrScripts/kerr_quasiperiods_bel.sage` | `sage` | no | symmetry of the second-kind quasi-period matrix κ, about 1e-12 |
| `KerrScripts/kerr_thirdkind_theta_closed.sage` | `sage` | yes | third-kind θ closed form, relative residual about 7e-5 |

`.sage` files are run with `sage file.sage` from their own directory; the `.sage.py` files are already
preparsed and are run with the environment's `python` (conda's `sage` has no `-python` option).
`abelfunctions` must be imported inside Sage and not from its source directory, which would shadow the
installed package.

The tolerances in `run_sage_checks.py` are acceptance margins around the orders of magnitude quoted in
the paper (1e-11 for the symmetries, 1e-4 for the residual), to absorb run-to-run and platform variability;
they are not the printed values themselves, and they do not certify every digit. By default each value is read from a versioned JSON result with exact metric and sample IDs,
independently of stdout; the eigenvalues of Im τ must be exactly two (genus two), finite and
positive. The runner checks only the six quantities in the table above, not every other property a script
may print (for instance the full Legendre relation).


## Structured results

The four producers optionally export finite binary64 measurements through
`NSB_RESULT_JSON`; the runner uses this interface by default and applies its own
bounds. `--report-json /tmp/nsb-sage.json` writes the measurements and verdicts.
`--format legacy` retains the former text parser for migration comparisons;
`test_run_sage_checks.py` explicitly tests that legacy mode. The shared adversarial
JSON tests are `../tests/test_verification_results.py` (standard library only).
See [VERIFICATION_JSON.md](../VERIFICATION_JSON.md) for the schema and failure rules.
