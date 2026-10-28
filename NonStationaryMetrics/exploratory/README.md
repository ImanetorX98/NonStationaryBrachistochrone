# exploratory/ — superseded and exploratory scripts

This folder keeps, for the record, scripts that were written while the results of
Papers I and II were being developed and that are no longer part of any
verification route. They were moved here in v1.8.4 so that the main folders
contain only the scripts that reproduce the published numbers.

- None of these files is cited by either manuscript.
- None is listed in `paper2/provenance/MANIFEST.tsv`.
- None is run by `run_regression.py`, `sage_env/run_sage_checks.py` or the tests.

Some still run and print what they printed when they were written; others
document an approach that failed, and say so in their header. Nothing here is
needed to reproduce a published result. Each subfolder mirrors the folder the
file came from, so `exploratory/VaidyaMetric/` holds former `VaidyaMetric/`
scripts. Run them from `NonStationaryMetrics/`, like every other script in the
repository.

The comments and printed labels were translated into English in v1.8.4; the
code was not changed, apart from the import paths of three files (marked below)
that had to follow the move. Some local variable names are still Italian words
(for example `nome`, `curva`, `esito`); they are identifiers only.

## Top level (from `NonStationaryMetrics/`)

| File | What it did | Status / current route |
|---|---|---|
| `brachistochrone_thrust_fuel.py` | Proper thrust `|a|`, characteristic velocity and Tsiolkovsky fuel fraction along a Schwarzschild τ-branch brachistochrone. | Exploratory estimate, not used by the papers. The thrust figures of Paper II come from `paper2/verification/observational/` (`probe_saturn_comparison.py`, `probe_solar_sail.py`, `probe_control_law.py`). |
| `no_inversion_schwarzschild_frozen.py` | Frozen-case reduction of the τ half-angle map Φ(r_min): two quadratures of Φ_τ compared, plus a comment block on the monotonicity criterion. | Superseded by `no_inversion_reduction.py`, `no_inversion_schwarzschild_closedform.py` and the computer-assisted proofs `no_inversion_schwarzschild_CAP*.py`. Its comments cite two scripts that are not distributed. |
| `no_inversion_schwarzschild_rigor.py` | Attempt to make the large-`r0` asymptotic of Φ_τ rigorous (remainder bounds, localisation of the critical point). | Superseded by `no_inversion_schwarzschild_asymptotic.py` and the CAP certificates; the localisation step only worked for `r0 ≳ 300`. |

## `FLRWmetric/`

| File | What it did | Status / current route |
|---|---|---|
| `flrw_perturbed_genus_jump.py` | Proof of concept: perturbing the FLRW index `n(r)` raises the genus of the spectral curve of the brachistochrone. | Kept as an illustration; not used by the papers. |

## `KerrSessionScripts/`

| File | What it did | Status / current route |
|---|---|---|
| `costate_finsler_derivation.py` | Legendre map of the rail brachistochrone: Pontryagin costate against the mechanical momentum. | Deprecated diagnostic (header explains why). The costate statement of the papers is checked by `VaidyaMetric/costate_vs_mechanical.py`. |
| `costate_vs_mechanical_L.py` | Tested whether the conserved costate `p_φ = J` equals the mechanical `u_φ`; the test fails, as expected. | Deprecated diagnostic; same current route as above. |
| `kerr_tau_Wij_naming.sage` (+ `.sage.py`) | Named the second-kind integrals `R_i` as bare Kleinian `ζ_i(u)`. | Failed at the θ-divisor (documented). Replaced by the difference form in `KerrSessionScripts/kerr_tau_Wij_diffform.sage` and `kerr_tau_Wij_diffform_integral.sage`. |
| `kerr_tau_Wij_naming2.sage` (+ `.sage.py`) | Second attempt: exact Riemann constant, Abel map based at the branch point. | Same as above. |
| `kerr_tau_Wij_naming3.sage` (+ `.sage.py`) | Third attempt: fit of a constant Jacobian offset `C` in `ζ_i(u + C)`. | Same as above. |
| `kerr_tau_dilog_qseries1.sage` (+ `.sage.py`) | q-series step 1: nomes from the Riemann matrix, θ nome series against `RiemannTheta`, ground-truth dilogarithm. | Intermediate step; the result is `KerrSessionScripts/kerr_tau_dilog_qseries5.sage`. |
| `kerr_tau_dilog_qseries3.sage` (+ `.sage.py`) | q-series step 3: Ω = log θ-ratio as a convergent nome series after the Fricke transformation. | Intermediate step; see `kerr_tau_dilog_qseries5.sage`. |
| `kerr_tau_dilog_qseries4.sage` (+ `.sage.py`) | q-series step 4: the dilogarithm Λ computed from the nome series of Ω. | Intermediate step; see `kerr_tau_dilog_qseries5.sage`. |
| `offshell_taubranch_closed_form.py` | Special-function closed form of the Thakurta–Kerr τ-branch off-shell dilation wrap. | Development of the off-shell closure; the closed form itself has no maintained equivalent. For the τ branch the kept route is the reduction in `ThakurtaMetric/adiabatic_offshell_reduction.py` (it treats both branches) and the validation against the true canonical flow in `ThakurtaMetric/fig_phi_validation_corrected.py` and `tests/test_adiabatic_noreg.py`. The symbolic coefficient table `adiabatic_offshell_coefficients.py` and the exact check `adiabatic_first_order_exact.py` are for the t/η branch. |
| `offshell_taubranch_physics_anchor.py` | Checked the closed form above against the true non-autonomous flow. **Import path changed** (`parents[2]`, package path `exploratory.KerrSessionScripts`) to follow the move. | As above: τ-branch validation in `fig_phi_validation_corrected.py` and `tests/test_adiabatic_noreg.py`. |
| `offshell_tbranch_closed_form.py`, `offshell_tbranch_wrap_assembly.py`, `offshell_tbranch_kernel_hermite.py`, `offshell_tbranch_thirdkind_Delta.py`, `offshell_tbranch_block3_assembly.py`, `offshell_tbranch_FULL_assembly.py` | t/η-branch off-shell wrap closed block by block: second-kind block in the canonical genus-2 basis, Hermite reduction of the kernel pole, third-kind letters at Δ = 0, full assembly against the direct double integral. | Development of the t/η closure. The kept route is `ThakurtaMetric/adiabatic_offshell_reduction.py`, `adiabatic_offshell_coefficients.py` (symbolic coefficients; the genus-2 quadratures are evaluated numerically there, not in θ functions) and `adiabatic_first_order_exact.py`; the t-branch on-shell assembly is `../KerrScripts/kerr_tbranch_psi_assembly.py`. |
| `physics_anchor_offshell_closed.py` | Chain of checks from the t-branch closed form to the true non-autonomous flow (slope 2). | As above. |

## `ThakurtaMetric/`

| File | What it did | Status / current route |
|---|---|---|
| `adiabatic_offshell_validation.py` | Early exploration with a *surrogate* flow (no dilation term), compared with half the on-shell coefficient. | **Withdrawn reading** (header explains). The validation against the true canonical flow is `ThakurtaMetric/fig_phi_validation_corrected.py`. |
| `adiabatic_offshell_kernel_wrap.py` | Showed that the off-shell kernel on the frozen shell is `A(r)/√S`, on the same genus-2 curve as the on-shell reduction. | Folded into `ThakurtaMetric/adiabatic_offshell_reduction.py`. |
| `adiabatic_offshell_wrap_assembly.py` | Weight-two assembly of the off-shell term into the genus-2 `W_jk` basis. | As above. |
| `adiabatic_offshell_physics_anchor.py` | Methodological anchor of the off-shell reduction against the physical flow. | As above; the slopes printed in Paper II come from `adiabatic_first_order_exact.py` through `paper2/provenance/`. |
| `colormaps_inversion_aJ_AJ.py` | Two colormaps of `Δr = r_min^t − r_min^τ` on the `(J, a)` and `(J, A)` planes. **Import path changed** (`sys.path` one level higher) to follow the move. | Not used by the papers; the inversion results are in `ThakurtaMetric/conformal_inversion_AJ.py` and `thakurta_kerr_inversion_t_tau.py`. |
| `psi_decomp_launch.py` | Split ψ into a second-kind (ζ) and a third-kind (dilogarithm) piece for several launch radii. | Intermediate check; the decomposition used is in `../KerrScripts/kerr_psi_explicit_verified.py` and `KerrSessionScripts/kerr_tau_Wij_assembly.py`. |
| `separatrix_offshell_audit.py` | Audit of the two separatrix degenerations (ergospheric τ, moving double root of the t branch). | Exploratory; the separatrix results are in `ThakurtaMetric/tk_sep_*.py`. |
| `separatrix_canonical_cubic.py` | Canonical cubic-pole term for the terminally anchored t/η separatrix problem. | Exploratory; as above. |

## `VaidyaMetric/`

| File | What it did | Status / current route |
|---|---|---|
| `evaporation_inversion.py` | Claimed an inversion for evaporating Vaidya. **Import path changed** (`sys.path` one level higher, for `paper_style`) to follow the move. | **Deprecated and wrong**: the inversion falls where the mass is already negative. Refuted by `VaidyaMetric/physical_inversion.py`; see `VaidyaMetric/VaidyaResults.md`. |
| `vaidya_offshell_tau_clean.py`, `vaidya_offshell_shellpoly.py`, `vaidya_offshell_both_numeric.py`, `vaidya_offshell_closed_form.py`, `vaidya_offshell_FULL_assembly.py` | Steps of the τ-branch (and early v-branch) off-shell closure: clean Hamiltonian, shell polynomials, numeric wrap, closed form, full assembly. | Development steps. The off-shell results used by Paper I come from `VaidyaMetric/vaidya_first_order_offshell.py`, `vaidya_tau_assembly.py` and `vaidya_fully_explicit.py`. |
| `vaidya_tau_physics_anchor.py` | Checked the τ-branch closed form against the true non-autonomous flow. | As above. |
| `vaidya_v_offshell_step1.py`, `vaidya_v_step1_clean.py`, `vaidya_v_kernel_probe.py`, `vaidya_v_offshell_structure.py`, `vaidya_v_offshell_assembly.py`, `vaidya_v_explicit_assembly.py`, `vaidya_v_FULL_closure.py`, `vaidya_v_physics_anchor.py` | The same programme for the advanced-time (v) branch: spectral curve, kernel and inner-letter split, assembly, physics anchor. | As above; the v-branch result of Paper I is `VaidyaMetric/vaidya_first_order_offshell.py`. |
| `vaidya_separatrix_Pih.py` | Second-kind horizon building block on the separatrix in Weierstrass ζ. | Folded into the separatrix closed forms (`VaidyaMetric/vaidya_sep_*.py`). |
| `vaidya_sep_rank_hp.py`, `vaidya_sep_beta_hp.py` | High-precision (mpmath) rank of the dilogarithm system and the coefficients β. `vaidya_sep_beta_hp.py` stops with an `IndexError` with the pinned library versions; it did so before the move as well. | Superseded by `VaidyaMetric/vaidya_sep_fay_rank.py`, `vaidya_sep_fay_fast.py` and `vaidya_sep_5term.py`. |
| `vaidya_sep_fay_analytic.wl` | Mathematica attempt at an analytic Fay reduction through the Weierstrass addition theorems. | Superseded by `VaidyaMetric/vaidya_sep_fay_symbolic.wl`. |
| `vaidya_thirdkind_fay.sage` (+ `.sage.py`) | Fay closure of the horizon third-kind letter as a θ-ratio. Needs `abelfunctions`. | Exploratory; the third-kind θ closure checked by `sage_env/run_sage_checks.py` is `../KerrScripts/kerr_thirdkind_theta_closed.sage`. |

## `PaperFigures/`

| File | What it did | Status / current route |
|---|---|---|
| `make_fig_brachistochrone_comparison.py` (was `PaperFigures/genera_fig_confronto.py`) | Side-by-side figure of the brachistochrones in FLRW, Vaidya and Thakurta–Kerr (`fig_brachistochrone_comparison`). | The figure is no longer in either manuscript. It writes into this folder. |

## `observational/` (from `paper2/verification/observational/`)

| File | What it did | Status / current route |
|---|---|---|
| `probe_kerr_m87.py`, `probe_kerr_m87.json` | Kerr thought experiment scaled to M87* (rails at 120M and 38M). | Removed from Paper II on 24 October 2026; the thought experiment is now around Gaia BH3 (`paper2/verification/observational/probe_gaia_bh3.py`). Kept because the archived releases up to v1.8.3 cite it. |
