# Constrained-worldline brachistochrones in non-stationary spacetimes

Reproducibility package for two companion papers on the fastest constrained **"rail"
worldline** from a fixed departure event to a fixed spatial target, with free arrival
time: a timelike curve held at a maintained energy by a controlled thrust, whose arrival
time (or proper time) is minimised. The study extends
the stationary-Kerr optical and genus-two machinery to **non-stationary** spacetimes:
FLRW, Vaidya (a dynamical black hole) and Thakurta–Kerr (the conformal Kerr geometry
`g = A(eta)^2 g_Kerr`).

The closed forms come with scripts that derive them symbolically and check them
numerically, and many results are cross-checked independently in Mathematica. The
evidence is not uniform: some checks are exact symbolic identities, some are numerical
or sampled, some use a single computer-algebra system, and a few long certificates are
archived rather than re-run by default. Each paper has its own verification map
(below) stating, claim by claim, where a result appears, which script reproduces it,
and at what level.

> Generated binaries (`*.png`, `*.pdf`, videos, raw `.npz` data) are **not** tracked;
> regenerate them from the scripts.

> **Structure.** This folder holds the non-stationary study and both manuscripts. The
> stationary-Kerr scripts the papers use (`KerrScripts/`) sit one level up, as a sibling
> (`../KerrScripts/`); the genus-two closure scripts used here
> are copied under `KerrSessionScripts/`.

---

## The two papers and the code

| | Title | DOI (always resolves to the latest version) |
|---|---|---|
| **Paper I** | *Controlled-rail brachistochrones in non-stationary spacetimes: conformal symmetry, Kodama energy, and Vaidya dynamics* | [10.5281/zenodo.21739997](https://doi.org/10.5281/zenodo.21739997) |
| **Paper II** | *Brachistochrones in conformal Kerr spacetimes: control domain, exterior separatrices, and adiabatic response* | [10.5281/zenodo.21739999](https://doi.org/10.5281/zenodo.21739999) |
| **Code** | this repository | [10.5281/zenodo.21707377](https://doi.org/10.5281/zenodo.21707377) |

These are Zenodo concept DOIs. Each deposited version also has its own DOI, listed on
the Zenodo record; cite a version DOI to fix the exact text or code you used (the
manuscripts cite the code release that reproduces their printed numbers).

- **Paper I** (`paper1/paper1_JMP.tex`) is the spherically symmetric sector: the
  controlled-rail formalism as a time-optimal control problem (existence, normality of
  interior minimizers, Hamilton–Jacobi verification), FLRW, Kodama energy and Vaidya.
- **Paper II** (`paper2/paper2.tex`) is the axisymmetric sector, Thakurta–Kerr. It takes
  its control-theory foundations from Paper I and adds: the compact-control domain and
  its end at the conformal stationary limit `r = 2M`; a normality lemma for interior
  free-arrival minimisers; equatorial separatrices in closed Weierstrass form, including
  an exterior retrograde one; the first-order adiabatic response, validated against the
  true non-autonomous flow; the classification of rail-swept surfaces (umbilicity and
  geodesicity under a conformal factor); and an observational section with physical
  scales, rotating benchmarks, a proposed Solar-System two-probe comparison (Saturn
  flyby), a Kerr thought experiment scaled to M87*, an ephemeris bound on a locally
  cosmological conformal factor, and the conjugate-point conditions under which a rail
  stops being a local minimum.

The combined master `paper/main.tex` (iopart) and `paper/main_prd_revtex.tex` (PRD) are
the **historical** single-paper version, kept for reference; the two companion papers
supersede them.

### Which scripts belong to which paper

| Paper | Directories / representative scripts |
|---|---|
| **Paper I** (rail / Kodama / Vaidya, spherical) | `FLRWmetric/`, `VaidyaMetric/`, `paper1/verification/`; shared control-theory checks `KerrSessionScripts/perlick_equivalence*.py`; `VaidyaMetric/costate_vs_mechanical.py`, `VaidyaMetric/universal_SD_source_vaidya.py` |
| **Paper II** (Thakurta–Kerr, axisymmetric) | `ThakurtaMetric/`, `KerrSessionScripts/` (genus-2 closure), `paper2/verification/` (including `observational/`), `paper2/provenance/`, `../KerrScripts/` |
| **Shared** | `paper/` (bibliography, historical master, `crosscheck_*.wl`), `PaperFigures/`, `run_regression.py`, `requirements*.txt`, `paper_style.py` |

---

## Layout

| Path | Contents |
|------|----------|
| `paper1/` | Paper I: `paper1_JMP.tex` (the live manuscript), `build_paper1.sh` (figures + compilation from a clean checkout), `verification/` (symbolic verification suite), `submitted_JMP_2026-08/` (frozen copy of the version submitted to JMP). |
| `paper2/` | Paper II: `paper2.tex`; `verification/` (Mathematica and Python checks of every numbered statement, and `observational/` for the observational section and the probe experiments); `provenance/` (generated slopes, hashes and cross-references; see below). |
| `paper/` | Shared bibliography (`refs.bib`, `companionI.bib`, `companionII.bib`), the historical master, and the Mathematica cross-checks `crosscheck_*.wl`. |
| `FLRWmetric/` | FLRW (degenerate base) analysis and the perturbed-FLRW genus-jump proof of concept. |
| `VaidyaMetric/` | Vaidya (dynamical BH, parameter = mass `m`): analysis, closed forms, figures. |
| `ThakurtaMetric/` | Thakurta–Kerr (conformal rotating, parameter `A` via `E_eff = Ê/A`): analysis, closed forms, figures. |
| `KerrSessionScripts/` | Genus-2 machinery and the generic-J (genus-2) closure pipeline. Two files marked DEPRECATED are exploratory diagnostics kept for the record. |
| `Experimental/` | Laboratory design studies (proofs of concept, not measured data). |
| `Fonti/` | `README.md` only: the bibliographic audit (what each source says and where it is cited). The PDFs are not distributed. |
| `reproduce_reductions.py` | One-file check of the algebraic heart of the adiabatic reductions. |
| `paper_style.py` | Shared matplotlib style helper. |

### Reference documents (read these first)
| File | What it collects |
|------|------------------|
| `paper1/verification/README.md` | Verification map of Paper I: which claim, where it appears, how it is checked (symbolic identities in two CAS for the formal core, plus a numerical quadrature of the §4.4 closed form). |
| `paper2/verification/README.md` | Verification map of Paper II: each numbered statement, its script, and the evidence level (full proof, sampled, or numerical). |
| `paper2/verification/observational/README.md` | Observational section: scales, rotating benchmarks, probe experiments (Saturn flyby, solar sails, M87*), conjugate points, ephemeris bound. |
| `paper1/verification/REPRODUCTION_section4.2.md` | Re-run of the trajectory scripts behind section 4.2 of Paper I, with the numbers obtained. |
| `ThakurtaMetric/ThakurtaResults.md`, `VaidyaMetric/VaidyaResults.md`, `FLRWmetric/FLRWresults.md` | Per-geometry results: what was derived, with the numbers obtained. *Historical working notes, mostly in Italian; they may refer to scripts or sibling folders that are not distributed. The verification maps above are the current reference.* |
| `ThakurtaMetric/phi_adiabatic_closed_form.md`, `VaidyaMetric/vaidya_adiabatic_setup.md` | The adiabatic closed form, set up and assembled (historical working notes, as above). |
| `Fonti/README.md` | Bibliographic audit of both papers. |

Private working notes, editorial correspondence and separate ongoing research are
not included in this release and are not required by any of the verification routes
listed here.

### Provenance of the printed numbers (Paper II)

The adiabatic slopes, the validation window and the SHA-256 prefixes of the scripts
printed in Paper II are not typed by hand: `paper2/provenance/make_provenance.py`
generates them as LaTeX macros from one dataset, and `--check` regenerates everything
in a scratch directory and fails on any difference. The comparison is exact in the
reference environment of `requirements-lock.txt` (Python 3.13, NumPy 2); in other
environments floating-point differences of order 1e-16 in the raw data can change the
hashes while the results agree. On a fresh clone, generate the raw dataset once first:

```bash
python3 ThakurtaMetric/fig_phi_validation_corrected.py      # writes the raw .npz (~15 s)
python3 paper2/provenance/make_provenance.py --check        # regenerate and compare
python3 paper2/provenance/make_provenance.py --help         # options
```

`paper2/provenance/collect_figures.py --check` compares each compiled figure with its
generator's last output (it does not re-run the generators; two figures written directly
into `paper/Immagini` and one without a listed generator are reported, not checked).

### Fixed-endpoint no-inversion (Lemma B)

> **Status (tiered, not absolute).** The fixed-endpoint no-inversion is **proved**
> for `r0 <= R*(E)`; **machine-certified** (interval arithmetic) at representative `r0`;
> supported, in the static case, by a large-`r0` asymptotic analysis checked numerically;
> and **conjectural** in general. The static/frozen branch ordering does not invert
> under the outgoing-clock replacement; the full outgoing-Vaidya optimal-control
> boundary-value problem **remains open**.

| File | What it does |
|------|--------------|
| `no_inversion_reduction.py` | Reduction + Lemma A closed-form proof (`symbolic_full_proof`). |
| `no_inversion_schwarzschild_closedform.py` | Closed form for `Phi_tau'` (IBP); elementary grazing/quarter regimes; single-crossing criterion. |
| `no_inversion_schwarzschild_asymptotic.py` | Large-`r0` asymptotic: `Phi=arctan(...)-b^{-3/2}sqrt(...)`, single non-degenerate max, `Phi''(x_pk)=-1/2 b^{9/4} V0^{-5/4}`. |
| `no_inversion_schwarzschild_CAP.py`, `..._CAP_r0_10.py`, `..._CAP_grid.py` | Computer-assisted proof (mpmath interval arithmetic): single-crossing certificates. Complete at `r0=10`; grid `E in {1.2,1.6,2.5} x r0 in {8,12}` (`CAP_grid_certificates.txt`). |
| `verify_lemmaB_mathematica.wls` | Independent Mathematica cross-check of all no-inversion derivations. |
| `ThakurtaMetric/adiabatic_offshell_validation.py` | **Superseded, kept for the record.** Early exploration: it integrates a *surrogate* flow (radial momentum evolved with `-H_r` only, without the dilation term `-eps P_r`) and compares it with half the on-shell Euler coefficient. It is not a validation of the true canonical dynamics, and its old "~98% / ~2%" reading is withdrawn; use `fig_phi_validation_corrected.py` below. |
| `ThakurtaMetric/adiabatic_tk_geodesic_check.py` | Integrates the time-dependent Thakurta–Kerr geodesic equations and checks what drifts along the orbit (`E_eff`, `J_eff`). |

### Complete first-order extended-Hamiltonian correction (Eq. 40, on-shell + off-shell)
| File | What it does |
|------|--------------|
| `ThakurtaMetric/adiabatic_first_order_exact.py` | Exact canonical PT of the non-autonomous optical Hamiltonian: first-order coefficient `P_r^(1) = (s Theta Hbar - S_D)/Hbar_Pr`, with the complete source `S_D = int (Theta Hbar + P_r Hbar_Pr) du` and `Theta = E_eff d_E_eff + J_eff d_J_eff`; the increment is `eps P_r^(1)`; closes the **Thakurta–Kerr** first order to `O(eps^2)` (fitted residual slope about 2.15, robust fit 2.21, on the stated window). |
| `VaidyaMetric/vaidya_first_order_offshell.py` | **Vaidya** `v`-branch analogue (`Theta=m d_m`, sign flip for accreting `m`, terminal anchoring of `p_v`); closes to `O(eps^2)` (slope **2.00** full vs **1.00** on-shell only). Emits `fig_vaidya_offshell`. |
| `ThakurtaMetric/fig_phi_validation_corrected.py` | True-dynamics validation figure (`fig_phi_validation_true_dynamic`): on-shell coefficient alone, `-int lam Theta F dr` with `Theta F = Theta G - (G_Pr/H_Pr) Theta H` (checked against finite differences of the explicit frozen integrand), slope ~1, vs the complete Eq. (40) correction (slope about 2.12, approaching 2 as the epsilon window shrinks), both against the true non-autonomous flow. It also prints the size of the two parts: on the fit windows the on-shell coefficient exceeds the complete one by a factor of about 19 to 41 (sampled sup norms), the off-shell part `-int (G_Pr/H_Pr) S_D dr` cancelling most of it. |

---

## Install / run

```bash
python3 -m pip install -r requirements.txt        # numpy scipy sympy matplotlib mpmath Pillow
python3 -m pip install -r requirements-lock.txt   # exact reference environment
```

`requirements-lock.txt` pins the reference environment (Python 3.13, NumPy 2). The
figure generators need NumPy 2 (`np.trapezoid`); the verification suites were also
re-run in October 2026 with Python 3.11, NumPy 1.24 and SciPy 1.11.

- **Python (sympy/mpmath/numpy/scipy)** — algebraic reductions, genus-1 closed forms
  (Weierstrass via `mpmath.jtheta`), numerical checks, figures. Standalone:
  `python3 <script>.py`.
- **SageMath 10.9 + `abelfunctions`** — genus-2 objects (Riemann matrix,
  `RiemannTheta`, Abel map): the `.sage` scripts. Run inside `sage` from a neutral cwd.
  The environment is rebuilt by `sage_env/setup_sage_env.sh` (conda-forge Sage, fixed
  `abelfunctions` revision and patch), and `sage_env/run_sage_checks.py` checks the four
  Sage results cited in Paper II; see `sage_env/SAGE_SETUP.md`.
- **Mathematica / `wolframscript`** (13.3) — independent cross-checks
  (`paper/crosscheck_*.wl`, `paper1/verification/*.wls`, `paper2/verification/*.wls`).

Verification scripts print their residuals or PASS/FAIL lines; the newer ones also exit
with a non-zero status on failure. Read the printed lines in any case: an exit code
alone is not the evidence.

> **Run all commands from inside this folder** (`.../NonStationaryMetrics/`), so paths like
> `VaidyaMetric/…`, `KerrSessionScripts/…`, `paper/…` resolve. The stationary-Kerr base lives
> one level up, so reference it as `../KerrScripts/…`.

### Quick start — key commands (explicit paths, run from `NonStationaryMetrics/`)
```bash
# --- symbolic coefficients (Python, fast) ---
python3 SEP_COEFF_SYMBOLIC.py                       # separatrix b1,b2,b3 (all branches), (M,a,E,r_d,Jc)
python3 sep_tracking_coeff.py                       # tracking coeffs + triple-pole cancellation
python3 VaidyaMetric/vaidya_generic_coeff.py        # Vaidya generic c_k,Q_kj,g_i,P in (m,E,J)
python3 ThakurtaMetric/tk_t_generic_coeff.py        # TK-t generic c_k,Q_kj,g_i in (a,E,J)
python3 KerrSessionScripts/kerr_psi_explicit_verified.py     # TK-τ source reduction (c_k)
python3 KerrSessionScripts/kerr_tau_Wij_principalparts.py    # g_i (principal parts)
python3 KerrSessionScripts/kerr_tau_Talg_explicit.py         # T_alg explicit (poly + logs)
python3 VaidyaMetric/vaidya_sep_residui_analitici.py         # closed-form residues (foundational)
python3 VaidyaMetric/sep_v_clock_residui.py                  # v-branch clock residues (4M)
python3 ThakurtaMetric/sep_t_clock_residui.py               # t-branch clock residues (2M)

# --- genus-2 assembly / naming / q-series (SageMath) ---
sage KerrSessionScripts/kerr_tau_Wij_assembly.py            # NB: pure-sympy, also runs with python3
sage KerrSessionScripts/kerr_tau_Wij_diffform_integral.sage # naming U_k in θ[δ] at e_±
sage KerrSessionScripts/kerr_tau_dilog_qseries5.sage        # dilog q-series (geometric convergence)
sage VaidyaMetric/vaidya_dilog_qseries.sage                 # Vaidya dilog q-series
sage ThakurtaMetric/tk_t_dilog_qseries.sage                 # TK-t dilog q-series
sage KerrSessionScripts/kerr_quasiperiods_bel.sage          # genus-2 quasi-periods κ

# --- independent cross-checks (Mathematica) ---
wolframscript -file paper/crosscheck_generic.wl     # Vaidya+TK-t generic coeffs
wolframscript -file paper/crosscheck_sep_bi.wl      # separatrix b_i (τ) via Laurent
wolframscript -file paper/crosscheck_tkt_bi.wl      # separatrix b_i (TK t±)
wolframscript -file paper/crosscheck_tracking.wl    # tracking theorem N_tot(r_d)=0
wolframscript -file paper/crosscheck_clock_res.wl   # clock residues (native Weierstrass)

# --- stationary-Kerr base machinery (sibling at repo root) ---
python3 ../KerrScripts/pipeline_completa_deltaphi.py   # full base (M,a,E,J)→δφ chain
python3 ../KerrScripts/kerr_psi_explicit_verified.py   # base TK-τ reduction
sage    ../KerrScripts/kerr_quasiperiods_bel.sage      # base genus-2 quasi-periods
```

---

## Script guide (by pipeline stage)

Naming: `sep_` = separatrix (genus-1); `generic`/`Wij`/`dilog` = generic-J (genus-2);
`crosscheck_*.wl` = Mathematica. **★ = new in the closed-form / verification effort.**

### 0. Metric setup & frozen orbits
- `*/*_metric_sympy.py`, `*_optical_metric_sympy.py`, `*_hamiltonians_sympy.py` — build
  each metric, its optical (Randers) form, and the branch Hamiltonians.
- `vaidya_kodama_sympy.py`, `kodama_conservation.py` — Kodama energy = rail invariant.
- `KerrSessionScripts/kerr_tau_general_genus2.py` — **generic-J frozen orbit** φ₀ on the
  τ sextic (genus-2 Kleinian), symbolic 1st/3rd-kind coefficients.

### 1. Separatrix closed forms (genus-1, `|J|=Jc`)
Frozen orbit + first-order adiabatic δφ in Weierstrass σ,ζ,℘ + Brown–Levin Γ̃.
- `VaidyaMetric/vaidya_sep_residui_analitici.py` — **closed-form residues** b1,b2,b3 (source at the triple pole) + clock residues. Foundational.
- `vaidya_sep_G_partialfractions.py`, `vaidya_sep_deltaphi*.py`, `vaidya_sep_weight2_*.py` — pole-adapted reduction, δφ block assembly, weight-2 (elliptic dilog) closure.
- `vaidya_sep_C0Ce_closed.py` — additive constants Ce, C0 in closed form (ζ,℘ at marked points).
- `vaidya_sep_vbranch*.py`, `vaidya_sep_v_*assembly.py`, `vaidya_sep_v_track*.py` — advanced-time (v) branch + tracking.
- `vaidya_sep_track.py` — Jc-tracking (moving separatrix).
- `vaidya_sep_fay_*.py`, `vaidya_sep_compact*.py`, `vaidya_sep_5term.py` — Fay compression 63→5 dilogs.
- `brown_levin_gamma.py` — the Γ̃ elliptic dilogarithm + its q-series (Kronecker–Eisenstein).
- `ThakurtaMetric/tk_sep_*.py` — TK τ separatrix (baseline, block assembly, tracking).
- `ThakurtaMetric/tk_t_sep_*.py` — TK t separatrix, both `Jc±` (prograde/retrograde), incl. tracking.

**★ Reusable symbolic coefficients (this effort):**
- `SEP_COEFF_SYMBOLIC.py` — **b1,b2,b3 symbolic in (M,a,E,r_d,Jc) for ALL branches**, via S-derivatives at r_d. `VaidyaMetric/sep_coeff_symbolic.py` = Vaidya-only version.
- `VaidyaMetric/sep_v_clock_residui.py` — v-branch clock residues (E r_d³/s, horizon **4M**).
- `ThakurtaMetric/sep_t_clock_residui.py` — t-branch clock residues (horizon r±, invariant 2M).
- `sep_tracking_coeff.py` — **tracking coefficients**; theorem: tracking cancels the triple pole (b3=0 because N_tot(r_d)=0).
- `VaidyaMetric/sep_periodlevel_test.py` — period-level quantities change with parameters (not universal).

### 2. Generic-J closed forms (genus-2, off-separatrix) — `KerrSessionScripts/`
The `W_ij` closure pipeline for TK-τ (archetype), brick by brick:
- `kerr_psi_explicit_verified.py` — **★ source reduction** ∂_E F → Σ c_k r^k/√S (symbolic c_k).
- `kerr_tau_Wij_reduction.sage` — ★ U_0..U_4 independent; U_2 = unique 3rd-kind (dilog source).
- `kerr_tau_Wij_oddmodel_reduce.py` — ★ 2nd-kind reduction of U_k to canonical Abelian integrals (odd model), symbolic; verified diff/integral.
- `kerr_tau_Wij_assembly.py` — ★ assembles ψ = ½Ê Σ P_ab w_ab + T_alg; isolates the single dilog (end-to-end 1e-15).
- `kerr_tau_Wij_principalparts.py` — ★ symbolic pole coefficients g_i = Taylor of q6^{-1/2}.
- `kerr_tau_Talg_explicit.py` — ★ T_alg in explicit elementary form (polynomial + Σ log over roots of S).
- `kerr_tau_Wij_naming{,2,3}.sage` — ★ naming U_k in Kleinian ζ (θ-divisor wall documented).
- `kerr_tau_Wij_diffform.sage`, `_diffform_integral.sage` — ★ **difference-form naming** that avoids the θ-divisor wall (θ[δ] at e_±), verified 1e-6.
- `kerr_tau_Wij_aperiods.sage` — ★ holomorphic coefficients α_k,β_k from a-periods (period-level).
- `kerr_tau_Wij_holomorphic.sage` — ★ diagnostic: holomorphic part not orbit-extractable.
- `kerr_tau_siegel_reduce.py` — ★ Siegel reduction of τ (nomes <1 via S = −τ⁻¹).
- `kerr_tau_dilog_qseries{1,3,4,5}.sage` — ★ **dilog q-series**: θ nome series, split Q(elem)+L(dilog), geometric convergence (N=4→1e-13).
- `kerr_quasiperiods_bel.sage`, `kerr_psi_forward_abel.sage`, `kerr_thirdkind_theta_closed.sage`, `kerr_holo_component_check.sage` — genus-2 quasi-periods (κ), Abel map, 3rd-kind θ closure, Hodge component check.

**★ Other generic-J branches (this effort, same template):**
- `VaidyaMetric/vaidya_generic_coeff.py` — Vaidya generic (param m): c_k, Q_kj, g_i, P symbolic in (m,E,J), verified 1e-15.
- `VaidyaMetric/vaidya_dilog_qseries.sage` — Vaidya generic dilog q-series (geometric convergence).
- `ThakurtaMetric/tk_t_generic_coeff.py` — TK-t generic: c_k, Q_kj, g_i symbolic in (a,E,J).
- `ThakurtaMetric/tk_t_dilog_qseries.sage` — TK-t generic dilog q-series.

### 3. Independent cross-checks — `paper/crosscheck_*.wl` (Mathematica) ★
| File | Checks |
|------|--------|
| `crosscheck_genus2.wl`, `crosscheck_Ialg_symbolic.wl`, `crosscheck_fullsymbolic.wl`, `crosscheck_P_params.wl` | TK-τ generic coefficients (g_i, c_k, Q_kj, T_alg elementary), full symbolic in (a,E,J). |
| `crosscheck_generic.wl` | Vaidya + TK-t generic reductions (identity=0, g_i, c_k match Python). |
| `crosscheck_sep_bi.wl`, `crosscheck_tkt_bi.wl` | Separatrix b_i (Vaidya τ, TK τ, TK t±) via independent Laurent extraction (1e-16…1e-38). |
| `crosscheck_clock_res.wl` | Clock residues (v: 4M via native Weierstrass; t: formula). |
| `crosscheck_tracking.wl` | Tracking theorem: N_tot(r_d)=0 mod {S,S'} (triple pole cancelled). |
| `crosscheck_identities.wl`, `crosscheck_tk_numeric.wl`, `crosscheck_vt_numeric.wl`, `crosscheck_tkt_flat.wl` | Separatrix identities and v/t numeric (found+fixed the original t-branch bug). |

Run any of them: `wolframscript -file paper/crosscheck_<name>.wl`.

### 4. Phenomenology, inversion, figures
- `*_penetration_*.py`, `*_asymmetry.py`, `plunge_*_t_tau.py`, `inversione_*.py`,
  `thakurta_kerr_inversione_t_tau.py`, `kerr_adiabatic_phi_hybrid*.py` — penetration/scattering
  trichotomy, accretion↔evaporation asymmetry, plunge inversion, hybrid adiabatic orbit.
- `genera_figure_*.py`, `fig_*.py` — regenerate paper figures (need `paper_style.py` on the path).

### Compile the papers
```bash
# from NonStationaryMetrics/; each line runs in its own subshell
(cd paper1 && ./build_paper1.sh)          # Paper I: regenerates its figures, then compiles
(cd paper2 && pdflatex paper2.tex && bibtex paper2 && pdflatex paper2.tex && pdflatex paper2.tex && pdflatex paper2.tex)
```
Paper II reads its figures from `paper/Immagini` (`\graphicspath`); since figures are not
tracked, regenerate them with their scripts and run
`python3 paper2/provenance/collect_figures.py --collect` before compiling.
Each script is run from its own directory (paths from the repository root):

| Script | Figures of Paper II |
|---|---|
| `NonStationaryMetrics/ThakurtaMetric/fig_phi_validation_corrected.py` | `fig_phi_validation_true_dynamic` (also writes the raw NPZ data) |
| `NonStationaryMetrics/ThakurtaMetric/fig_phi_closed_validation.py` | `fig_phi_closed_validation` |
| `NonStationaryMetrics/ThakurtaMetric/fig_phi_penetrating.py` | `fig_phi_penetrating` |
| `NonStationaryMetrics/ThakurtaMetric/genera_indicatrici_separate_codex.py` | `fig_indicatrici_thakurta_kerr` |
| `NonStationaryMetrics/ThakurtaMetric/genera_figure_thakurta_rami.py` | `fig_thakurta_kerr_rami` |
| `NonStationaryMetrics/ThakurtaMetric/cuspide_ergosfera.py` | `fig_thakurta_cuspide_ergosfera` |
| `NonStationaryMetrics/ThakurtaMetric/thakurta_kerr_plunge_t_tau.py` | `fig_thakurta_kerr_plunge_t_tau` |
| `NonStationaryMetrics/ThakurtaMetric/thakurta_kerr_plunge_levels.py` | `fig_thakurta_kerr_plunge_levels` |
| `NonStationaryMetrics/ThakurtaMetric/inversione_conformale_AJ.py` | `fig_thakurta_kerr_inversione_AJ` |
| `KerrScripts/fig_atlas_tau_t.py` | `fig_atlas_t`, `fig_atlas_tau` |
| `KerrScripts/fig_master_penetration_taut.py` | `fig_master_penetration_taut` |
| `KerrScripts/fig_penetrate_bounce_orbit.py` | `fig_penetrate_bounce_orbit` |
| `KerrScripts/kerr_penetration_threshold.py` | `fig_penetration_threshold` |
| `KerrScripts/kerr_separatrix_trajectories.py` | `fig_separatrix_3traiettorie` |
| `KerrScripts/kerr_tbranch_separatrix_figure.py` | `fig_tbranch_separatrix` |
| `KerrScripts/kerr_jcm_capture_figure.py` | `fig_jcm_capture` |
| `KerrScripts/kerr_jpt_genus2_figure.py` | `fig_jpt_genus2` |
| `KerrScripts/bvp_estremi_fissi.py` | `fig_bvp_estremi_fissi` |
| `KerrScripts/colormap_conforme_asimmetrico.py` | `fig_colormap_conforme_asimmetrico` (about 5 min) |
| `KerrScripts/colormap_conforme_estremi_fissi.py` | `fig_colormap_conforme_estremi_fissi` (about 5 min) |
| `KerrScripts/colormap_spin_asimmetrico.py` | `fig_colormap_spin_asimmetrico` (slow: about 2 hours on one core) |
| `KerrScripts/colormap_spin_estremi_fissi.py` | `fig_colormap_spin_estremi_fissi` (slow: about 2 hours on one core) |

---

## How to reproduce a result
1. Open `paper1/verification/README.md` or `paper2/verification/README.md`, find the result and its script.
2. Run the script (`python3 …` or `sage …`); each prints its own numeric verification.
3. For independent confirmation, run the matching `paper/crosscheck_*.wl` with `wolframscript`.

The closed forms are: **symbolic coefficients** (rational in M,a,E,J and r_d,Jc — universal
formulas) **×** tabulated special functions (Weierstrass/Kleinian σ,ζ,℘; Brown–Levin Γ̃;
genus-2 elliptic polylog via θ nome series) **+** period-level data (τ, marked points, ζ,℘/θ
at those points, Ce,C0, α,β) evaluated per-curve — as with `K(m)` in the pendulum.
