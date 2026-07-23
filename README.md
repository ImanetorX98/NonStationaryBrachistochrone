# Constrained-worldline brachistochrones in non-stationary spacetimes

Source code, closed-form derivations and manuscript for the fastest constrained
"rail" worldline between two events, extended from stationary to **non-stationary**
spacetimes — FLRW, Vaidya (dynamical black hole), and Thakurta–Kerr (conformal
rotating) — built on the stationary-Kerr optical/genus-2 machinery.

Every closed form has a script that **derives it symbolically and verifies it
numerically**; the substantive results are additionally **cross-checked
independently in Mathematica** (`paper/crosscheck_*.wl`).

> Generated binaries (`*.png`, `*.pdf`, videos) are **not** tracked; regenerate from
> the figure scripts.

> **Note on structure.** This folder gathers the non-stationary study. The broader
> **stationary-Kerr base machinery** (`KerrScripts/`, `KerrMetric/`) remains at the
> repository root as **siblings** (`../KerrScripts/`, `../KerrMetric/`); the paper
> builds on it. The genus-2 closure scripts used in this study are copied here under
> `KerrSessionScripts/`.

---

## Layout

| Path | Contents |
|------|----------|
| `paper/` | LaTeX sources (`main.tex` CQG/iopart, `main_prd_revtex.tex` PRD) **and** all Mathematica cross-checks `crosscheck_*.wl`. |
| `FLRWmetric/` | FLRW (degenerate base) analysis + the perturbed-FLRW genus-jump proof of concept. |
| `VaidyaMetric/` | Vaidya (dynamical BH, parameter = mass `m`) analysis, closed forms, figures. |
| `ThakurtaMetric/` | Thakurta–Kerr (conformal rotating, parameter = `A` via `E_eff=Ê/A`) analysis, closed forms, figures. |
| `KerrSessionScripts/` | Genus-2 machinery + the **generic-J (genus-2) closure pipeline** used in this study. |
| `reproduce_reductions.py` | One-file check of the algebraic heart. |
| `paper_style.py` | Shared matplotlib style helper. |

### Reference documents (read these first)
| File | What it collects |
|------|------------------|
| `VERIFICATION_STATUS.md` | **Master table** of every result and its verification (SYM / NUM / SAGE / Mathematica). Start here. |
| `SEP_COEFF_SYMBOLIC.md` | Separatrix (genus-1) closed forms: symbolic `b1,b2,b3`, clock residues, `Ce,C0`, period-level hierarchy. |
| `GENUS2_CLOSED_FORM.md` | Generic-J (genus-2) closed form: `δφ = ½Ê Σ Q_kj W_kj + G_alg`, symbolic coefficients, tabulated special functions. |
| `progress.md` | Chronological log of the whole derivation, brick by brick. |
| `EXPLICIT_FORMS_PROGRESS.md`, `RESULTS_ROADMAP.md`, `ALGORITHM_adiabatic_closed_forms.md`, `UNIFIED_ADIABATIC.md`, `sumUp.md`, `paperOutline.md` | Earlier progress notes, the reusable algorithm, figure index, outline. |

### Fixed-endpoint no-inversion (Lemma B) + referee response (recent)

> **Status (tiered, not absolute).** The fixed-endpoint no-inversion is **proved**
> for `r0 <= R*(E)` and, in the static case, asymptotically as `r0 -> infinity`;
> **machine-certified** (interval arithmetic) at representative `r0`; and
> **conjectural** in general. The static/frozen branch ordering does not invert
> under the outgoing-clock replacement; the full outgoing-Vaidya optimal-control
> boundary-value problem **remains open**.

| File | What it does |
|------|--------------|
| `RESULTS_lemmaB_frozen_schwarzschild.md` | Master summary of the frozen-Schwarzschild no-inversion: Lemma A closed; Lemma B via elementary (`r0<=R*`) + closed-form large-`r0` asymptotic + CAP. |
| `no_inversion_reduction.py` | Reduction + Lemma A closed-form proof (`symbolic_full_proof`). |
| `no_inversion_schwarzschild_closedform.py` | Closed form for `Phi_tau'` (IBP); elementary grazing/quarter regimes; single-crossing criterion. |
| `no_inversion_schwarzschild_asymptotic.py` | Large-`r0` asymptotic: `Phi=arctan(...)-b^{-3/2}sqrt(...)`, single non-degenerate max, `Phi''(x_pk)=-1/2 b^{9/4} V0^{-5/4}`. |
| `no_inversion_schwarzschild_CAP.py`, `..._CAP_r0_10.py`, `..._CAP_grid.py` | Computer-assisted proof (mpmath interval arithmetic): single-crossing certificates. Complete at `r0=10`; grid `E in {1.2,1.6,2.5} x r0 in {8,12}` (`CAP_grid_certificates.txt`). |
| `verify_lemmaB_mathematica.wls` | Independent Mathematica cross-check of all no-inversion derivations. |
| `REFEREE_RESPONSE_main8.md` | Point-by-point response to the main8 referee report (DONE vs PLANNED). |
| `ThakurtaMetric/adiabatic_offshell_validation.py`, `adiabatic_tk_geodesic_check.py` | Validation of the adiabatic correction against the true optical-metric geodesic: `1/2`-Euler captures ~98%; off-shell (`H2!=0`, costate `p_eta`) term is the ~2% residual (referee 4.6). |

### Complete first-order extended-Hamiltonian correction (Eq. 40, on-shell + off-shell)
| File | What it does |
|------|--------------|
| `ThakurtaMetric/adiabatic_first_order_exact.py` | Exact canonical PT of the non-autonomous optical Hamiltonian: `delta p_r=(lambda ThetaH - S)/H_pr`, `S=int Theta H dlambda`; closes the **Thakurta–Kerr** first order to `O(eps^2)` (residual slope ~1.86). |
| `VaidyaMetric/vaidya_first_order_offshell.py` | **Vaidya** `v`-branch analogue (`Theta=m d_m`, sign flip for accreting `m`, terminal anchoring of `p_v`); closes to `O(eps^2)` (slope **2.00** full vs **1.00** on-shell only). Emits `fig_vaidya_offshell`. |
| `ThakurtaMetric/fig_phi_validation_corrected.py` | True-dynamics validation figure (`fig_phi_validation_true_dynamic`): leading on-shell (slope ~1, ~2% physical error) vs exact Eq. (40) (slope ~1.86) vs the true non-autonomous flow. |

---

## Install / run

```bash
python3 -m pip install -r requirements.txt      # numpy scipy sympy matplotlib mpmath Pillow
```

- **Python (sympy/mpmath)** — algebraic reductions + genus-1 (Weierstrass via
  `mpmath.jtheta`) closed forms. Most `.py` scripts. Standalone: `python3 <script>.py`.
- **SageMath + `abelfunctions`** — genus-2 objects (Riemann matrix, `RiemannTheta`,
  Abel map): the `.sage` scripts. Run inside `sage` from a neutral cwd.
- **Mathematica / `wolframscript`** — independent cross-checks `paper/crosscheck_*.wl`
  (native `WeierstrassP`, `GroebnerBasis`, `Series`, `Resultant`).

Verification scripts print numerical residuals; ~machine precision means the identity holds.

> **Run all commands from inside this folder** (`.../NonStationaryMetrics/`), so paths like
> `VaidyaMetric/…`, `KerrSessionScripts/…`, `paper/…` resolve. The stationary-Kerr base lives
> one level up, so reference it as `../KerrScripts/…` and `../KerrMetric/…`.

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

### Compile the paper
```bash
cd paper && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

---

## How to reproduce a result
1. Open `VERIFICATION_STATUS.md`, find the result and its script.
2. Run the script (`python3 …` or `sage …`); each prints its own numeric verification.
3. For independent confirmation, run the matching `paper/crosscheck_*.wl` with `wolframscript`.

The closed forms are: **symbolic coefficients** (rational in M,a,E,J and r_d,Jc — universal
formulas) **×** tabulated special functions (Weierstrass/Kleinian σ,ζ,℘; Brown–Levin Γ̃;
genus-2 elliptic polylog via θ nome series) **+** period-level data (τ, marked points, ζ,℘/θ
at those points, Ce,C0, α,β) evaluated per-curve — as with `K(m)` in the pendulum.
