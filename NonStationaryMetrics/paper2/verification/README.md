# Verification of Paper II

Scripts that check the explicit computations of `paper2/paper2.tex`: symbolic identities,
polynomial certificates, interval enclosures and high-precision evaluations. Most scripts
exit non-zero when a check fails; read the printed PASS/FAIL lines in any case, since an exit
code alone is not the evidence. The adiabatic pipeline and the computer-assisted
certificates are indexed separately, in the result-to-script map (`tab:script-map`) and in
`paper2/provenance/MANIFEST.tsv` (`python3 paper2/provenance/make_provenance.py --check`; on a
fresh clone run `python3 ThakurtaMetric/fig_phi_validation_corrected.py` once first, from
`NonStationaryMetrics/`, since its raw
`.npz` output is not distributed; the comparison is exact in the reference environment of
`requirements-lock.txt`).

```
python3 <script>.py                 # from this folder
wolframscript -file <script>.wls    # from this folder
```

**Coverage column:**
- **full**: the checks test the claim as the manuscript states it;
- **partial**: the script tests identities or samples, and the rest of the argument is in
  the text;
- **sampled**: a scan at chosen parameters, not a proof.

## Randers reduction, separation, conjugate points

| Script | Claims | Manuscript | Coverage |
|---|---|---|---|
| `verify_randers_convexity.py` | arrival cost vs Jacobi–Maupertuis; convexity and its boundary; gauge dependence of `b`, invariance of `db` | `prop:randers-reduction`, `eq:FT-vs-JM`, `eq:db` | full |
| `verify_randers_finsler.wls` | Randers/Zermelo data of Thakurta–Kerr; null Finsler boundary as the control-domain edge | `sec:randers-reduction`, `lem:compact` | full for the null data; the finite-energy cost is in `verify_randers_convexity.py` |
| `verify_geodesic_separation.py` | static certificate `K_t < 0` for `Ehat^2 >= 3/2`, `r >= 2M` | `eq:optical-gauss`, `cor:static-separation` | full |
| `verify_rotating_separation.py` | the closed form of `B` and the geodesic-curvature identity (symbolic); the magnetic Jacobi equation against the full linearised flow, conservation of `p_phi`, first and second variation of the arrival cost (numerical, at sampled parameters) | `prop:rot-separation`, `eq:magnetic-jacobi`, `eq:magnetic-B` | partial: identities symbolic, flow and variation checks sampled |
| `verify_rail_surfaces.wls` | extrinsic geometry of rail-swept surfaces in `A(eta)^2 g_Kerr` with arbitrary `A`: conformal change on the same immersion, seed and unit-normal umbilicity minor, radial family and the `A r = c` (areal-radius) criterion, isolated instants only for `A_eta/A` constant, common-launch leading terms; elimination of the tangential component of the magnetic variation field | `lem:rail-surfaces`, `eq:rail-surface-minor`, `app:rail-surfaces`, introduction | full for the identities; the `o(tau)` estimates of (iii) are in the text. A run is archived in `verify_rail_surfaces_run.txt` |
| `verify_rail_surfaces_sympy.py` | independent SymPy route for the same identities (equatorial Christoffel symbols rebuilt from the definitions; polar component in 4D) | `lem:rail-surfaces` | full for (i), (ii). A run is archived in `verify_rail_surfaces_sympy_run.txt` |
| `verify_optical_planes.py` | central planes of the frozen optical data `n α_K + β_K` off the equator: derivation from the two conditions, equatorial reduction, parity in `a`, equator totally geodesic, inclined planes totally geodesic for `n α_K` up to `O(a²)` but not for the arrival-time data (normal force `eq:inclined-plane-force`), proper-time data even in `a`; ordinary umbilicity for the Riemannian metrics `h_t`, `h_τ` (conformal to one metric `h̄`): equatorial and meridional planes totally geodesic exactly, inclined planes with no umbilical point at order `a²` (`eq:inclined-plane-umbilicity`), with a mutation test and an exact-metric Cartesian check; independent numerical integration of the extremals as magnetic geodesics | appendix paragraph *Central planes of the optical data*, introduction | full for the symbolic statements at frozen `A`; sampled for the numerical check. A run is archived in `verify_optical_planes_run.txt` |
| `verify_prograde_asymmetry.py` | sign of `B'` on the exterior; `W - K > 0` where `T^phi > 0` | `prop:prograde-asym` | full |
| `verify_wsup_window.py` | direction-independent bound `W <= W_sup(r)`; the two interval-certified windows | `eq:Wsup` | full on the stated windows |
| `verify_myers_optical.wls` | optical Gauss and scalar curvature; sign certificates (Reduce) excluding the Myers and curvature-operator positivity hypotheses; the trace identity `R_opt = 2(2K_t + K_tan)` | `eq:optical-scalar`, the curvature-operator paragraph of `sec:submersion` | full for these identities and certificates; the spectrum and the out-of-plane formula are derived in the text; one diagnostic scan. A full run is archived in `verify_myers_optical_run.txt` |
| `verify_3d_sectional.wls` | tangential sectional curvature `K_tan`; its sign at the circular orbit and in the far field | `eq:ktan` | partial: formulas exact, energy samples |
| `verify_outofplane_jacobi.wls` | rotational Killing field as out-of-plane Jacobi field; first zero at azimuth `pi` | `eq:outofplane-J` | full for out-of-plane conjugate points; not a 3D minimality statement |
| `verify_conjugate_maxwell.wls` | sign of `K` on the exterior at `Ehat^2 = 2`; energy threshold | `prop:noconj`, Maxwell discussion | partial: the Maxwell consequences are argued in the text, not computed |
| `verify_optical_completeness.wls` | completeness of the arrival-time optical surface; residue `2M` of `ds/dr` at `r = 2M` | proof of `prop:globalmin` | partial: the limits are exact, the global conclusion is geometric |
| `verify_sufficient_conditions.wls` | radial eikonal saturated; inner negative-curvature region below threshold | `prop:hjb-frozen` | partial: radial identity only; energies sampled |
| `verify_two_branch_optics.wls`, `verify_which_clock.wls` | which optical surface belongs to which clock; curvatures and limits per branch | `sec:randers-reduction` | full for formulas; positivity scans sampled |

## Submersion, selectors, drift

| Script | Claims | Manuscript | Coverage |
|---|---|---|---|
| `verify_tension_sympy.py` | fibre acceleration; Lorentzian tension from the connections; not a harmonic morphism; sign change for `Ehat^2 < 3/2` | `eq:fibre-accel`, `eq:tension-residual`, `sec:submersion` | full |
| `verify_submersion_link.wls` | horizontal optical factor; Vaidya rail-drift identity | `eq:dilation`, `sec:submersion` | full for those; the fibre geometry is in `verify_tension_sympy.py` |
| `verify_dilation_fibre.wls` | the Perlick dilation varies along the fibres exactly when `A` runs | `eq:dilation-fibre` | full; no Clairaut conclusion drawn |
| `verify_tk_ckv_rail.wls` | conformal-Killing selector of Thakurta–Kerr; drift rate `A'/A` | `eq:tk-drift` | full; non-Einstein check at `a = 0` |
| `verify_soliton_rail.wls` | drift identity for a Ricci-soliton selector, `lambda + Ric(u,u)` | `eq:soliton-drift` | full for the identity |
| `verify_weighted_soliton.wls` | generalised m-Bakry–Émery tensor at `m = -2`, `phi = 2 ln A` (not the metric weighting of the cited soliton papers): `Ric_phi^m = mu g`; `(1/2) L_W g + Ric_phi^m = lambda g`; `lambda = c` for `A = A0 exp(c eta)`; geodesic drift; equatorial `gbar^{eta eta}` depends on `r` for `M > 0`; full Ricci tensor of Thakurta–Kerr and of its `a = 0` case, with a control that must fail | `eq:weighted-soliton` | full for the tensor identities; non-constancy of `lambda` sampled at one non-exponential `A`; the necessity of `A = A0 exp(c eta)` is the analytic argument in the text, with its ingredient (`r`-dependence of `gbar^{eta eta}`, sampled point) checked. A run is archived in `verify_weighted_soliton_run.txt` |
| `verify_weighted_seed.wls` | weighted seed `A^{-2} g = g_Kerr`; divergence and gradient of the weight | "The weighted seed" | full for the identities; the cited framework's other hypotheses are not checked |
| `verify_thrust_bound.wls` | minimum thrust; finite at the stationary limit, divergent at freezing | `eq:thrust-bound` | full |
| `verify_gp2002_fork.wls` | where the controlled rail departs from Giannoni–Piccione (2002) | introduction | full for the drift identity (about 4 min; a run is archived in `verify_gp2002_fork_run.txt`) |
| `verify_conformal_killing_horizon.wls` | twist of the seed; equatorial criterion at the stationary limit | `sec:domain` | full; no global horizon claimed |
| `verify_vaidya_tube.wls` | Kodama selector, divergence, Misner–Sharp mass, signature of the apparent-horizon tube | Vaidya discussion | full |
| `verify_clock_comparison.wls` | scalar comparison lemma for the clock (Paper I) | Paper I Lemma I.B | partial: sampled comparisons |

## Separatrices, classification, adiabatic response

| Script | Claims | Manuscript | Coverage |
|---|---|---|---|
| `verify_cusp_corner.wls` | local shape at `r_e`: semicubical cusp for `0 < |J| < J_c` (radial at `J = 0`); cusp amplitude divergent as `|J| -> J_c`; corner slope at the marginal values; for `|J| > J_c` the extremal reflects outside `r_e` | `prop:classification` | full for the shape |
| `verify_marginal_hamilton.wls` | marginal dynamics: `J = -J_c` reaches `r_e` only asymptotically, `J = +J_c` crosses it at a finite rate | `prop:classification` | full for the rates |
| `verify_exterior_separatrix.wls` | exterior retrograde separatrix `r_d = 3.5139M`, `J_c^- = -8.0535` | `prop:exterior-sep` | full at the stated parameters |
| `verify_same_launch_ordering.wls` | same-launch depth ordering at `J_t = J_tau`: turning conditions, identity `J_t - J_tau = -(A/(E Q)) cond_t`, `cond_t < 0` on the whole prograde domain (exact `Reduce`) | `prop:same-launch` | full on the stated domain (CAD) |
| `verify_tau_separatrix_dictionary.py` | marked-point dictionary of the rotating tau-separatrix | `eq:sep-phi` | full |
| `verify_tau_branch_drift_clock.py` | tau-branch drift clock against the Hamiltonian flow; weight-two assembly with its horizon term | `eq:clock-eta-tau`, `eq:psi-split` | full at the stated parameters |
| `verify_adiabatic_raw_codex.py` | refits the adiabatic slopes from the raw NPZ data (not tracked: on a fresh clone run `ThakurtaMetric/fig_phi_validation_corrected.py` first) | `tab:adiab-valid` | partial: not a check of the whole pipeline, nor an automatic comparison with the printed table |
| `verify_independent_rederivation.py`, `verify_independent_rederivation_codex.wls` | re-derived from the printed equations: the quartic `Q_2` and the marginal factorisation; the exterior separatrix digits and spin scan; the optical-curvature coefficients; the conformal transfer map; the marginal rates and action identities (Props. 3.1–3.2, section 2.2, Paper I Vaidya) | several | full for the listed identities; section D of the `.wls` is the Riemannian assembly, not the tension |

## Fixed-endpoint no-inversion

| Script | Claims | Manuscript | Coverage |
|---|---|---|---|
| `verify_lemma_B.py` | monotonicity regions (i)–(iv), single peak for `E^2 >= 3/2`, `W' < 0` for `r >= 5M/2`; polynomial certificates; exact `C_c`; the two interval-certified counterexamples | `prop:lemmaB`, `app:lemmaB` | full for identities and certificates; the analytic steps are in the appendix |
| `../../no_inversion_reduction.py` | Lemma (A) factorisation, symbolic in `M` | proof of (A) | full; its Lemma B blocks are sampled diagnostics |
| `../../no_inversion_schwarzschild_CAP_grid.py` | interval certificates at seven static configurations, windows to the half-potential radius | CAP paragraph | full at the stated binary parameters (about 90 min) |

## Observational context and probe experiments

The scripts of `sec:observational-context` and `app:tk-energy` are in `observational/`, with their own
README: physical scales of the catalogued objects, Komar and Hawking–Hayward quantities, energy
conditions, the rotating benchmarks, and the `probe_*` scripts of the proposed Solar-System two-probe
comparison, of the solar-sail comparison and of the clock thought experiment around Gaia BH3.

| Script | Claims | Manuscript | Coverage |
|---|---|---|---|
| `verify_rule_of_two.py` | relativistic rule of two: in a static spacetime the maintained-charge arrival rail has proper acceleration normal to the velocity with magnitude `2 g_perp` (symbolic identity; exact four-acceleration from the Christoffel symbols along Schwarzschild rails at `r_min = 10M, 50M, 1e5 M`, `Ê = 1, 1.05, 1.3`, to 40 digits; mutation of the factor rejected) | `prop:rule-of-two` | full for the identity; sampled for the kinematic check. A run is archived in `verify_rule_of_two_run.txt` |
| `observational/probe_control_law.py` | weak-field rail thrust `2 mu sin(alpha)/r^2`; chord thrust; convexity; unique turning radius | Solar-System paragraph of `sec:observational-context` | full (symbolic) |
| `observational/probe_saturn_comparison.py` | Saturn flyby, Hill exit, rails to 4 and 3 AU with prescribed planets; daughter lags `∝ d^2` | Solar-System paragraph | sampled: one model phase, circular planets; a concept, not a mission certification |
| `observational/probe_conjugate_points.py`, `observational/verify_conjugate_points.wls` | curvature `K = Lap(Omega) - 2 grad(Omega)^2/f` of the rail optical metric; conjugate points from L4 (Earth–Moon) and on a hyperbolic rail around Jupiter, with the sign change of the daughter lag; no Jacobi-field zero on the Saturn rails within 9 years | Solar-System paragraph | sampled: two L4 rays, one Jupiter rail, one model phase |
| `observational/probe_solar_sail.py` | solar-sail clock comparison: arrival- and onboard-cost rails at Ê = 1 (separation 0.9215 GM/c² to first post-Newtonian order, inversions), Newtonian sail law (rule of 2, modulation q = 2/(β_max sin α), lightness needed by mother and daughters, absorption), daughters with X = t − τ linear in d, energy sensitivity of X | Solar-sail paragraph of `sec:observational-context` | full at the stated parameters; ideal actuator, a concept, not a mission design |
| `observational/probe_gaia_bh3.py` | Gaia BH3 clock thought experiment (M = 32.70 M☉, companion from the Gaia Collaboration 2024, a* = 0.5 illustrative): encounter clock effect `6πMa cos i/√(r(r−3M))` (symbolic first-order derivation and exact equatorial closed form); single-station circuit (rail plus arc) with coordinate signal `2∮β_K` and the station's proper-time reading, at `r_min = 1e5, 1.5e5` (benchmark) `and 2e5 M`; first-order shift of the station arrival difference by the companion's quadrupole tide for free sliding (law B) and for an ideal compensation of the modelled tide with a mass-only model error (law A), with the static cancellation, a convergence check and a nonlinear integration of the sliding motion; Gaia BH1 and BH2 at the same `r_min/M` for comparison; thrust budget (rail peak `2 g_s`, arc, turns at maintained charge) and ideal photon-rocket ratio; static station (proper acceleration, budget per circuit and per year, ideal transitions from and to the free orbit, advance of a free station); transfer for scale | `par:obs-kerr-thought`, `eq:encounter-clock` | full for the formulas; sampled for the companion (48 orientations, two epochs) |
| `observational/probe_kerr_m87.py` | earlier Kerr thought experiment scaled to M87* (rails at 120M and 38M), replaced in the manuscript by Gaia BH1 on 24 October 2026; kept for the archived v1.8.x releases | none in the current text | full at its stated parameters |
| `observational/ephemeris_conformal_bound.py` | Newtonian parameter `m A` of `A^2 g_Schw`; `H0` against the EPM2019 bound on `d(GM_sun)/dt` | drift paragraph of `sec:observational-context` | full for the limit; the bound is the cited value |
