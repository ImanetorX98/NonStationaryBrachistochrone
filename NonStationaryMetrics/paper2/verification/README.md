# Verification of Paper II

Scripts that check the explicit computations of `paper2/paper2.tex`: symbolic identities,
polynomial certificates, interval enclosures and high-precision evaluations. Most scripts
exit non-zero when a check fails; read the printed PASS/FAIL lines in any case, since an exit
code alone is not the evidence. The adiabatic pipeline and the computer-assisted
certificates are indexed separately, in the result-to-script map (`tab:script-map`) and in
`paper2/provenance/MANIFEST.tsv` (`python3 paper2/provenance/make_provenance.py --check`).

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
| `verify_gp2002_fork.wls` | where the controlled rail departs from Giannoni–Piccione (2002) | introduction | full for the drift identity |
| `verify_conformal_killing_horizon.wls` | twist of the seed; equatorial criterion at the stationary limit | `sec:domain` | full; no global horizon claimed |
| `verify_vaidya_tube.wls` | Kodama selector, divergence, Misner–Sharp mass, signature of the apparent-horizon tube | Vaidya discussion | full |
| `verify_clock_comparison.wls` | scalar comparison lemma for the clock (Paper I) | Paper I Lemma I.B | partial: sampled comparisons |

## Separatrices, classification, adiabatic response

| Script | Claims | Manuscript | Coverage |
|---|---|---|---|
| `verify_cusp_corner.wls` | local shape at `r_e`: semicubical cusp for `0 < |J| < J_c` (radial at `J = 0`); cusp amplitude divergent as `|J| -> J_c`; corner slope at the marginal values; for `|J| > J_c` the extremal reflects outside `r_e` | `prop:classification` | full for the shape |
| `verify_marginal_hamilton.wls` | marginal dynamics: `J = -J_c` reaches `r_e` only asymptotically, `J = +J_c` crosses it at a finite rate | `prop:classification` | full for the rates |
| `verify_exterior_separatrix.wls` | exterior retrograde separatrix `r_d = 3.5139M`, `J_c^- = -8.0535` | `prop:exterior-sep` | full at the stated parameters |
| `verify_tau_separatrix_dictionary.py` | marked-point dictionary of the rotating tau-separatrix | `eq:sep-phi` | full |
| `verify_tau_branch_drift_clock.py` | tau-branch drift clock against the Hamiltonian flow; weight-two assembly with its horizon term | `eq:clock-eta-tau`, `eq:psi-split` | full at the stated parameters |
| `verify_adiabatic_raw_codex.py` | refits the adiabatic slopes from the archived NPZ data | `tab:adiab-valid` | partial: not a check of the whole pipeline, nor an automatic comparison with the printed table |
| `verify_independent_rederivation.py`, `verify_independent_rederivation_codex.wls` | re-derived from the printed equations: the quartic `Q_2` and the marginal factorisation; the exterior separatrix digits and spin scan; the optical-curvature coefficients; the conformal transfer map; the marginal rates and action identities (Props. 3.1–3.2, section 2.2, Paper I Vaidya) | several | full for the listed identities; section D of the `.wls` is the Riemannian assembly, not the tension |

## Fixed-endpoint no-inversion

| Script | Claims | Manuscript | Coverage |
|---|---|---|---|
| `verify_lemma_B.py` | monotonicity regions (i)–(iv), single peak for `E^2 >= 3/2`, `W' < 0` for `r >= 5M/2`; polynomial certificates; exact `C_c`; the two interval-certified counterexamples | `prop:lemmaB`, `app:lemmaB` | full for identities and certificates; the analytic steps are in the appendix |
| `../../no_inversion_reduction.py` | Lemma (A) factorisation, symbolic in `M` | proof of (A) | full; its Lemma B blocks are sampled diagnostics |
| `../../no_inversion_schwarzschild_CAP_grid.py` | interval certificates at seven static configurations, windows to the half-potential radius | CAP paragraph | full at the stated binary parameters (about 90 min) |
