# Observational context of Paper II (section `sec:observational-context`, appendix `app:tk-energy`)

| Script | Claims | Output |
|---|---|---|
| `compute_observational_scales.py` | physical scales r_g, t_g, ε₀ = H₀t_g of the catalogued objects; static radial benchmark Δt, Δτ; conditional birth-mass exercise | `observational_scales.json` |
| `audit_rotating_geometry.py` | Komar J = A²Ma on BL spheres; Hawking–Hayward energy M_HH = A[m_K + (A′/A)²C] and the table of m_K, C | `rotating_geometry.json` |
| `audit_rotating_asymptotics.py` | large-r expansions of m_K and C; spin deformation of flux and shear | `rotating_asymptotics.json` |
| `audit_rotating_energy.py` | Einstein tensor of A²g_Kerr in the ZAMO frame; Hawking–Ellis type-IV and NEC/WEC boundaries | `rotating_energy.json` |
| `audit_rotating_bvp.py` | rotating common-endpoint benchmark (t and τ branches), converted with the Gaia BH3 mass | `rotating_bvp.json` |
| `audit_final_E1_benchmark.py` | rotating common-endpoint benchmark at $\hat E=1$ (t and τ branches, ZAMO speeds, Gaia BH3 units); exact K_t at $\hat E=1$ | `audit_final_E1_benchmark.json` |
| `verify_tk_komar.wls` | independent Mathematica check of J = A²Ma | printed PASS lines |
| `probe_control_law.py` | weak-field control law: rail thrust 2μ sin α/r² (twice the normal gravity), chord thrust μ sin α/r², strict convexity of the monotone-r arrival functional, unique Snell turning radius | `probe_control_law.json` |
| `probe_saturn_comparison.py` | proposed two-probe comparison: Earth–Saturn transfer (C3), flyby at −85° (pericentre, ring-plane crossings), exit from Saturn's Hill sphere, rail to 4 and 3 AU with prescribed circular planets; daughter paths sharing endpoints and initial tangent, lag ∝ d² | `probe_saturn_comparison.json` |
| `probe_conjugate_points.py` | conjugate points of maintained-energy rails: curvature of the optical metric \|dx\|²/(2Ω − C) (symbolic, and K = Eμ/(r²(Er + μ)) for Kepler); Earth–Moon rays from L4 and a hyperbolic rail around Jupiter, with daughter lags changing sign across the conjugate point (least favourable fixed-endpoint deformation, initial direction free); negative heliocentric energy, and no sign change of the transverse endpoint fields on the Saturn rails within 9 years (centred differences, convergence checked; a sampled numerical check, not a certificate) | `probe_conjugate_points.json` |
| `verify_conjugate_points.wls` | independent Mathematica check of the curvature formulas and of the L4 and Jupiter conjugate times (agreement with Python within 1e-6 d); exits with status 1 on failure | `verify_conjugate_points_run.txt` |
| `probe_solar_sail.py` | solar-sail clock comparison: optical metrics of the arrival and onboard costs at maintained Ê (symbolic), their separation at Ê = 1 on arcs symmetric about r_min = 0.1 and 0.05 au, the two inversions, thrust Δv and sail lightness β = 2/sin α, critical areal density and sail temperature; daughters r = r_t(φ) + d sin²(πu): t(d) quadratic, vertex of τ(d), X = t − τ linear in d; energy sensitivity of X; mpmath at 45 digits, agreeing with an independent Mathematica calculation | `probe_solar_sail.json` |
| `probe_kerr_m87.py` | Kerr thought experiment scaled to M87*: Ê = 1 arrival rails at 120M and 38M, a/M = −0.8, 0, 0.8; turning radius, peak proper acceleration in g, scalar budget and photon-rocket ratio, spin shift of the arrival | `probe_kerr_m87.json` |
| `ephemeris_conformal_bound.py` | slow-motion limit of A²g_Schw: Newtonian parameter m A, growing at H; comparison of H₀ with the EPM2019 value of d(GM☉)/dt/GM☉ (Pitjeva et al. 2021, eq. (1)) | `ephemeris_conformal_bound.json` |

`python3 <script>.py` from this folder (Python 3, mpmath, SymPy, NumPy, SciPy); `audit_rotating_bvp.py` reuses
`KerrScripts/colormap_spin_estremi_fissi.py` of the repository. Run `compute_observational_scales.py` first.

The `probe_*` and `ephemeris_*` scripts are self-contained (about 30 s in total); `probe_conjugate_points.py` imports `probe_saturn_comparison.py` from this folder. The Solar-System model uses planets on
prescribed planar circular orbits with J2000 mean longitudes from JPL *Approximate Positions of the Planets*
(Table 1) and the Saturn pole of the JPL kernel PCK00011; its dates are model phases, not launch windows, and
the concept is not a launch, ring-safety, power or navigation certification.
