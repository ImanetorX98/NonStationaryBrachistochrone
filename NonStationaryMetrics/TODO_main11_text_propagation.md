# TODO — propagate the main11 fix to the manuscript text

The **scripts** are fixed and verified (commit `1ee9794`; oracle `tests/test_adiabatic_noreg.py`,
reduction `ThakurtaMetric/adiabatic_offshell_reduction.py`). The **text** of `paper/main.tex`
and `paper/main_prd_revtex.tex` still reflects the pre-main11 state and must be updated.

## Core technical fix (Thakurta-Kerr adiabatic correction)
- [ ] **Eq. (40) / eq:adiab-exact (body + App. C):** Euler operator on the source becomes the
      FULL `D = E_eff d_E_eff + J_eff d_J_eff + P_r d_P_r`; the anchored integral is
      `S_D = int (Theta H + P_r H_Pr) dlambda`. Keep `delta P_r = (lambda Theta H - S_D)/H_Pr`,
      `delta phi = int [ G_Pr delta P_r - lambda Theta G ] dr`.
- [ ] **App. C — normalized dynamics:** state the factorization `H_br = A * Hbar(r, p_r/A, p_phi/A; Ehat/A)`
      and the true canonical flow `dP_r/ds = -Hbar_r - (A'/A) P_r` (the dilation term). Derive `S_D`.
- [ ] **App. C — new closed-form block:** add the reduction `int P_r H_Pr dlambda = int p_r dr`
      (radial action) `= const [ sum a_k U_k + third-kind at r_pm ]`, exact rational coefficients,
      SAME spectral curve `R6` as the on-shell reduction. Cite `adiabatic_offshell_reduction.py`.
      State clearly: level-B reduction done; the outer length-two wrap (genus-2 dilog) is level-C, open.
- [ ] **Fig. 10 caption (fig:phi-true-dynamic):** now vs the TRUE canonical flow; on-shell term
      slope ~1.0, exact `S_D` slope ~2.1 (t, tau, second set). Replace old 1.95/1.97/2.16.
- [ ] **App. C slope sentence:** replace `1.95 +/- 0.01 (t), 1.97 +/- 0.02 (tau), 2.16 (set2)` with
      the true-flow values (leading ~1.0, exact ~2.1); Vaidya still 2.00.
- [ ] **Tab. A3 (tab:script-map):** update slopes and SHA-256:
      `adiabatic_first_order_exact.py` = `3cd72fd798e3`,
      `fig_phi_validation_corrected.py` = `92909d1d6a50`,
      `vaidya_first_order_offshell.py` = `5dc9fe177f2f` (unchanged);
      add `adiabatic_offshell_reduction.py` and `tests/test_adiabatic_noreg.py` rows.

## Presentation qualifications (main11 sec. 5)
- [ ] **Eq. (32):** rename `phi = phi0 + eps(C+psi)` -> `phi_on` (C+psi is on-shell only; Fig. 10
      shows it has O(eps) physical error). Or add the costate term explicitly.
- [ ] **Eq. (38):** call it "complete **on-shell** coefficient" in the sentence itself.
- [ ] **Eq. (51):** the unified theorem must be the COMPLETE (on- + off-shell) form; the on-shell
      formula is a corollary/component. `+O(lambda_dot)` for the on-shell-only trajectory.
- [ ] **App. B title:** "fully explicit first-order adiabatic brachistochrone" ->
      "explicit **on-shell** separatrix adiabatic forms".
- [ ] **Abstract / Tab. 1:** keep the three-claim split; ensure "on-shell closed form" vs
      "complete first-order (with off-shell S_D)" vs "explicit off-shell assembly (open)" is exact.
- [ ] **Slow variable (sec 5.2):** disambiguate eta (conformal time) vs accumulated clock t(r)/tau(r);
      write `dA/dlambda = (dA/deta)(deta/dlambda)` per branch.
- [ ] **Protocol:** Eq. (40) is fixed-charge/launch, not fixed-endpoint; fixed spatial endpoints
      need `J = J0 + eps J1` from the angular condition (state, or compute J1).

## Other analytic/physical (main11 sec. 6)
- [ ] **PMP "well-posed":** soften to "normal PMP extremals are well defined with a unique maximizing
      direction in the timelike-selector, nondegenerate domain" (not global existence/minimality).
- [ ] **J costate vs mechanical momentum:** present the full branch-by-branch canonical/Legendre map
      including P_r=p_r/A and the time-dependent generator (not just the axial costate).
- [ ] **Outgoing Vaidya:** replace "evaporation" for m'<0 in the ingoing chart with "negative ingoing
      mass slope" / "formal ingoing continuation" (Fig. 6 and phenomenology).
- [ ] **Ergosphere / horizon nomenclature (see memory [[thakurta-kerr-not-black-hole]]):** TK central
      object is a COMPACT OBJECT, not a BH; `Delta=0` is the "seed Kerr null surface" / "conformal
      Killing-horizon candidate", not a horizon. Inner curves = "analytically continued exterior
      extremal". Keep BH/horizon language for Vaidya only.
- [ ] **Fig. 14:** separate the two residuals (ODE-closed ~1e-5 near instability vs 2.8e-10 global)
      and their domains.

## Bibliography (main11 sec. 8) — titles now needed (switch bib style if iopart-num hides them)
- [ ] [15] Caponio-Javaloyes-Sanchez: add title, Mem. AMS 300 (2024) no. 1501, 121 pp., DOI 10.1090/memo/1501.
- [ ] [31] Brown-Levin: add title "Multiple Elliptic Polylogarithms" + arXiv:1110.6917.
- [ ] [43] Baune et al.: add title + DOI 10.1088/1751-8121/ad8197 (authors already fixed).
- [ ] [44] D'Hoker-Schlotterer: add title "Fay identities for polylogarithms on higher-genus Riemann
      surfaces" + arXiv:2407.11476 with version/date.
- [ ] [45] D'Hoker-Hidding-Schlotterer: add title + DOI 10.4310/CNTP.250531031558 (CNTP 19(2) 355-413, 2025).
- [ ] [46] D'Hoker-Enriquez-Schlotterer-Zerbini: add title + DOI 10.1007/s00220-025-05540-x.
- [ ] [49] DLMF: current release at access date is 1.2.7 (2026-06-15), not 1.2.4; state which was used.
- [ ] [54] "Nario" -> **Natario** (Jose Natario / Jose Natario), add title, pp. 2579-2586, DOI 10.1007/s10714-009-0781-2.
- [ ] [58] SageMath 10.9 released 2026-05-04/05 (NOT 2026-02-22); fix date; cite tag/commit or SW Heritage.
- [ ] [59] abelfunctions fork: give exact commit + persistent archive (not "et al. 2026").
- [ ] [64] repo: remove `<SHA-TBD>`; fill Zenodo DOI once minted. Referee saw NSB commit
      `52522c77ebdac6d367178b311330ff8cdee0d9c9`.

## Reproducibility (main11 sec. 9)
- [ ] Tagged/archived release tied exactly to the revised PDF; complete commit in Data Availability.
- [ ] Environment file with exact versions (Python, SymPy, NumPy, SciPy, SageMath, Mathematica, abelfunctions fork).
- [ ] Machine-readable symbolic coefficients (not only decimal examples).
- [ ] Per slope plot: tolerances, precision, residual norm, sub-arc, refinement.
- [ ] Automated test of the ORIGINAL H_t/H_tau flow vs the normalized representation with -alpha P_r.
      (DONE: `tests/test_adiabatic_noreg.py` — reference it in the paper.)
