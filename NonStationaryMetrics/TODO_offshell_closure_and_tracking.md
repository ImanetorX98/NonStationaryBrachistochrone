# TODO — off-shell closed form: remaining items + TK separatrix tracking counterterm

Status snapshot after the generic t-branch off-shell wrap closure (commits 250258b, e903544,
f15c2c4). The generic (J != Jc) off-shell wrap for the **t-branch** is now assembled in explicit
special functions and physics-anchored:
  d phi_extra = A (2nd-kind W_jk -> Kleinian zeta,sigma + one s=0 genus-2 dilog)
              + B (Delta=0 genus-2 dilogs D_{j,root} at seed-Kerr null r_pm)
              + C (elementary, from Hermite kernel-pole reduction),
verified vs the direct double integral to 1e-14 (two configs) AND against the TRUE non-autonomous
flow (closed == flow-validated dilation sub-piece to 6e-9; total slope ~2).
Scripts: KerrSessionScripts/offshell_tbranch_{wrap_assembly,thirdkind_Delta,block3_assembly,
kernel_hermite,FULL_assembly}.py, physics_anchor_offshell_closed.py.

## DONE this session
- [x] Second-kind block closed in canonical basis (symbolic in E).
- [x] Delta=0 third-kind block: Pi = sum_root rho_root Pi_root (rho_E=0), genus-2 dilogs D_{j,root}.
- [x] Kernel pole N3/Q2 via Hermite reduction -> second-kind U_k + elementary boundary.
- [x] FULL t-branch wrap assembled, verified 1e-14, physics-anchored (slope ~2 vs true flow).

## OPEN — off-shell closure remainder
- [ ] **Theta-nome naming of the genus-2 dilogs** (D_{j,root} at Delta=0, and the s=0 dilog in
      block A). They are isolated with symbolic coefficients but NOT yet identified with the
      canonical tabulated hyperelliptic-polylog class (Baune et al. arXiv:2306.08644 /
      D'Hoker-Schlotterer arXiv:2407.11476). This is the genuine FRONTIER piece (referee Issue 7):
      needs the explicit theta-function representation + a naming theorem, not just numerics.
      Existing partial machinery: KerrSessionScripts/kerr_tau_Wij_diffform.sage (Kleinian naming of
      the weight-1 w_ab), kerr_tau_dilog_qseries*.sage (nome q-series for the s=0 dilog).
- [ ] **Vaidya generic off-shell full closure.** The Vaidya spectral curve is S_V = r*Emu*Q2(a=0)
      -- genus-2, SAME structure as the t-branch with a=0, third-kind at the Schwarzschild horizon
      r=2m (not r_pm). The source is Theta=m d_m (NO dilation letter), so the inner letter is
      int Theta H dlambda, not int p_r dr. Mechanical port of the t-branch machinery (curve + source
      reduction + Hermite + third-kind at r=2m + assembly + physics anchor). Not frontier, just work.
- [ ] **All-parameter (M,a,J) symbolic reduction coefficients.** The coefficient TABLES
      (a_k, p_j, Q_jk, N3, R) are ALREADY symbolic in (M,a,E,J) in
      ThakurtaMetric/adiabatic_offshell_coefficients.py. The Hermite rem_k and the Delta=0 rho were
      shown E-symbolic (a,J,M fixed); the method is general but the modular inverse
      P = -2 N3 (Q2' T)^{-1} mod Q2 over QQ(M,a,E,J) is a SymPy performance wall (gcdex hangs).
      Do it with a better CAS/domain (e.g. Singular, or reduce over QQ(a,E,J)[M] tower).
- [ ] **Cosmetic: large near-cancellation** in the Hermite decomposition (block A ~ -75, block C ~
      +73, total ~ -2). Valid closed form but numerically delicate; a more natural (residue-adapted)
      decomposition would avoid the ~50x cancellation.
- [ ] **Paper update**: Tab. 1 row "Explicit off-shell special-function closed form: open" can be
      upgraded to "t-branch assembled + physics-anchored; theta-naming and Vaidya deferred". Sync
      main.tex + PRD once naming is settled.

## OPEN — TK t-branch SEPARATRIX tracking counterterm (dr_d/dlambda)  [carried over]
- [ ] The generic-J off-shell wrap develops a POWER divergence 1/(r-r_d)^2 at the physical on-path
      separatrix double root (E=1.2,a=0.9: Jc=2.9364, r_d=1.5123 in the ergosphere), residue
      proportional to DeltaS(r_d)=34.75 -- strictly worse than the on-shell log winding. The charge
      tracking dJc/dE=-0.115 moves the root (dr_d/dE=+0.0515) but does NOT cancel the dilation kernel
      power pole. Physics: the non-uniformity is LOCALIZED at r_d (clean arc matches true flow to
      floor even at 99.4% Jc; overshoot grows only near r_d and saturates) -> a moving-double-root /
      boundary-layer counterterm dr_d/dlambda is required and NOT constructed.
      Scripts: ThakurtaMetric/tk_sep_offshell_divergence.py, tk_sep_physical_separatrix.py,
      tk_sep_tracking_vs_trueflow.py. Both manuscripts App C state this honestly.
      NEXT: construct the dr_d/dlambda counterterm (differentiate the separatrix condition along the
      full dilation ray (E,J)->(E/A,J/A), get the moving-root contribution, check it subtracts the
      34.75 power pole) OR do the inner/outer matched-asymptotic (separatrix-layer) analysis.
      See memory [[offshell-third-kind-dilog]].
