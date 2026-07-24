# -*- coding: utf-8 -*-
"""
Elliptic assembly of the OFF-SHELL separatrix term (main11 -> separatrix).

STATUS: the elliptic BASIS is closed; the separatrix TRACKING of the dilation pole is OPEN.
(This file previously OVER-CLAIMED that dJc/dE closes it. Corrected below and cross-checked
against tk_sep_offshell_tracking_setup.py and tk_sep_offshell_divergence.py.)

The separatrix off-shell correction IS the generic off-shell term at J=Jc (the a_k(M,a,E,J)
are polynomials in J, finite at Jc). The PHYSICAL t-branch separatrix comes from the GENUINE
branch quartic Q2 (kerr/wrap_assembly.py), NOT the simplified r*Delta-J^2*DE:
    E=1.2, a=0.9  ->  Jc=2.93635,  r_d=+1.51229   (ON the physical arc [r+=1.436, r0]).
At J=Jc the genus-2 sextic S=r*Emu*Q2 DEGENERATES: Q2=(r-r_d)^2*(quadratic), S=(r-r_d)^2 Q4,
Q4=r*Emu*(quadratic). Each genus-2 letter U_k -> elliptic (Weierstrass P, zeta, sigma), and the
weight-two W_jk -> ELLIPTIC dilogarithm (Brown-Levin, tabulated). So AWAY from r_d the off-shell
term closes in classical elliptic special functions (verified against the true flow on sub-arcs:
slope ~2.1, residual ~4e-5).

THE OPEN PIECE (separatrix TRACKING of the dilation pole), verified in tk_sep_offshell_divergence.py:
  - p_r = (r-r_d) sqrt(Q4)/(Delta*Emu) is regular and VANISHES linearly at r_d, but the dilation
    letter value DeltaS(r_d) = int_{r_d}^{r0} p_r dr = 34.75 is FINITE and NONZERO.
  - the wrap kernel A/sqrt(S) ~ 1/(r-r_d)^3 (A ~ 1/Q2 double pole, sqrt(S) ~ (r-r_d)), so the
    off-shell integrand ~ DeltaS(r_d)/(r-r_d)^3 -> the wrap has a POWER divergence 1/(r-r_d)^2 at
    r_d, STRICTLY WORSE than the on-shell log winding K/sqrt(S) ~ 1/(r-r_d).
  - the charge tracking dJc/dE=-0.115 DOES move the double root (dr_d/dE=+0.0515 on the family),
    but it was built to cancel the on-shell SOURCE pole (N_tau,N_J), not the dilation KERNEL pole.
    Whether the FULL tracked total derivative (charge shift dJc/dE + root shift dr_d/dlambda)
    cancels the residue-34.75 power pole is NOT constructed here.

VERDICT: off-shell separatrix term = elliptic dilogarithm basis (closed) + an OPEN dr_d/dlambda
tracking counterterm for the dilation power-pole at the on-path double root. Physically this is the
non-uniformity of first-order adiabatic theory across the separatrix (power vs log singularity),
requiring a moving-double-root (boundary-layer) tracking term, not a plain algebraic subtraction.
Vaidya is immune (no conformal rescaling => no dilation letter => no dilation pole).
"""
