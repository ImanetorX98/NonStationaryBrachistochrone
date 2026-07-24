# -*- coding: utf-8 -*-
"""
SETUP of the off-shell extension of the separatrix-tracking theorem (main11 -> separatrix).

The on-shell tracking (tk_t_sep_track_blockassembly.py) cancels the triple pole of the
adiabatic source via N_tot = N_t + (dJc/dE) N_J, dJc/dE = -Q2_E/Q2_J at (r_d, Jc): the
charge-tracking follows the moving double root in the (E,J) plane.

The main11 off-shell source is D H = (E d_E + J d_J + P_r d_Pr) H. The NEW piece P_r d_Pr
(the conformal dilation letter, = int p_r dr) injects a triple-pole residue at the separatrix
double root r_d:
    kernel  A/sqrt(S) ~ c3/(r-r_d)^3  (A ~ 1/Q2 double pole, sqrt(S) ~ (r-r_d)),
    source dilation part -> Delta S(r_d) (finite radial action to the separatrix),
    residue ~ c3 * Delta S(r_d)  =  14.39 != 0  (E=1.2, a=0.9, Jc=2.936, r_d~1.512, inside ergosphere).
The charge-only tracking dJc/dE does NOT cancel this (it acts on E,J, not on P_r). Hence the
tracking theorem must be EXTENDED with a term tied to the moving double root dr_d/dlambda to
cancel the dilation triple-pole residue -- OR the residual survives and the separatrix
first-order form gains a new term. Vaidya is immune (no rescaling, no dilation letter).

This module: computes c3, Delta S(r_d) and the residue for a t-branch (ergosphere) separatrix.
"""

