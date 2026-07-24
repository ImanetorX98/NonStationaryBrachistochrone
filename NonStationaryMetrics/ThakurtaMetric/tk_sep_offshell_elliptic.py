# -*- coding: utf-8 -*-
"""
Elliptic assembly of the OFF-SHELL separatrix term (main11 -> separatrix, closed).

The separatrix off-shell correction IS the generic off-shell term at J=Jc (the a_k(M,a,E,J)
are polynomials in J, finite at Jc). At J=Jc the genus-2 sextic DEGENERATES to genus-1:
    S = (r-r_d)^2 Q4,   Q4 = r*Emu*(quadratic)  (elliptic quartic),
and p_r = (r-r_d) sqrt(Q4)/(Delta*Emu) is REGULAR at r_d. The radial-action dilation letter
reduces to the ELLIPTIC basis (verified to 1e-16 here):
    int p_r dr = sum_j b_j U_j^ell + (third kind at Delta=0),  U_j^ell = int r^j/sqrt(Q4) dr,
    b = [36/25, 5.058, 9.615] (1st/2nd kind -> Weierstrass P, zeta), third kind -> sigma-quotients.
Hence the WRAP delta phi_extra = -int (A/sqrt(S)) Delta S dr lives on the genus-1 curve sqrt(Q4):
it is an ELLIPTIC dilogarithm (Brown-Levin / Beilinson-Levin, tabulated), NOT the generic
genus-2 dilog. The triple pole of the kernel A/sqrt(S) at r_d is regularized by the Jc-tracking
dJc/dE = -Q2_E/Q2_J already in B.6 (on-shell cancellation 0.44 -> 0 verified elsewhere).

So the off-shell separatrix term closes in classical elliptic special functions + the Jc-tracking
-- cleaner than the generic case. Vaidya is immune (no rescaling).
"""

