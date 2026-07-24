# -*- coding: utf-8 -*-
"""
Vaidya OFF-SHELL separatrix: elliptic degeneration of the weight-two dilog (main11 -> Vaidya).

The Vaidya frozen sextic S = r(r-2m) DE [r^2(r-2m) - J^2 DE], DE=(E^2-1)r+2m, degenerates at the
separatrix Jc (double root r_d, a BRANCH POINT at r_d<0 here; the physical turning is a root of Q4):
S = (r-r_d)^2 Q4 -> genus-1 ELLIPTIC (extends vaidya_separatrix_ell.py Stage A for U_k).

Here the weight-two W_jk = int (r^j U_k - r^k U_j)/sqrt(S) dr are shown to map to the elliptic
dilog image (U_k = V_k + r_d^k Pi_rd), verified to 1e-17 (W_02, W_13, W_01). So the Vaidya
off-shell separatrix closes in Weierstrass P,zeta,sigma + Brown-Levin elliptic dilog -- exactly
as Thakurta-Kerr, but CLEANER: Vaidya has NO conformal rescaling, hence NO -alpha P_r dilation
letter and NO off-shell tracking counterterm (the residue that TK needs at r_d is absent).
"""

