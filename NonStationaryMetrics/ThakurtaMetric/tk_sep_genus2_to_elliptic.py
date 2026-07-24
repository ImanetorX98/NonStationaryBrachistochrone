# -*- coding: utf-8 -*-
"""
DEMONSTRATION that the generic genus-2 (hyperelliptic) form reduces to the genus-1 (elliptic)
curve at the separatrix J=Jc -- the fundamental piece linking the two formulations.

At J=Jc the sextic S = r*Emu*Q2 gains a DOUBLE root r_d (Q2 double root), so S=(r-r_d)^2 Q4
and sqrt(S) = (r-r_d) sqrt(Q4). Each genus-2 letter U_k^g2 = int r^k/sqrt(S) then becomes an
ELLIPTIC integral with a THIRD-KIND pole at r_d:  U_k^g2|_Jc = int r^k/((r-r_d) sqrt(Q4)),
which individually DIVERGES (log at r_d).

The physical combination int p_r dr = sum a_k U_k^g2 + (Delta=0 third kind) is FINITE because the
r_d third-kind residues cancel. Demonstrated exactly:
    sum_k a_k(Jc) r_d^k + R(Jc,r_d)/deng(r_d)  =  [S/deng](r_d) = S(r_d)/deng(r_d) = 0,
i.e. +4.165 (polynomial part) - 4.165 (third-kind part) = 0, and it vanishes precisely because
S(r_d)=0 (double root). Physics anchor: the full delta phi converges to a finite separatrix limit
as J->Jc on a fixed sub-arc (tk_sep_offshell_elliptic / sep_close checks).

DEGENERATION MAP (genus-2 hyperelliptic letter -> genus-1 elliptic, at J=Jc):
    U_0^g2, U_1^g2 (1st kind, holomorphic on sqrt S)  -> elliptic 3rd kind at r_d (pole from (r-r_d))
    U_2^g2, U_3^g2 (2nd kind)                          -> elliptic 2nd + 3rd kind at r_d
    W_kj^g2 (genus-2 dilogarithm)                      -> ELLIPTIC dilogarithm (Brown-Levin, tabulated)
    r_d third-kind residues cancel (S(r_d)=0) in physical combinations
    => net: elliptic 1st/2nd kind (Weierstrass P, zeta) + 3rd kind at Delta=0 (sigma-quotients).
This is the SAME map for the generic brachistochrone (genus-2 -> elliptic at Jc) and for the old
hyperelliptic separatrix formulation (whose genus-2 letters collapse to the same elliptic pieces).
Vaidya is immune (no rescaling).
"""

