# -*- coding: utf-8 -*-
# main14 memo 4.4: the dilation term -alpha*P_r in the normalized flow (eq:normflow) is a
# CONSEQUENCE of the canonical flow in the ORIGINAL variables (r, p_r) plus the chain rule
# for the time-dependent normalization P_r = p_r/A(s) -- NOT the generator of a canonical
# transformation (P_r = p_r/A is conformally symplectic: dr ^ dp_r = A dr ^ dP_r).
#
# We verify symbolically, for the factorized branch Hamiltonian
#   H_br(r, p_r; A) = A * Hbar(r, p_r/A),
# that the canonical equation dp_r/ds = -dH_br/dr, followed by P_r = p_r/A and the chain
# rule, reproduces exactly  dP_r/ds = -Hbar_r - alpha P_r  with alpha = A'/A.
import sympy as sp

r, pr, s = sp.symbols('r p_r s')
A = sp.Function('A')(s)
Hbar = sp.Function('Hbar')          # Hbar(r, P_r), abstract
Pr = pr/A
H_br = A*Hbar(r, Pr)                 # branch Hamiltonian in ORIGINAL canonical variables

# canonical Hamilton equations in (r, p_r):
dr_ds = sp.diff(H_br, pr)           # dr/ds =  dH_br/dp_r
dpr_ds = -sp.diff(H_br, r)          # dp_r/ds = -dH_br/dr   (canonical)

# pass to P_r = p_r/A and differentiate along the flow (chain rule):
dPr_ds = sp.diff(Pr, pr)*dpr_ds + sp.diff(Pr, s)

# claimed RHS of eq:normflow:  -Hbar_r|_{P_r} - alpha P_r
alpha = sp.diff(A, s)/A
X = sp.Symbol('X')
Hbar_r = sp.diff(Hbar(r, X), r).subs(X, Pr)   # d/dr Hbar at fixed P_r
claim = -Hbar_r - alpha*Pr

print("dP_r/ds (canonical + chain rule) - eq:normflow RHS =", sp.simplify(dPr_ds - claim))
print("   [expect 0 -> -alpha P_r is canonical in (r,p_r), not a canonical transf.]")
print("dH_br/dr|_{p_r} - A*Hbar_r|_{P_r}                  =",
      sp.simplify(sp.diff(H_br, r) - A*Hbar_r), "  [expect 0]")
