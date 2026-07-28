# -*- coding: utf-8 -*-
# Vaidya tau-branch off-shell correction (CLEAN non-surrogate Hamiltonian, no "-1" cost).
# Source Theta=m d_m (metric modulus, NO dilation). Off-shell wrap:
#   Phi = -int (d_pr G / H_pr) * Sigma dr,  Sigma(r)=int (m H_m / H_pr) dr  (inner letter).
# Clean shell: H_tau=0 => P_r0=sqrt(S)/(r(r-2m)DE), S=r(r-2m)DE*Q3, Q3=r^2(r-2m)-J^2 DE, DE=(E^2-1)r+2m.
# STEP 1: build H_tau, verify clean P_r0, kernel=A/sqrt(S), reduce inner letter to poly/sqrtS + 3rd kind @ r=2m.
import sympy as sp, numpy as np
from scipy.integrate import quad

r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
Dl=r*(r-2*m); f=(r-2*m)/r; DE=(E**2-1)*r+2*m; v=DE/(E**2*r)   # a=0 Vaidya/Schwarzschild
# clean tau Hamiltonian (a=0: frame-drag b=0, ptilde=J):
Htau=sp.sqrt(Dl*v/r**2)*sp.sqrt((Dl/r**2)*pr**2+J**2/r**2)-f/E
# on-shell P_r0
pr0=sp.sqrt(r*(r-2*m)*DE*(r**2*(r-2*m)-J**2*DE))/(r*(r-2*m)*DE)
S=sp.expand(r*(r-2*m)*DE*(r**2*(r-2*m)-J**2*DE))
chk=sp.simplify(Htau.subs(pr,pr0))
print("H_tau(P_r0) on shell =",chk," (=0 => clean P_r0 verified)")

Hpr=sp.diff(Htau,pr); Hm=sp.diff(Htau,m); G=sp.diff(Htau,J)/Hpr; Gpr=sp.diff(G,pr)
# on-shell substitutions
def onshell(expr): return sp.simplify(expr.subs(pr,pr0))
kernel=onshell(Gpr/Hpr)          # should be A/sqrt(S)
inner=onshell(m*Hm/Hpr)          # inner-letter integrand (Theta H / H_pr)
print("\nkernel (d_pr G/H_pr) on shell:"); print("  ",sp.simplify(kernel*sp.sqrt(S)),"/ sqrt(S)")
print("\ninner letter integrand (m H_m/H_pr) on shell:"); print("  ",sp.nsimplify(sp.simplify(inner)))
# structure of inner: rational; split poly/sqrtS + third-kind (poles at r=2m and DE=0)
print("\ninner as (rational): factor ="); print("  ",sp.factor(inner))
