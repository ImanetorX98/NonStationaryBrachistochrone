# -*- coding: utf-8 -*-
# Vaidya v-branch off-shell STEP 1: identify the spectral curve D_v and split the kernel.
# H_v = pr(f-E^2) - 1 + sqrt(w) Y,  Y=sqrt(E^2 pr^2+J^2/r^2), w=E^2-f, f=1-2m/r.
# shell quadratic in pr: K2 pr^2 + 2 K1 pr + (w J^2/r^2 - 1)=0, K1=f-E^2, K2=wE^2-K1^2.
#   K1+pr K2 = +/- sqrt(D_v),  D_v = K1^2 - K2(wJ^2/r^2-1)  (discriminant/4).
# kernel d_pr G/H_pr = -wJ(K2+K1^2)/(r^2 D_v^{3/2}) +/- wJK1/(r^2 D_v)  (2nd-kind + ELEMENTARY).
# Compare D_v to the tau-branch S=r(r-2m)DE Q3.
import sympy as sp, numpy as np
r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
f=1-2*m/r; w=E**2-f; DE=(E**2-1)*r+2*m
K1=sp.simplify(f-E**2); K2=sp.simplify(w*E**2-K1**2)
D_v=sp.simplify(K1**2-K2*(w*J**2/r**2-1))
print("K1 =",sp.factor(K1),"  K2 =",sp.factor(K2))
print("D_v (v-branch spectral curve, /4 discriminant) =",sp.factor(D_v))
# tau curve
S_tau=sp.expand(r*(r-2*m)*DE*(r**2*(r-2*m)-J**2*DE))
ratio=sp.simplify(D_v/S_tau)
print("\nD_v / S_tau =",sp.factor(ratio),"  (constant/rational => same curve up to factor)")
# clear denominators: D_v as polynomial numerator
num,den=sp.fraction(sp.together(D_v))
print("\nD_v numerator (poly) =",sp.factor(num))
print("D_v denominator =",sp.factor(den))
# kernel split coefficients
A_2nd=sp.simplify(-w*J*(K2+K1**2)/r**2)     # over D_v^{3/2}
A_elem=sp.simplify(w*J*K1/r**2)             # over D_v (elementary, sign = physical branch)
print("\nkernel 2nd-kind numerator (over D_v^{3/2}) =",sp.factor(A_2nd))
print("kernel elementary numerator (over D_v) =",sp.factor(A_elem))
