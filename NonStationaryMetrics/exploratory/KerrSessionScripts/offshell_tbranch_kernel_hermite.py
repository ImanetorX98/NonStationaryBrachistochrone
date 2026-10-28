# -*- coding: utf-8 -*-
# HERMITE reduction of the kernel-pole part N3/(Q2 sqrt S) (t-branch off-shell).
# S = Q2*T,  T=r*Emu.  Claim:  N3/(Q2 sqrt S) = d/dr[ P sqrt S / Q2 ] + Rem/sqrt S,
# with P a polynomial (deg<=3) and Rem polynomial, both SYMBOLIC in E. Derivation:
#   2 Rem Q2 = 2 N3 - 2 P' Q2 T - P Q2 T' + P Q2' T   => need P Q2' T = -2 N3 (mod Q2),
#   i.e. P = -2 N3 (Q2' T)^{-1} mod Q2 ; then Rem = [...]/(2 Q2) exact.
# => int N3/(Q2 sqrt S) dr = [P sqrt S / Q2] (algebraic) + sum_k rem_k U_k  (second kind).
# Folds the DOMINANT kernel-pole blocks into the U_k basis. Verify numerically.
import sympy as sp, numpy as np
from scipy.integrate import quad

r,E=sp.symbols('r E',positive=True)
M,a,J=sp.Rational(1),sp.Rational(9,10),sp.Integer(6); E0=sp.Rational(7,5)
Emu=(E**2-1)*r+2*M; Dl=r**2-2*M*r+a**2
Q2=(2*E**2*J**2*M*r-E**2*J**2*r**2-4*E**2*J*M*a*r+2*E**2*M*a**2*r+E**2*a**2*r**2
    +E**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2-8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
T=sp.expand(r*Emu); S=sp.expand(Q2*T)
Qp,N3=sp.div(sp.Poly(sp.expand(E**2*J*r**4*Emu),r),sp.Poly(sp.expand(Q2),r)); N3=N3.as_expr()
Q2p=sp.diff(Q2,r); Tp=sp.diff(T,r)

# P = -2 N3 (Q2' T)^{-1}  mod Q2
Q2poly=sp.Poly(Q2,r)
inv=sp.invert(sp.Poly(sp.expand(Q2p*T),r), Q2poly)          # (Q2' T)^{-1} mod Q2
P=sp.rem(sp.expand(-2*N3*inv.as_expr()), Q2, r); P=sp.expand(P)
# Rem = (2 N3 - 2 P' Q2 T - P Q2 T' + P Q2' T)/(2 Q2)   -- must be exact polynomial
num=sp.expand(2*N3 - 2*sp.diff(P,r)*Q2*T - P*Q2*Tp + P*Q2p*T)
Rem,rr=sp.div(sp.Poly(num,r), sp.Poly(sp.expand(2*Q2),r))
print("exact division remainder (must be 0):", sp.simplify(rr.as_expr()))
Rem=Rem.as_expr()
rem_k=[sp.simplify(c) for c in sp.Poly(Rem,r).all_coeffs()[::-1]]
print("P(r) deg =",sp.degree(sp.Poly(P,r)),"  Rem deg =",sp.degree(sp.Poly(Rem,r)))
print("rem_k (symbolic, second-kind U_k coeffs of the kernel-pole part):")
for k,c in enumerate(rem_k): print(f"  rem_{k} =",c)

# ---- numeric verification: int N3/(Q2 sqrt S) dr == [P sqrt S/Q2] + sum rem_k U_k ----
Sn=sp.lambdify(r,S.subs(E,E0),'numpy'); sqn=lambda x:np.sqrt(max(Sn(x),0.0))
Q2n=sp.lambdify(r,Q2.subs(E,E0),'numpy'); N3n=sp.lambdify(r,N3.subs(E,E0),'numpy')
Pn=sp.lambdify(r,P.subs(E,E0),'numpy'); remn=[float(c.subs(E,E0)) for c in rem_k]
r0=11.0
def Uk(x,k): return quad(lambda t:t**k/sqn(t),r0,x,limit=200)[0]
def PhiN_direct(x): return quad(lambda t:N3n(t)/(Q2n(t)*sqn(t)),r0,x,limit=200)[0]
def alg(x): return Pn(x)*sqn(x)/Q2n(x)
def PhiN_asm(x): return (alg(x)-alg(r0)) + sum(remn[k]*Uk(x,k) for k in range(len(remn)))
print("\nverify Phi_N = int N3/(Q2 sqrtS) dr  vs  [P sqrtS/Q2] + sum rem_k U_k:")
for x in [10.0,8.0,6.5]:
    d=PhiN_direct(x); asm=PhiN_asm(x)
    print(f"  r={x}: direct={d:+.8f}  assembly={asm:+.8f}  diff={abs(d-asm):.1e}")
print("\n=> kernel-pole part reduced to second-kind U_k + algebraic boundary, SYMBOLIC in E.")
print("   => full kernel A/sqrtS = sum_k (p_k+rem_k) U_k-differential + d[P sqrtS/Q2].")
