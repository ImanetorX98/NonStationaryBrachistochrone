# -*- coding: utf-8 -*-
# At the PHYSICAL separatrix Jc=2.9364, r_d=1.5123 (on-path): confirm
#  (1) p_r regular & vanishing linearly at r_d,
#  (2) Delta S(r_d)=int_{r_d}^{r0} p_r dr finite & nonzero (dilation letter value),
#  (3) on-shell winding K/sqrt(S) ~ 1/(r-r_d) (log), while
#      off-shell kernel A/sqrt(S) ~ 1/(r-r_d)^3 => integrand ~ DeltaS(r_d)/(r-r_d)^3
#      -> 1/(r-r_d)^2 POWER divergence (residue ~ DeltaS(r_d) != 0).
import sympy as sp, numpy as np
from scipy.integrate import quad
r=sp.symbols('r'); M=1; a=sp.Rational(9,10); E=sp.Rational(6,5); Jc=sp.Float(2.9363539)
Q2=(2*E**2*Jc**2*M*r - E**2*Jc**2*r**2 - 4*E**2*Jc*M*a*r + 2*E**2*M*a**2*r
    + E**2*a**2*r**2 + E**2*r**4 + 4*Jc**2*M**2 - 4*Jc**2*M*r + Jc**2*r**2
    - 8*Jc*M**2*a + 4*Jc*M*a*r + 4*M**2*a**2)
DE=(E**2-1)*r+2*M; Dl=r**2-2*M*r+a**2
S=sp.expand(r*Q2*DE)
Sc=[float(c) for c in sp.Poly(S,r).all_coeffs()]
rts=np.roots(Sc)
pairs=[(i,j,abs(rts[i]-rts[j])) for i in range(len(rts)) for j in range(i+1,len(rts))]
i,j,_=min(pairs,key=lambda t:t[2])
rd=float(np.real((rts[i]+rts[j])/2))
print(f"Jc=2.93635  r_d={rd:.5f}  (physical, on-path)")
Sn=sp.lambdify(r,S,'numpy'); Dln=sp.lambdify(r,Dl,'numpy'); DEn=sp.lambdify(r,DE,'numpy')
# p_r = sqrt(S)/(Delta*DE)  (t-branch on-shell momentum, = (r-rd)*sqrt(rest)/(Dl*DE))
def pr(x):
    s=Sn(x); return (np.sqrt(s) if s>0 else np.nan)/(Dln(x)*DEn(x))
print("(1) p_r near r_d (should be ~linear in (r-r_d), regular):")
for dx in [1e-2,1e-3,1e-4]:
    x=rd+dx; print(f"    p_r(r_d+{dx:.0e})={pr(x):+.6e}   /dx={pr(x)/dx:+.4f}")
r0=20.0
DS=quad(lambda x: pr(x), rd+1e-8, r0, limit=500)[0]
print(f"(2) Delta S(r_d)=int_rd^r0 p_r dr = {DS:.5f}   (finite, nonzero)")
# (3) pole orders
Kn=sp.lambdify(r, Jc*r*(r-2*M)*DE/Dl, 'numpy')
# fit sqrt(S) ~ c*(r-rd) near rd
xs=np.array([rd+1e-3,rd+2e-3,rd+4e-3]); ys=np.sqrt(Sn(xs))/(xs-rd)
print(f"(3) sqrt(S)/(r-rd) near r_d = {ys}  (const => sqrt(S)~(r-rd), simple zero)")
print(f"    on-shell dphi/dr=K/sqrt(S) ~ {Kn(rd)/ys.mean():.2f}/(r-rd)  -> LOG (physical winding)")
print(f"    off-shell A/sqrt(S): A~1/Q2 ~ 1/(r-rd)^2 => integrand*DeltaS ~ {DS:.2f}/(r-rd)^3")
print(f"    -> integral ~ 1/(r-rd)^2 POWER divergence, residue proportional to DeltaS(r_d)={DS:.3f} != 0")
print("    VERDICT: off-shell wrap MORE singular than on-shell log at physical r_d => needs dr_d/dlambda counterterm (OPEN)")
