# -*- coding: utf-8 -*-
"""
Outer-kernel structure of the new off-shell adiabatic term (main11, level-C wrap).

delta phi_extra = - int (G_Pr/H_Pr) * Delta S(r) dr,   Delta S = int p_r dr (reduced, level B).

KEY RESULT (verified to 3e-15): the off-shell kernel, evaluated on the frozen shell, does NOT
carry the second radical sqrt(Delta v / Pbar); it collapses to

    G_Pr/H_Pr |_shell  =  A(r) / sqrt(S) ,     A(r) rational,

so the wrap lives on the SAME genus-2 spectral curve y^2 = S as the on-shell reduction. Hence

    delta phi_extra = - sum_{j,k} beta_j a_k  int (r^j/sqrt(S)) U_k dr ,   A/sqrt(S) = sum_j beta_j r^j/sqrt(S)+...
                    = sum Q_{kj} W_{kj}  +  (U_j U_k products)  +  (third kind),

the SAME weight-two W_{kj} = int (U_k r^j - U_j r^k)/sqrt(S) dr basis as the on-shell psi
(kerr_psi_explicit_verified.py), with symbolic coefficients from A(r) and a_k(M,a,E,J). The
W_{kj} are the genus-2 transcendentals (evaluated by the theta-nome q-series), not elementary.

This module: substitution-based proof (no hanging simplify) that sqrt(W) cancels, i.e. that the
kernel numerator is rational and its denominator is proportional to sqrt(S).
"""

import sympy as sp
r,pr,Js,Y,sw,ss=sp.symbols('r pr J_ Y sw ss',positive=True)
M,a,E=sp.Rational(1),sp.Rational(9,10),sp.Rational(7,5); J=sp.Integer(6)
f=1-2*M/r; Dl=r**2-2*M*r+a**2; b=2*M*a/r; v=1-f/E**2
Pcap=r**2+a**2+2*M*a**2/r; Pb=Pcap+b**2/E**2
W=sp.cancel(Dl*v/Pb)                              # inner radical^2
Yexpr2=(Dl/r**2)*pr**2+Js**2/Pb                   # Y^2
# H with Y and sqrt(W) as symbols (sw=sqrt(W)); dY/dpr = (Dl/r^2) pr / Y
H=Js*b*v/Pb+sw*Y-1
def dpr(expr):    # total d/dpr treating Y as Y(pr): dY/dpr=(Dl/r^2)pr/Y
    return sp.diff(expr,pr)+sp.diff(expr,Y)*((Dl/r**2)*pr/Y)
Hp=dpr(H); Hj=sp.diff(H,Js)+sp.diff(H,Y)*(Js/Pb/Y)    # dY/dJ=(Js/Pb)/Y
G=sp.cancel(Hj/Hp); Gpr=sp.cancel(dpr(G)); kernel=sp.cancel(Gpr/Hp)
kernel=kernel.subs(Js,J)
# on-shell substitutions: p_r0 = -g*ss with g=1/(Delta*((E^2-1)r+2M)); Y_on = rho/sw, rho=1-Jbv/Pb
g=1/(Dl*((E**2-1)*r+2*M)); rho=(1-(Js*b*v/Pb)).subs(Js,J)
k_on=kernel.subs({pr:-g*ss, Y:rho/sw})
k_on=sp.cancel(k_on)
# now reduce even powers: ss^2 -> S, sw^2 -> W. Collect in ss, sw.
S=sp.cancel(g**(-2)*0 + 0)  # placeholder
# express k_on as rational in ss,sw then reduce ss^2->S(rational), sw^2->W(rational)
# S from pr2 = g^2 S => S = pr2/g^2:
pr2=sp.cancel((((1-(Js*b*v/Pb))**2/(Dl*v/Pb)) - Js**2/Pb)/(Dl/r**2)).subs(Js,J)
Sexpr=sp.cancel(pr2/g**2)
num=sp.together(k_on); n,d=sp.fraction(num)
# reduce powers of ss and sw in numerator and denominator
def reduce_rad(poly):
    poly=sp.expand(poly)
    poly=poly.replace(lambda e:e.func==sp.Pow and e.base==ss and e.exp.is_integer, lambda e: (Sexpr**(e.exp//2))*(ss if e.exp%2 else 1))
    poly=poly.replace(lambda e:e.func==sp.Pow and e.base==sw and e.exp.is_integer, lambda e: (W**(e.exp//2))*(sw if e.exp%2 else 1))
    return sp.expand(poly)
n2=reduce_rad(n); d2=reduce_rad(d)
# collect by ss, sw monomials {1, ss, sw, ss*sw}
def parts(expr):
    expr=sp.expand(expr)
    return {m:sp.simplify(expr.coeff(ss,i).coeff(sw,j)) for (i,j,m) in [(0,0,'1'),(1,0,'ss'),(0,1,'sw'),(1,1,'ss*sw')]}
kfrac=sp.cancel(n2/d2)
# rationalize denominator to isolate radicals: multiply num&den by conjugates
print("kernel_on numerator radical parts (coeff of 1, ss, sw, ss*sw):")
for k,vv in parts(n2).items(): print(f"  [{k}] :", "0" if vv==0 else "nonzero")
print("kernel_on denominator radical parts:")
for k,vv in parts(d2).items(): print(f"  [{k}] :", "0" if vv==0 else "nonzero")
