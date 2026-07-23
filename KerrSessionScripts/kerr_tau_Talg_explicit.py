# -*- coding: utf-8 -*-
# T_alg / G_alg ESPLICITO in forma chiusa elementare (TK tau, genus-2).
# G_alg = 2 int (A5/sqrtS) deta + boundary,  deta=(r^3-2M r^2)/sqrtS dr.
# int (A5/sqrtS)(r^3-2M r^2)/sqrtS dr = int A5(r^3-2M r^2)/S dr  (1/sqrtS * 1/sqrtS = 1/S)
#  = RAZIONALE in r (S polinomio) -> partial fractions -> polinomio + sum log(r-r_i) su zeri di S.
# Verifica: forma chiusa vs integrale diretto.
import sympy as sp, numpy as np
from scipy.integrate import quad

r,E=sp.symbols('r E',positive=True)
M,a,J=sp.Rational(1),sp.Rational(9,10),sp.Rational(5,2)
Dl=r**2-2*M*r+a**2; Emu=(E**2-1)*r+2*M
S=sp.expand(r*(r-2*M)*Emu*(r*Dl-J**2*Emu)); Sp=sp.diff(S,r)
K=J*r*(r-2*M)*Emu/Dl; F=K/sp.sqrt(S)
# riduzione dE F = d(A5/sqrtS) + sum c_k r^k/sqrtS  (2N = 2S A5' - A5 S' + 2 S M)
dEF=sp.diff(F,E); N=E*J*r**4*(r-2*M)**2*Emu
Ac=[sp.Symbol(f'A{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
A5=sum(Ac[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
sol=sp.solve(sp.Poly(sp.expand(2*N-(2*S*sp.diff(A5,r)-A5*Sp+2*S*Mp)),r).all_coeffs(),Ac+ck,dict=True)[0]
A5e=sp.expand(A5.subs(sol))
print("A5(r) (deg 5, coeff razionali in E) trovato.")

# --- pezzo elementare: I_el = int A5 (r^3 - 2M r^2)/S dr  (denom S -> poli SEMPLICI) ---
E0=sp.Rational(7,5)
Ae=sp.expand((A5e*(r**3-2*M*r**2)).subs(E,E0))   # numeratore (poly)
Se=sp.expand(S.subs(E,E0)); Spe=sp.diff(Se,r)
# parte polinomiale + resto proprio
quo,rem=sp.div(sp.Poly(Ae,r),sp.Poly(Se,r))
polypart=quo.as_expr()
Ipoly=sp.integrate(polypart,r)     # integrale polinomio (elementare)
# residui ai 6 zeri di S (numerici): res_i = numeratore(r_i)/S'(r_i)  (partial fraction del resto proprio)
Anum=sp.lambdify(r,rem.as_expr(),'numpy'); Spn=sp.lambdify(r,Spe,'numpy')
# NB il resto rem e' proprio (deg<6); res_i = rem(r_i)/S'(r_i)
remc=[complex(c) for c in sp.Poly(rem.as_expr(),r).all_coeffs()]
Sc=[complex(c) for c in sp.Poly(Se,r).all_coeffs()]
roots=np.roots(Sc)
res=[np.polyval(remc,ri)/np.polyval(np.polyder(Sc),ri) for ri in roots]
print("\n=== G_alg pezzo elementare I_el = int A5(r^3-2M r^2)/S dr (FORMA CHIUSA) ===")
print("  I_el(r) = P(r) + Sum_{i=1..6} res_i * log(r - r_i)")
print("  P(r) (polinomio) =",sp.nsimplify(sp.expand(Ipoly),rational=True) if False else sp.expand(Ipoly))
print("  zeri r_i di S e residui res_i:")
for ri,re_ in zip(roots,res): print(f"    r_i={ri:+.5f}  res_i={re_:+.6f}")
Ipolyn=sp.lambdify(r,Ipoly,'numpy')
def I_el_closed(rv):
    return complex(Ipolyn(rv)+sum(res[i]*np.log(complex(rv-roots[i])) for i in range(6)))
intn=sp.lambdify(r,(Ae/Se),'numpy'); r0=11.5
print("\n  VERIFICA forma chiusa vs integrale diretto:")
for rv in [10.0,8.0,6.5]:
    closed=(I_el_closed(rv)-I_el_closed(r0)).real
    direct=quad(lambda x: float(np.real(intn(x))),r0,rv,limit=200)[0]
    print(f"    r={rv}: chiuso={closed:+.8f}  diretto={direct:+.8f}  diff={abs(closed-direct):.1e}")
print("\n=> I_el ELEMENTARE: polinomio P(r) + Sum res_i log(r-r_i) sui 6 zeri di S. FORMA CHIUSA.")
print("   G_alg = 2 I_el + [eta*(cost - A5/sqrtS)] boundary (eta=U3-2M U2 tabulato x algebrico).")
print("   => T_alg/G_alg e' ELEMENTARE (polinomio+log) + boundary esplicito, NON irriducibile.")
