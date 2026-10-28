# -*- coding: utf-8 -*-
# EXPLICIT T_alg / G_alg in elementary closed form (TK tau, genus-2).
# G_alg = 2 int (A5/sqrtS) deta + boundary,  deta=(r^3-2M r^2)/sqrtS dr.
# int (A5/sqrtS)(r^3-2M r^2)/sqrtS dr = int A5(r^3-2M r^2)/S dr  (1/sqrtS * 1/sqrtS = 1/S)
#  = RATIONAL in r (S a polynomial) -> partial fractions -> polynomial + sum log(r-r_i) over the zeros of S.
# Check: closed form vs direct integral.
import sympy as sp, numpy as np
from scipy.integrate import quad

r,E=sp.symbols('r E',positive=True)
M,a,J=sp.Rational(1),sp.Rational(9,10),sp.Rational(5,2)
Dl=r**2-2*M*r+a**2; Emu=(E**2-1)*r+2*M
S=sp.expand(r*(r-2*M)*Emu*(r*Dl-J**2*Emu)); Sp=sp.diff(S,r)
K=J*r*(r-2*M)*Emu/Dl; F=K/sp.sqrt(S)
# reduction dE F = d(A5/sqrtS) + sum c_k r^k/sqrtS  (2N = 2S A5' - A5 S' + 2 S M)
dEF=sp.diff(F,E); N=E*J*r**4*(r-2*M)**2*Emu
Ac=[sp.Symbol(f'A{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
A5=sum(Ac[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
sol=sp.solve(sp.Poly(sp.expand(2*N-(2*S*sp.diff(A5,r)-A5*Sp+2*S*Mp)),r).all_coeffs(),Ac+ck,dict=True)[0]
A5e=sp.expand(A5.subs(sol))
print("A5(r) (deg 5, coeffs rational in E) found.")

# --- elementary piece: I_el = int A5 (r^3 - 2M r^2)/S dr  (denominator S -> SIMPLE poles) ---
E0=sp.Rational(7,5)
Ae=sp.expand((A5e*(r**3-2*M*r**2)).subs(E,E0))   # numerator (poly)
Se=sp.expand(S.subs(E,E0)); Spe=sp.diff(Se,r)
# polynomial part + proper remainder
quo,rem=sp.div(sp.Poly(Ae,r),sp.Poly(Se,r))
polypart=quo.as_expr()
Ipoly=sp.integrate(polypart,r)     # integral of the polynomial (elementary)
# residues at the 6 zeros of S (numerical): res_i = numerator(r_i)/S'(r_i)  (partial fractions of the proper remainder)
Anum=sp.lambdify(r,rem.as_expr(),'numpy'); Spn=sp.lambdify(r,Spe,'numpy')
# NB the remainder rem is proper (deg<6); res_i = rem(r_i)/S'(r_i)
remc=[complex(c) for c in sp.Poly(rem.as_expr(),r).all_coeffs()]
Sc=[complex(c) for c in sp.Poly(Se,r).all_coeffs()]
roots=np.roots(Sc)
res=[np.polyval(remc,ri)/np.polyval(np.polyder(Sc),ri) for ri in roots]
print("\n=== G_alg elementary piece I_el = int A5(r^3-2M r^2)/S dr (CLOSED FORM) ===")
print("  I_el(r) = P(r) + Sum_{i=1..6} res_i * log(r - r_i)")
print("  P(r) (polynomial) =",sp.nsimplify(sp.expand(Ipoly),rational=True) if False else sp.expand(Ipoly))
print("  zeros r_i of S and residues res_i:")
for ri,re_ in zip(roots,res): print(f"    r_i={ri:+.5f}  res_i={re_:+.6f}")
Ipolyn=sp.lambdify(r,Ipoly,'numpy')
def I_el_closed(rv):
    return complex(Ipolyn(rv)+sum(res[i]*np.log(complex(rv-roots[i])) for i in range(6)))
intn=sp.lambdify(r,(Ae/Se),'numpy'); r0=11.5
print("\n  CHECK closed form vs direct integral:")
for rv in [10.0,8.0,6.5]:
    closed=(I_el_closed(rv)-I_el_closed(r0)).real
    direct=quad(lambda x: float(np.real(intn(x))),r0,rv,limit=200)[0]
    print(f"    r={rv}: closed={closed:+.8f}  direct={direct:+.8f}  diff={abs(closed-direct):.1e}")
print("\n=> I_el ELEMENTARY: polynomial P(r) + Sum res_i log(r-r_i) over the 6 zeros of S. CLOSED FORM.")
print("   G_alg = 2 I_el + [eta*(const - A5/sqrtS)] boundary term (eta=U3-2M U2 tabulated x algebraic).")
print("   => T_alg/G_alg is ELEMENTARY (polynomial+log) + explicit boundary term, NOT irreducible.")
