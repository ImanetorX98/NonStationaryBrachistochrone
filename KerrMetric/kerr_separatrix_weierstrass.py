# -*- coding: utf-8 -*-
"""
phi(r) in forma chiusa per la SEPARATRICE J = J_c = a/E della
brachistocrona tau equatoriale di Kerr (l'orbita penetrante, doranTau).

Dalla fattorizzazione doranTau  Dl - J_c^2 w = f (r^2 + c^2),  c = a/E:

    dphi/dr = (a/E) sqrt(Q4(r)) / ( Dl(r) P2(r) )
    Q4(r) = r (2M + (E^2-1) r) (r^2 + c^2)      [QUARTICA => ellittico]
    P2(r) = r^2 + c^2 ,   Dl = (r - r_+)(r - r_-)

PASSI:
 1. frazioni parziali:  Q4/(Dl P2) = (E^2-1) + sum_k alpha_k/(r - c_k),
    poli c_k in {r_+, r_-, +ic, -ic}  (residui espliciti)
 2. uniformizzazione di Weierstrass dalla radice r=0 di Q4:
       r(z) = A/(P(z) - B),   A = Q4'(0)/4,  B = Q4''(0)/24
       z(r) = P^{-1}(A/r + B),  invarianti g2, g3 dalla quartica
    [verifica simbolica: (dP/dz)^2 = 4P^3 - g2 P - g3 IDENTICAMENTE]
 3. integrali di terza specie (Whittaker-Watson 20.53):
       int dz/(P(z) - P(v)) = (1/P'(v)) [ 2 zeta(v) z + ln sigma(z-v)/sigma(z+v) ]

FORMA FINALE:

    phi(r) = (a/E) [ Lam0 * z(r)
                     + sum_k lam_k * ln( sigma(z - v_k) / sigma(z + v_k) ) ]

    P(v_k) = A/c_k + B
    lam_k  = - A alpha_k / ( c_k^2 P'(v_k) )
    Lam0   = (E^2-1) - sum_k alpha_k/c_k + sum_k 2 zeta(v_k) lam_k

(z lineare + 4 rapporti di sigma; la coppia complessa +-ic si combina in
forma reale). Equivalente Legendre: 1 F + 3 Pi.
"""

import numpy as np
import sympy as sp

r, M, E, a = sp.symbols('r M E a', positive=True)
c = a / E
f = 1 - 2 * M / r
w = E**2 - f
Dl = r**2 - 2 * M * r + a**2
P2 = r**2 + c**2
Q4 = sp.expand(r * (2 * M + (E**2 - 1) * r) * (r**2 + c**2))

print("=" * 72)
print("[1] separatrice e quartica")
print("=" * 72)
# dphi/dr = J_c r sqrt(w f)/(Dl sqrt(Dl - J_c^2 w)) con la fattorizzazione
lhs = (a / E) * r * sp.sqrt(w * f) / (Dl * sp.sqrt(f * P2))
rhs = (a / E) * sp.sqrt(Q4) / (Dl * P2)
chk = sp.simplify((lhs / rhs)**2)
print("  (dphi/dr)^2 / [(a/E)^2 Q4/(Dl P2)^2] =", chk)
print("  Q4 =", sp.factor(Q4))

print()
print("=" * 72)
print("[2] frazioni parziali con residui espliciti")
print("=" * 72)
rp, rm = sp.symbols('r_+ r_-', positive=True)     # radici di Dl
ic = sp.I * c
alpha_p = Q4.subs(r, rp) / ((rp - rm) * P2.subs(r, rp))
alpha_m = Q4.subs(r, rm) / ((rm - rp) * P2.subs(r, rm))
alpha_ic = Q4.subs(r, ic) / (Dl.subs(r, ic) * 2 * ic)
alpha_mic = Q4.subs(r, -ic) / (Dl.subs(r, -ic) * (-2 * ic))
ricomp = ((E**2 - 1) + alpha_p / (r - rp) + alpha_m / (r - rm)
          + alpha_ic / (r - ic) + alpha_mic / (r + ic))
resto = sp.simplify(sp.together(Q4 / (Dl * P2) - ricomp)
                    .subs([(rp * rm, a**2), (rp + rm, 2 * M)]))
# verifica numerica robusta (l'identita' simbolica con radici vincolate)
num_ok = True
rng = np.random.default_rng(2)
for _ in range(6):
    Mv, av, Ev = 1.0, rng.uniform(0.3, 0.95), rng.uniform(1.05, 1.6)
    rpv = Mv + np.sqrt(Mv**2 - av**2)
    rmv = Mv - np.sqrt(Mv**2 - av**2)
    rv = rng.uniform(2.5, 12.0)
    sub = {M: Mv, a: av, E: Ev, r: rv, rp: rpv, rm: rmv}
    lhs_n = complex((Q4 / (Dl * P2)).subs(sub))
    rhs_n = complex(ricomp.subs(sub))
    num_ok = num_ok and abs(lhs_n - rhs_n) < 1e-10
print("  Q4/(Dl P2) = (E^2-1) + somma residui:", num_ok)
print("  alpha_+ =", sp.simplify(alpha_p))
print("  alpha_(+ic) =", sp.simplify(alpha_ic))

print()
print("=" * 72)
print("[3] uniformizzazione di Weierstrass dalla radice r = 0")
print("=" * 72)
A_ = sp.diff(Q4, r).subs(r, 0) / 4
B_ = sp.diff(Q4, r, 2).subs(r, 0) / 24
print("  A = Q4'(0)/4  =", sp.simplify(A_))
print("  B = Q4''(0)/24 =", sp.simplify(B_))
# invarianti della quartica a0 r^4 + 4a1 r^3 + 6a2 r^2 + 4a3 r + a4
a0 = sp.Poly(Q4, r).coeff_monomial(r**4)
a1 = sp.Poly(Q4, r).coeff_monomial(r**3) / 4
a2 = sp.Poly(Q4, r).coeff_monomial(r**2) / 6
a3 = sp.Poly(Q4, r).coeff_monomial(r) / 4
a4 = sp.Integer(0)
g2 = sp.simplify(a0 * a4 - 4 * a1 * a3 + 3 * a2**2)
g3 = sp.simplify(a0 * a2 * a4 + 2 * a1 * a2 * a3 - a2**3 - a0 * a3**2
                 - a1**2 * a4)
print("  g2 =", g2)
print("  g3 =", g3)
# VERIFICA: P(z) = A/r + B soddisfa (P')^2 = 4P^3 - g2 P - g3
Pw = A_ / r + B_
lhs_w = (A_ / r**2)**2 * Q4          # (dP/dz)^2 = (A/r^2)^2 (dr/dz)^2, (dr/dz)^2=Q4
rhs_w = 4 * Pw**3 - g2 * Pw - g3
print("  (P')^2 - (4P^3 - g2 P - g3) =",
      sp.simplify(sp.expand(lhs_w - rhs_w)))

print()
print("=" * 72)
print("[4] FORMA FINALE (Weierstrass sigma/zeta)")
print("=" * 72)
print("""
  z(r) = P^{-1}( A/r + B ; g2, g3 )        [P = Weierstrass p]

  phi(r) = (a/E) [ Lam0 * z(r)
           + SUM_k lam_k * ln( sigma(z(r) - v_k) / sigma(z(r) + v_k) ) ]
           + cost

  poli:      c_k in { r_+, r_-, +i a/E, -i a/E }
  immagini:  P(v_k) = A/c_k + B
  pesi:      lam_k  = - A alpha_k / ( c_k^2 P'(v_k) )
  lineare:   Lam0   = (E^2-1) - SUM_k alpha_k/c_k + SUM_k 2 zeta(v_k) lam_k

  con A, B, g2, g3, alpha_k espliciti sopra. I termini k = +-i a/E sono
  complessi coniugati: la loro somma e' reale (2 Re[...]).
  Equivalente in forma di Legendre: 1 integrale di prima specie F
  + 3 di terza specie Pi (caratteristiche r_+, r_-, coppia +-ic).
""")

print("=" * 72)
print("[5] validazione numerica: quadratura in r vs quadratura in z")
print("=" * 72)
# cambio di variabili + decomposizione (tutto tranne l'antiderivata
# classica sigma/zeta, che e' Whittaker-Watson 20.53)
from scipy.integrate import quad
Mv, av, Ev = 1.0, 0.9, 1.2
cv = av / Ev
rpv = Mv + np.sqrt(Mv**2 - av**2)
rmv = Mv - np.sqrt(Mv**2 - av**2)
Q4n = sp.lambdify(r, Q4.subs([(M, Mv), (a, av), (E, Ev)]), 'numpy')
integ_r = lambda rv: (av / Ev) * np.sqrt(Q4n(rv)) / (
    (rv**2 - 2 * Mv * rv + av**2) * (rv**2 + cv**2))
r_a, r_b = 3.0, 9.0
I_r = quad(integ_r, r_a, r_b, limit=200)[0]

al = {}
for ck in (rpv, rmv, 1j * cv, -1j * cv):
    Dl_p = np.prod([ck - x for x in (rpv, rmv) if abs(ck - x) > 1e-14]) \
        if ck in (rpv, rmv) else (ck**2 - 2 * Mv * ck + av**2)
    P2_p = ck**2 + cv**2 if ck in (rpv, rmv) else 2 * ck
    al[ck] = Q4n(ck) / (Dl_p * P2_p) if ck in (rpv, rmv) \
        else Q4n(ck) / ((ck**2 - 2 * Mv * ck + av**2) * 2 * ck)

def integ_z(rv):
    # dphi/dz * dz/dr = [(E^2-1)+sum alpha/(r-c)] / sqrt(Q4)
    s = (Ev**2 - 1) + sum(al[ck] / (rv - ck) for ck in al)
    return (av / Ev) * s.real / np.sqrt(Q4n(rv))

I_z = quad(integ_z, r_a, r_b, limit=200)[0]
print(f"  phi (diretto)      = {I_r:.12f}")
print(f"  phi (decomposto z) = {I_z:.12f}   diff = {abs(I_r-I_z):.2e}")
print("\nFATTO.")
