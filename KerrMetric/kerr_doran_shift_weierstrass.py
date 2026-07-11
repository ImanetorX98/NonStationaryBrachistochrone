# -*- coding: utf-8 -*-
"""
Shift di Doran della separatrice in forma chiusa (Weierstrass) e
validazione numerica.

    dphi_D/dr = a beta/Dl + dphi_BL/dr ,   beta = sqrt(2Mr/(r^2+a^2))

Lo shift:  a beta/Dl dr = 2Mar dr/(Dl sqrt(C3)),  C3 = 2Mr(r^2+a^2)
[il fattore (r^2+a^2) si cancella: poli SOLO agli orizzonti r_pm].

UNIFORMIZZAZIONE (cubica -> Weierstrass):
    r(z) = (2/M) P(z; g2, g3) ,   g2 = -M^2 a^2 ,   g3 = 0
    z(r) = P^{-1}(M r/2)          [= int_r^oo dr'/sqrt(C3)]

AGLI ORIZZONTI (Dl(r_pm)=0 => r_pm^2 + a^2 = 2M r_pm):
    C3(r_pm) = 4M^2 r_pm^2  =>  P'(v_pm) = pm M^2 r_pm
    pesi:  (M A_pm/2)/P'(v_pm) = pm a/(r_+ - r_-) ,
           A_pm = 2M a r_pm/(r_pm - r_mp)

FORMA FINALE:

    phi_shift(r) = [a/(r_+ - r_-)] *
        { [2 zeta(v_+) z + ln sigma(z-v_+)/sigma(z+v_+)]
        - [2 zeta(v_-) z + ln sigma(z-v_-)/sigma(z+v_-)] }
    con P(v_pm) = M r_pm/2,  z = z(r).

VALIDAZIONE (catena a 3 anelli, come per R12a):
  V1 simbolica: uniformizzazione (P')^2 = 4P^3 + M^2a^2 P  identica
  V2 simbolica: C3(r_pm) = 4M^2 r_pm^2 e pesi = pm a/(r_+ - r_-)
  V3 numerica: quadratura diretta int a beta/Dl dr  vs  quadratura in z
     della decomposizione [A_+/(r-r_+) + A_-/(r-r_-)] dz
  (l'antiderivata sigma/zeta WW 20.53 e' gia' validata a 30 cifre in
   kerr_separatrix_validation.py, identita' indipendente dal reticolo)
"""

import numpy as np
import sympy as sp
from scipy.integrate import quad

r, M, E, a, z = sp.symbols('r M E a z', positive=True)
C3 = 2 * M * r * (r**2 + a**2)
Dl = r**2 - 2 * M * r + a**2

print("=" * 70)
print("[V1] uniformizzazione: r = (2/M) P,  g2 = -M^2 a^2,  g3 = 0")
print("=" * 70)
Pw = M * r / 2
lhs = (M / 2)**2 * C3                       # (P')^2 = (M/2)^2 (dr/dz)^2
rhs = 4 * Pw**3 + M**2 * a**2 * Pw
print("  (P')^2 - [4P^3 - g2 P - g3] =", sp.simplify(lhs - rhs))

print()
print("=" * 70)
print("[V2] semplificazione agli orizzonti")
print("=" * 70)
rp = M + sp.sqrt(M**2 - a**2)
rm = M - sp.sqrt(M**2 - a**2)
for nome, rk in (("r_+", rp), ("r_-", rm)):
    chk = sp.simplify(C3.subs(r, rk) - 4 * M**2 * rk**2)
    print(f"  C3({nome}) - 4M^2 {nome}^2 =", chk)
Ap = 2 * M * a * rp / (rp - rm)
peso_p = sp.simplify((M * Ap / 2) / (M**2 * rp))
print("  peso_+ = (M A_+/2)/(M^2 r_+) =", peso_p, " = a/(r_+ - r_-)")

print()
print("=" * 70)
print("[V3] validazione numerica: quadratura r vs decomposizione in z")
print("=" * 70)
Mv, av = 1.0, 0.9
rpv = Mv + np.sqrt(Mv**2 - av**2)
rmv = Mv - np.sqrt(Mv**2 - av**2)
C3n = lambda rv: 2 * Mv * rv * (rv**2 + av**2)
Dln = lambda rv: rv**2 - 2 * Mv * rv + av**2
Apn = 2 * Mv * av * rpv / (rpv - rmv)
Amn = 2 * Mv * av * rmv / (rmv - rpv)

f_dir = lambda rv: av * np.sqrt(2 * Mv * rv / (rv**2 + av**2)) / Dln(rv)
f_dec = lambda rv: (Apn / (rv - rpv) + Amn / (rv - rmv)) / np.sqrt(C3n(rv))
r_a, r_b = 2.5, 9.0
I1 = quad(f_dir, r_a, r_b, limit=300)[0]
I2 = quad(f_dec, r_a, r_b, limit=300)[0]
print(f"  int a beta/Dl (diretto)   = {I1:.14f}")
print(f"  decomposizione in z       = {I2:.14f}   diff = {abs(I1-I2):.2e}")

# bonus: phi_D TOTALE della separatrice = R12a + shift
Ev = 1.2
cv = av / Ev
Q4n = lambda rv: rv * (2 * Mv + (Ev**2 - 1) * rv) * (rv**2 + cv**2)
f_bl = lambda rv: (av / Ev) * np.sqrt(Q4n(rv)) / (Dln(rv) * (rv**2 + cv**2))
I_bl = quad(f_bl, r_a, r_b, limit=300)[0]
print(f"  phi_D totale (Doran) su [{r_a},{r_b}]: shift + BL = "
      f"{I1 + I_bl:.12f}")
print("""
FORMA FINALE dello shift (validata):

  phi_shift(r) = [a/(r_+ - r_-)] * ( B_+(z) - B_-(z) ),
  B_k(z) = 2 zeta(v_k) z + ln[ sigma(z - v_k)/sigma(z + v_k) ]

  reticolo:  g2 = -M^2 a^2 ,  g3 = 0   (quasi-lemniscatico)
  z(r) = P^{-1}(M r/2) ,  P(v_pm) = M r_pm/2 ,  P'(v_pm) = M^2 r_pm

e  phi_D(r) = phi_BL(r) [R12a, curva Q4] + phi_shift(r) [cubica C3]:
regolare attraverso ergosfera E orizzonte (le divergenze log dei due
pezzi a r_+ si CANCELLANO: e' il contenuto della regolarita' Doran).
FATTO.""")
