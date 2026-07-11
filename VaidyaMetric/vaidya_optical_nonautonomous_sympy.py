# -*- coding: utf-8 -*-
"""
"Metriche ottiche" non autonome per Vaidya generica — senso debole.

Invertendo l'indicatrice (Legendre-inversa delle Hamiltoniane §3d),
per OGNI m(v) si ottengono forme ottiche chiuse, con f = f(v,r):

    ramo tau:  dtau = sqrt( dr^2 + f r^2 dphi^2 ) / sqrt(w)
               [Riemann-Jacobi tempo-dipendente]
    ramo v:    dv   = dr/f + (E/(f sqrt(w))) sqrt( dr^2 + f r^2 dphi^2 )
               [Randers tempo-dipendente: deriva dr/f + norma]

    w = E^2 - f,  f = 1 - 2m(v)/r.

Con m' = 0 sono ESATTAMENTE le due metriche di Perlick 1991 in EF
(doranTau: L_tau = sqrt((B^2+Af)/w) con B=1, A = r^2 phi'^2, g_rr=0):
Perlick generalizzato verbatim con f -> f(v,r).

MA: senza CKV nulla e' rigido — deriva E norma dipendono entrambe da v.
Non esiste geometria statica sottostante (tranne m = mu*v, omotetia,
R5-R6): "ottico" solo in senso non autonomo. La gerarchia:

    FLRW           indice separabile -> orologio, geometria banale
    Thakurta-Schw  geometria statica x indice n(eta,r)
    Thakurta-Kerr  geometria statica + vento RIGIDO x indice
    Vaidya mu*v    geometria statica sul quoziente autosimilare
    Vaidya m(v)    tutto respira: solo senso debole (= duale delle H)

Verifiche:
  V1  entrambe le forme soddisfano il vincolo/indicatrice per m(v)
      qualsiasi (residui numerici ~1e-14 su punti casuali)
  V2  m'=0: dtau riproduce L_tau di doranTau (B=1, A=r^2 phi'^2, a=0)
  V3  EL della forma dv (m'=0, r come parametro): dphi/dr = forma
      chiusa doranT a=0 (K = fJ/E)
"""

import numpy as np
import sympy as sp
from scipy.optimize import brentq

v, r, E, J = sp.symbols('v r E J', positive=True)
dr_, dphi_, dv_ = sp.symbols('dr dphi dv', real=True)
m = sp.Function('m', positive=True)

f = 1 - 2 * m(v) / r
w = E**2 - f

dl2 = dr_**2 + f * r**2 * dphi_**2
dtau_opt = sp.sqrt(dl2) / sp.sqrt(w)
dv_opt = dr_ / f + (E / (f * sp.sqrt(w))) * sp.sqrt(dl2)

print("=" * 72)
print("[V1] le forme ottiche soddisfano il vincolo (m(v) qualsiasi)")
print("=" * 72)
# worldline: -ds^2 = f dv^2 - 2 dv dr - r^2 dphi^2 = dtau^2  e
# vincolo -u_v = E  <=>  dtau = (f dv - dr)/E
# => eliminando: f w dv^2 - 2 w dr dv - (dr^2 + E^2 r^2 dphi^2) = 0
quadrica = f * w * dv_**2 - 2 * w * dr_ * dv_ \
    - (dr_**2 + E**2 * r**2 * dphi_**2)
res_v = quadrica.subs(dv_, dv_opt)
res_tau = (f * dv_opt - dr_) / E - dtau_opt      # coerenza dtau = (f dv-dr)/E

rng = np.random.default_rng(11)
w1 = w2 = 0.0
mv = sp.Symbol('mv', positive=True)
for _ in range(10):
    sub = {E: float(rng.uniform(1.05, 2.5)), r: float(rng.uniform(2.5, 12)),
           m(v): float(rng.uniform(0.5, 1.4)),
           dr_: float(rng.uniform(-1, 1)), dphi_: float(rng.uniform(-1, 1))}
    w1 = max(w1, abs(float(res_v.subs(sub))))
    w2 = max(w2, abs(float(res_tau.subs(sub))))
print(f"  quadrica(dv = Randers):        max residuo = {w1:.2e}")
print(f"  dtau_opt = (f dv_opt - dr)/E:  max residuo = {w2:.2e}")
print("  => valgono per OGNI m(v): f = f(v,r) dentro le forme")

print()
print("=" * 72)
print("[V2] m'=0: dtau = L_tau di doranTau (B=1, A=r^2 phi'^2, a=0)")
print("=" * 72)
php = sp.Symbol("php", real=True)
L_tau_doran = sp.sqrt((1 + (r**2 * php**2) * f) / w)     # (B^2+Af)/w, B=1
L_tau_nostra = dtau_opt.subs([(dr_, 1), (dphi_, php)])   # per unita' di dr
diff = L_tau_doran - L_tau_nostra
w3 = 0.0
for _ in range(8):
    sub = {E: float(rng.uniform(1.05, 2.5)), r: float(rng.uniform(2.5, 12)),
           m(v): float(rng.uniform(0.5, 1.4)),
           php: float(rng.uniform(-1, 1))}
    w3 = max(w3, abs(float(diff.subs(sub))))
print(f"  L_tau(doranTau) - dtau_opt/dr: max residuo = {w3:.2e}")

print()
print("=" * 72)
print("[V3] m'=0: EL di dv_opt vs forma chiusa doranT (a=0)")
print("=" * 72)
M_n, E_n, J_n = 1.0, 1.2, 1.3
L_v = dv_opt.subs([(dr_, 1), (dphi_, php), (m(v), M_n), (E, E_n)])
pJ = sp.diff(L_v, php)
ok = True
for rv in (4.0, 6.0, 9.0):
    pJ_r = sp.lambdify(php, pJ.subs(r, rv), 'numpy')
    php_star = brentq(lambda x: pJ_r(x) - J_n, -3, 3)
    fv = 1 - 2 * M_n / rv
    wv = E_n**2 - fv
    Dv = rv**2 * fv
    Kv = fv * J_n / E_n
    cf = Kv * rv * np.sqrt(wv * fv) / (Dv * np.sqrt(Dv - Kv**2 * wv))
    ok = ok and abs(php_star - cf) < 1e-9 * cf
    print(f"  r={rv}:  dphi/dr da dL/dphi'=J: {php_star:.10f}"
          f"   forma chiusa: {cf:.10f}")
print("  match:", ok)
print("\n  => Perlick 1991 generalizzato verbatim con f -> f(v,r);")
print("     senza CKV pero' deriva E norma respirano entrambe:")
print("     'ottico' solo in senso non autonomo (duale delle Hamiltoniane).")
print("\nFATTO.")
