# -*- coding: utf-8 -*-
"""
Metrica ottica dipendente dal tempo conforme per Thakurta-Kerr.

Dall'indicatrice (equatoriale, rotaia Ehat = -u_eta):

    (r^2/Dl) dr^2 + Pb (dphi - phi'_0 deta)^2 = R^2 deta^2 ,
    R^2 = Dl vb^2/Pb ,  phi'_0 = (2Ms/r) vb^2/Pb ,  vb^2 = 1 - A^2 f/Ehat^2

risolvendo la quadratica in deta (radice futura):

    deta = n(eta, r) * alpha_K  +  beta_K                    [RANDERS]

    alpha_K^2 = r^2 dr^2/(f Dl) + Dl dphi^2/f^2      (parte Riemanniana)
    beta_K    = - (2Ms/(r f)) dphi                   (1-forma gravitomagnetica)
    n(eta,r)  = 1/vb = Ehat / sqrt(Ehat^2 - A(eta)^2 f)

TEOREMA (verificato sotto): alpha_K e beta_K sono ESATTAMENTE i dati di
Randers ottici della LUCE in Kerr (Fermat nullo, Gibbons-Werner):
tutta la fisica di massa/rotaia/espansione sta nell'indice scalare
n(eta,r), che moltiplica SOLO la parte Riemanniana. La 1-forma
(trascinamento) e' rigida: invariante conforme e indipendente da Ehat.

Limiti:
    Ehat -> oo : n -> 1: metrica ottica nulla di Kerr (inv. conforme ok)
    A = 1      : n = E/sqrt(w): Perlick stazionario su Kerr
    s = 0      : beta = 0, alpha = ottica di Schwarzschild: Thakurta-Schw
    f = 1(s=0,M=0): n = Ehat/sqrt(Ehat^2-A^2) = 1/v FLRW

Singolarita':
    f = 0 (ergosfera): alpha e beta divergono — la RIDUZIONE ottica
        muore all'ergosfera (come in Kerr statico: 1/F intrinseco);
        l'indicatrice/Hamiltoniana resta regolare fino all'orizzonte.
    vb = 0 (congelamento conforme): n -> oo, "istante-superficie" mobile.

Verifiche:
  V1  deta = n alpha + beta soddisfa la quadrica dell'indicatrice
      (identita' sympy esatta, radice futura)
  V2  beta e alpha coincidono coi dati di Randers del Fermat NULLO di
      Kerr (derivati indipendentemente da ds^2 = 0)
  V3  limiti: s=0, A=1, Ehat->oo, e f->1 (FLRW)
  V4  numerico: dphi/dr dall'EL di F = n alpha + beta (a eta fisso piu'
      termine non autonomo nullo per A=1) vs forma chiusa doranT
"""

import numpy as np
import sympy as sp
from scipy.optimize import brentq

r, M, s_, E, A = sp.symbols('r M s Ehat A', positive=True)
dr_, dphi_, deta_ = sp.symbols('dr dphi deta', real=True)

f = 1 - 2 * M / r
Dl = r**2 - 2 * M * r + s_**2
P = r**2 + s_**2 + 2 * M * s_**2 / r
vb2 = 1 - A**2 * f / E**2
Pb = P + A**2 * (2 * M * s_ / r)**2 / E**2
php0 = (2 * M * s_ / r) * vb2 / Pb
R2 = Dl * vb2 / Pb

alpha2 = r**2 * dr_**2 / (f * Dl) + Dl * dphi_**2 / f**2
beta = -(2 * M * s_ / (r * f)) * dphi_
n_idx = 1 / sp.sqrt(vb2)
F = n_idx * sp.sqrt(alpha2) + beta

print("=" * 72)
print("[V1] F = n*alpha + beta soddisfa l'indicatrice")
print("=" * 72)
quadrica = (r**2 / Dl) * dr_**2 + Pb * (dphi_ - php0 * deta_)**2 \
    - R2 * deta_**2
res_expr = quadrica.subs(deta_, F)
# radicali annidati: simplify esplode — verifica numerica esatta
rng = np.random.default_rng(7)
worst = 0.0
for _ in range(10):
    sub = {M: 1.0, s_: float(rng.uniform(0.1, 0.95)),
           E: float(rng.uniform(1.05, 2.5)), A: float(rng.uniform(0.4, 1.3)),
           r: float(rng.uniform(2.5, 12.0)),
           dr_: float(rng.uniform(-1, 1)), dphi_: float(rng.uniform(-1, 1))}
    if float(vb2.subs(sub)) <= 0:
        continue
    worst = max(worst, abs(float(res_expr.subs(sub))))
print(f"  quadrica(deta = n*alpha + beta): max residuo su 10 punti"
      f" casuali = {worst:.2e}")

print()
print("=" * 72)
print("[V2] alpha, beta = dati di Randers del Fermat NULLO di Kerr")
print("=" * 72)
# ds^2 = 0 equatoriale Kerr (base, A qualsiasi si cancella):
# -f dt^2 - (4Ms/r) dt dphi + (r^2/Dl) dr^2 + P dphi^2 = 0
dt = sp.Symbol('dt', real=True)
null_cond = -f * dt**2 - (4 * M * s_ / r) * dt * dphi_ \
    + (r**2 / Dl) * dr_**2 + P * dphi_**2
# radice futura: dt = -(2Ms/(rf))dphi + sqrt(alpha_K^2)
F_null = beta + sp.sqrt(alpha2)
res2_expr = null_cond.subs(dt, F_null)
worst2 = 0.0
for _ in range(10):
    sub = {M: 1.0, s_: float(rng.uniform(0.1, 0.95)),
           r: float(rng.uniform(2.5, 12.0)),
           dr_: float(rng.uniform(-1, 1)), dphi_: float(rng.uniform(-1, 1))}
    worst2 = max(worst2, abs(float(res2_expr.subs(sub))))
print(f"  condizione nulla( dt = beta + alpha ): max residuo = {worst2:.2e}")
print("  => luce: n=1; massiva+espansione: SOLO n cambia, la 1-forma")
print("     gravitomagnetica beta = -(2Ms/rf) dphi e' RIGIDA")
print("     (invariante conforme e indipendente da Ehat).")

print()
print("=" * 72)
print("[V3] limiti")
print("=" * 72)
F_s0 = F.subs(s_, 0)
target_s0 = (1 / sp.sqrt(vb2.subs(s_, 0))) \
    * sp.sqrt(dr_**2 / f**2 + r**2 * dphi_**2 / f)
worst3 = 0.0
for _ in range(8):
    sub = {M: 1.0, E: float(rng.uniform(1.05, 2.5)),
           A: float(rng.uniform(0.4, 1.3)), r: float(rng.uniform(2.5, 12.0)),
           dr_: float(rng.uniform(-1, 1)), dphi_: float(rng.uniform(-1, 1))}
    if float(vb2.subs(s_, 0).subs(sub)) <= 0:
        continue
    worst3 = max(worst3, abs(float((F_s0 - target_s0).subs(sub))))
print(f"  s=0:  F - n*alpha_Schw: max residuo = {worst3:.2e}")
n_A1 = sp.simplify(n_idx.subs(A, 1))
print("  A=1:  n =", n_A1, " = E/sqrt(w)  (Perlick su Kerr)")
n_null = sp.limit(n_idx, E, sp.oo)
print("  Ehat->oo:  n =", n_null, " (ottica nulla di Kerr)")
n_flrw = sp.simplify(n_idx.subs([(M, 0), (s_, 0)]))
print("  M=s=0:  n =", n_flrw, " = 1/v FLRW  (orologio lam)")

print()
print("=" * 72)
print("[V4] numerico: EL di F (A=1) vs forma chiusa doranT")
print("=" * 72)
M_n, s_n, E_n, J_n, A_n = 1.0, 0.9, 1.2, 1.3, 1.0
# EL con r come parametro: F_r(phi', r) = n*sqrt(a_rr + a_pp phi'^2) + b phi'
php = sp.Symbol("php", real=True)
F_r = (n_idx * sp.sqrt(r**2 / (f * Dl) + Dl * php**2 / f**2)
       - (2 * M * s_ / (r * f)) * php)
F_rn = F_r.subs([(M, M_n), (s_, s_n), (E, E_n), (A, A_n)])
# momento di Fermat conservato: dF/dphi' = J
pJ = sp.diff(F_rn, php)
ok = True
for rv in (4.0, 6.0, 9.0):
    pJ_r = sp.lambdify(php, pJ.subs(r, rv), 'numpy')
    lo, hi = -3.0, 3.0
    php_star = brentq(lambda x: pJ_r(x) - J_n, lo, hi)
    fv = 1 - 2 * M_n / rv
    wv = E_n**2 - fv
    Dv = rv**2 - 2 * M_n * rv + s_n**2
    Kv = (fv * J_n + 2 * M_n * s_n / rv) / E_n
    cf = Kv * rv * np.sqrt(wv * fv) / (Dv * np.sqrt(Dv - Kv**2 * wv))
    ok = ok and abs(php_star - cf) < 1e-9 * cf
    print(f"  r={rv}:  dphi/dr da dF/dphi'=J: {php_star:.10f}"
          f"   forma chiusa doranT: {cf:.10f}")
print("  match (il ramo eta della metrica ottica = ramo t di Kerr):", ok)
print("\nFATTO.")
