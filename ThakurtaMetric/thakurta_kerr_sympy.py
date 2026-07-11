# -*- coding: utf-8 -*-
"""
Thakurta ORIGINALE (1981): Kerr conforme, g = A(eta)^2 * g_Kerr.
Brachistocrone equatoriali con rotaia conforme Ehat = -u_eta.

Base Kerr equatoriale (Boyer-Lindquist, spin s, ' = d/deta):

    ds^2 = A(eta)^2 [ -f deta^2 - (4Ms/r) deta dphi + (r^2/Dl) dr^2 + P dphi^2 ]
    f = 1-2M/r ,  Dl = r^2-2Mr+s^2 ,  P = r^2+s^2+2Ms^2/r

Vincolo -u_eta = Ehat  =>  Lam = dtau/deta = A^2 (f + (2Ms/r) phi')/Ehat.
Eliminando Lam, l'indicatrice e' una QUADRICA in (r', phi'):

    (r^2/Dl) r'^2 + Pb phi'^2 - (4Ms/r) vb^2 phi' - f vb^2 = 0

    vb^2 = 1 - A^2 f/Ehat^2          (il "v barrato": congelamento conforme)
    Pb   = P + (2Ms A/(r Ehat))^2 * A^2   [P + A^2(2Ms/r)^2/Ehat^2]

ELLISSE con centro spostato SOLO in phi' (vento angolare = frame dragging):

    phi'_0 = (2Ms/r) vb^2 / Pb ,   R^2 = f vb^2 + Pb phi'_0^2

Hamiltoniane (max PMP sull'ellisse, forma chiusa):

    H_eta = p_phi phi'_0 + R sqrt( (Dl/r^2) p_r^2 + p_phi^2/Pb ) - 1
    H_tau = pt_phi phi'_0 + R sqrt( (Dl/r^2) p_r^2 + pt_phi^2/Pb ) - A^2 f/Ehat
            con  pt_phi = p_phi - 2Ms A^2/(r Ehat)
    (il costo tau dipende dal controllo phi': shift gravitomagnetico
     CONFORME del momento angolare, ∝ A(eta)^2)

Verifiche:
  V1  quadrica e forma ellittica: identita' sympy esatte
  V2  s -> 0: si riducono alle Hamiltoniane Thakurta-Schwarzschild
  V3  A = 1 (Kerr statico): dphi_BL/dr dal flusso H=0 vs forme chiuse
      doranTau/doranT:  K r sqrt(wf) / ( Dl sqrt(Dl - K^2 w) ),
      K_tau = J,  K_t = (f J + 2Ms/r)/E   [w = E^2 - f]
  V4  superficie di luce che respira: R^2 = 0.
      Fuori dall'ergosfera (f>0): solo vb = 0 (congelamento conforme
      A^2 f = Ehat^2). Nell'ergoregione (f<0): R^2 = 0 con vb != 0 —
      la superficie di luce di Kerr, ora dipendente da A(eta).
"""

import numpy as np
import sympy as sp
from scipy.optimize import brentq

r, M, s_, E, J, pr, A = sp.symbols('r M s Ehat J p_r A', positive=True)
rp_, php_ = sp.symbols("r' phi'", real=True)

f = 1 - 2 * M / r
Dl = r**2 - 2 * M * r + s_**2
P = r**2 + s_**2 + 2 * M * s_**2 / r
w = E**2 - f

print("=" * 72)
print("[V1] quadrica dell'indicatrice e forma ellittica")
print("=" * 72)
Lam2 = A**2 * (f + (4*M*s_/r)*php_ - (r**2/Dl)*rp_**2 - P*php_**2)
Fk = f + (2*M*s_/r)*php_
quadrica = sp.expand(Lam2 - A**4 * Fk**2 / E**2)     # = 0 sul vincolo
vb2 = 1 - A**2 * f / E**2
Pb = P + A**2 * (2*M*s_/r)**2 / E**2
forma = -A**2 * ((r**2/Dl)*rp_**2 + Pb*php_**2
                 - (4*M*s_/r)*vb2*php_ - f*vb2)
print("  quadrica - forma dichiarata =", sp.simplify(quadrica - forma))
php0 = (2*M*s_/r) * vb2 / Pb
R2 = f * vb2 + Pb * php0**2
# completamento del quadrato
resto = sp.simplify(Pb*(php_ - php0)**2 + (r**2/Dl)*rp_**2 - R2
                    - ((r**2/Dl)*rp_**2 + Pb*php_**2
                       - (4*M*s_/r)*vb2*php_ - f*vb2))
print("  completamento quadrato (atteso 0):", resto)
print("  => ellisse: centro (0, phi'_0), phi'_0 = (2Ms/r) vb^2/Pb  (vento)")
# identita' chiave: f P + (2Ms/r)^2 = Dl  =>  R^2 = Dl vb^2 / Pb
idk = sp.simplify(f * P + (2*M*s_/r)**2 - Dl)
print("  identita'  f P + (2Ms/r)^2 - Dl =", idk)
print("  =>  R^2 = Dl * vb^2 / Pb  (verifica:",
      sp.simplify(R2 - Dl * vb2 / Pb), ")")
print("  degenerazioni: Dl = 0 (ORIZZONTE, conformemente invariante)")
print("                 vb = 0 (congelamento conforme A^2 f = Ehat^2)")

print()
print("=" * 72)
print("[V2] limite s -> 0: Thakurta-Schwarzschild")
print("=" * 72)
pphi = sp.Symbol('p_phi', real=True)
R = sp.sqrt(R2)
H_eta = pphi * php0 + R * sp.sqrt((Dl/r**2)*pr**2 + pphi**2/Pb) - 1
ptphi = pphi - 2*M*s_*A**2/(r*E)
H_tau = ptphi * php0 + R * sp.sqrt((Dl/r**2)*pr**2 + ptphi**2/Pb) \
    - A**2 * f / E
target = sp.sqrt(vb2) * sp.sqrt(f) * sp.sqrt(f*pr**2 + pphi**2/r**2) - 1
target_t = sp.sqrt(vb2)*sp.sqrt(f)*sp.sqrt(f*pr**2 + pphi**2/r**2) \
    - A**2*f/E
# i residui simbolici sono prodotti di radicali uguali (sympy non li
# combina senza assunzioni): verifica numerica su punti casuali
rng = np.random.default_rng(3)
ok_e, ok_t = True, True
for _ in range(8):
    sub = {M: 1.0, E: float(rng.uniform(1.05, 2.0)),
           A: float(rng.uniform(0.5, 1.2)), r: float(rng.uniform(3, 12)),
           pr: float(rng.uniform(-2, 2)), pphi: float(rng.uniform(-2, 2))}
    if float(vb2.subs(s_, 0).subs(sub)) <= 0:
        continue
    ok_e = ok_e and abs(float((H_eta.subs(s_, 0) - target).subs(sub))) < 1e-12
    ok_t = ok_t and abs(float((H_tau.subs(s_, 0) - target_t).subs(sub))) < 1e-12
print("  H_eta(s=0) == H_Thakurta (numerico):", ok_e)
print("  H_tau(s=0) == H_Thakurta (numerico):", ok_t)

print()
print("=" * 72)
print("[V3] A=1 (Kerr statico): flusso H=0 vs forme chiuse doranTau/doranT")
print("=" * 72)
M_n, s_n, E_n, J_n = 1.0, 0.9, 1.2, 1.3
subs_n = [(M, M_n), (s_, s_n), (E, E_n), (A, 1)]

def check_ramo(Hexpr, Kfun, nome):
    Hn = Hexpr.subs(subs_n).subs(pphi, J_n)
    dHdp = sp.lambdify((r, pr), sp.diff(Hn, pr), 'numpy')
    dHdJ = sp.lambdify((r, pr), sp.diff(Hexpr.subs(subs_n), pphi)
                       .subs(pphi, J_n), 'numpy')
    Hfun = sp.lambdify((r, pr), Hn, 'numpy')
    ok = True
    for rv in (4.0, 6.0, 9.0):
        pg = np.linspace(-80, 80, 320001)
        Hg = Hfun(rv, pg)
        roots = [brentq(lambda p: Hfun(rv, p), pg[i], pg[i+1])
                 for i in range(len(pg)-1)
                 if np.isfinite(Hg[i]) and np.isfinite(Hg[i+1])
                 and Hg[i]*Hg[i+1] <= 0]
        p_in = [p_ for p_ in roots if dHdp(rv, p_) < 0][0]
        num = abs(dHdJ(rv, p_in) / dHdp(rv, p_in))
        fv = 1 - 2*M_n/rv
        wv = E_n**2 - fv
        Dv = rv**2 - 2*M_n*rv + s_n**2
        Kv = Kfun(rv, fv)
        cf = Kv * rv * np.sqrt(wv*fv) / (Dv * np.sqrt(Dv - Kv**2*wv))
        ok = ok and abs(num - cf) < 1e-9 * cf
        print(f"  ramo {nome} r={rv}: dphi/dr flusso = {num:.10f}"
              f"   forma chiusa = {cf:.10f}")
    print(f"  ramo {nome}: match:", ok)

check_ramo(H_tau, lambda rv, fv: J_n, "tau")
check_ramo(H_eta, lambda rv, fv: (fv*J_n + 2*M_n*s_n/rv)/E_n, "t  ")

print()
print("=" * 72)
print("[V4] due superfici critiche:  R^2 = Dl vb^2/Pb = 0")
print("=" * 72)
R2n = sp.lambdify((r, A), R2.subs([(M, M_n), (s_, s_n), (E, E_n)]), 'numpy')
r_p = M_n + np.sqrt(M_n**2 - s_n**2)
print(f"  (a) Dl = 0: orizzonte r_+ = M+sqrt(M^2-s^2) = {r_p:.6f}")
print("      CONFORMEMENTE INVARIANTE: non dipende da A(eta) ne' da Ehat")
for Av in (1.0, 1.1, 1.3):
    rL = brentq(lambda rv: R2n(rv, Av), 1.2, 2.0 - 1e-9)
    print(f"      A = {Av:4.2f}:  radice interna di R^2 = {rL:.6f}"
          f"   (= r_+, diff {abs(rL-r_p):.1e})")
print("  (b) vb = 0: congelamento conforme A^2 f = Ehat^2 — RESPIRA:")
for Av in (1.25, 1.4, 1.6, 2.0, 3.0):
    fz = E_n**2 / Av**2
    if fz < 1:
        r_fr = 2*M_n/(1-fz)
        print(f"      A = {Av:4.2f}:  r_freeze = {r_fr:9.4f}")
    else:
        print(f"      A = {Av:4.2f}:  nessun congelamento "
              f"(Ehat^2/A^2 = {fz:.3f} > 1): tutto accessibile")
print("  => l'indicatrice sopravvive REGOLARE attraverso l'ergosfera")
print("     (f=0: R^2 = vb^2 Dl/Pb > 0) fino all'orizzonte — come il")
print("     worldline vincolato in Doran; il congelamento conforme scende")
print("     dall'infinito (A > Ehat) e comprime la regione accessibile")
print("     nella shell orizzonte-congelamento: versione rotante della")
print("     cattura da espansione.")
print("\nFATTO.")
