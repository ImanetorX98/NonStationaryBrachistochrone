# -*- coding: utf-8 -*-
"""
Brachistocrone t e tau in Vaidya entrante — worldline vincolato NON autonomo.

Metrica (equatoriale, c=G=1):
    ds^2 = -f dv^2 + 2 dv dr + r^2 dphi^2 ,   f = 1 - 2 m(v)/r

Vincolo rotaia (generalizzazione diretta di Kerr, -u_t = E):
    -u_v = E  fissata          (d_v NON e' Killing: la rotaia la impone,
                                 con forza f_v != 0 come f_eta = aH in FLRW)

Riduzione (parametro lambda = r, tratto con r crescente; ' = d/dr):
    Lam = dtau/dr = sqrt(f v'^2 - 2 v' - r^2 phi'^2)
    vincolo => Lam = (f v' - 1)/E
    eliminando Lam:   phi'^2 = W(v, v', r)
        W = [ f v'^2 - 2 v' - (f v' - 1)^2/E^2 ] / r^2

Funzionali ridotti (J = moltiplicatore di Fermat per Delta_phi, costante
perche' phi resta ciclica — la sfericita' sopravvive):

    ramo tau:  F_tau = (f v' - 1)/E - J*sqrt(W)     [ T_tau = int Lam dr ]
    ramo t:    F_t   =  v'          - J*sqrt(W)     [ T_v  = int v'  dr ]

Struttura non autonoma:
    F dipende da v SOLO attraverso m(v)  =>  dF/dv = (dF/dm) m'(v):
    forzante proporzionale a m'. In Schwarzschild (m'=0) p_v = dF/dv'
    e' conservato; l'arrivo spaziale con v libero da' la trasversalita'
    p_v(r_arrivo) = 0  =>  caso stazionario: p_v == 0 lungo tutta l'orbita.

Verifiche:
  V1  dF/dm stampato: il forzante non autonomo esplicito (chain rule esatta)
  V2  stazionario (m'=0, p_v=0): dphi/dr in forma chiusa =
        ramo tau:  J  sqrt(w f) r / ( Dl sqrt(Dl - J^2 w) )
        ramo t:    K  sqrt(w f) r / ( Dl sqrt(Dl - K^2 w) ),  K = f J / E
      con Dl = r^2 f, w = E^2 - f: le formule Kerr a=0 (doranTau.md,
      doranT.md con 2Ma/r -> 0). Check numerico multi-punto.
  V3  numerico: Vaidya lineare m(v) = m0 + mu*v, ramo tau, integrazione
      all'indietro dal punto d'arrivo con p_v(r1) = 0; confronto con
      Schwarzschild (mu=0) a parita' di (E, J, r0, r1, v_arrivo).

Tutto in simboli piani (niente Function(r)): l'equazione EL usa la
derivata totale  d/dr = d_r + v' d_v + v'' d_v' esplicita.
"""

import sympy as sp
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

r, vv, uu, q = sp.symbols('r vv uu q', real=True)     # r, v, v', v''
E, J, mm, mp, mu = sp.symbols('E J mm mp mu', positive=True)  # m(v), m'(v)

f = 1 - 2 * mm / r
W = (f * uu**2 - 2 * uu - (f * uu - 1)**2 / E**2) / r**2
phip = sp.sqrt(W)

F_tau = (f * uu - 1) / E - J * phip
F_t = uu - J * phip

def eq_moto(F):
    """EL: d/dr(dF/du) = dF/dv con m=m(v): ritorna v'' = RHS e p_v."""
    p = sp.diff(F, uu)
    # derivata totale lungo la curva: m -> m(v(r)), dm/dr = mp*uu
    dp_tot = (sp.diff(p, r) + sp.diff(p, vv) * uu + sp.diff(p, mm) * mp * uu
              + sp.diff(p, uu) * q)
    dFdv = sp.diff(F, mm) * mp                      # chain rule: F(v) via m(v)
    EL = dp_tot - dFdv
    vpp = sp.solve(sp.Eq(EL, 0), q)[0]              # lineare in q: veloce
    return vpp, p

print("=" * 72)
print("[V1] forzante non autonomo: dF/dv = (dF/dm) * m'(v)  (chain rule)")
print("=" * 72)
for nome, F in (("tau", F_tau), ("t  ", F_t)):
    dFdm = sp.radsimp(sp.diff(F, mm))
    print(f"\n  ramo {nome}:  dF/dm =", sp.factor_terms(dFdm))
print("""
  => m' = 0: p_v conservato; trasversalita' p_v(r_arrivo)=0 => p_v == 0.
     m' != 0: dp_v/dr ~ m'(v) != 0; p_v(r)!=0 lungo il percorso pur con
     p_v(r1)=0: la brachistocrona 'anticipa' l'accrescimento futuro.""")

print("=" * 72)
print("[V2] limite stazionario: p_v = 0  =>  forma chiusa Kerr a=0")
print("=" * 72)

w = E**2 - f
Dl = r**2 * f

def vp_radiale(f_val, E_val):
    """Radice grande di W=0 in v' (moto radiale, phi'=0)."""
    a2 = f_val - f_val**2 / E_val**2
    a1 = -2 + 2 * f_val / E_val**2
    a0 = -1 / E_val**2
    return (-a1 + np.sqrt(a1**2 - 4 * a2 * a0)) / (2 * a2)

vals = {mm: 1.0, E: 1.2, J: 1.3}
for nome, F, K in (("tau", F_tau, J), ("t  ", F_t, f * J / E)):
    pv = sp.diff(F, uu)
    target = K * sp.sqrt(w * f) * r / (Dl * sp.sqrt(Dl - K**2 * w))
    ok = True
    for rv in (4.0, 7.3, 12.0):
        sub = dict(vals);  sub[r] = rv
        pv_n = sp.lambdify(uu, pv.subs(sub), 'numpy')
        lo = vp_radiale(float(f.subs(sub)), vals[E]) * (1 + 1e-9)
        hi = lo * 2
        while pv_n(hi) < 0:
            hi *= 2
        up_star = brentq(pv_n, lo, hi)
        phi_num = float(phip.subs(sub).subs(uu, up_star))
        phi_cf = float(target.subs(sub))
        ok = ok and abs(phi_num - phi_cf) < 1e-9 * max(1, abs(phi_cf))
        print(f"  ramo {nome} r={rv:5.1f}:  dphi/dr num = {phi_num:.12f}"
              f"   forma chiusa = {phi_cf:.12f}")
    print(f"  ramo {nome}: MATCH forma chiusa Kerr a=0:", ok)

print()
print("=" * 72)
print("[V3] numerico: ramo tau in Vaidya lineare m(v) = m0 + mu*v")
print("=" * 72)

m0v, E_n, J_n = 1.0, 1.2, 1.3
r1, r0 = 10.0, 3.5   # arrivo a r1, partenza a r0 (periasse a r~3.01: sopra)

vpp_tau, pv_tau = eq_moto(F_tau)
lin = [(mm, m0v + mu * vv), (mp, mu), (E, E_n), (J, J_n)]
args = (r, vv, uu, mu)
vpp_fn = sp.lambdify(args, vpp_tau.subs(lin), 'numpy')
pv_fn = sp.lambdify(args, pv_tau.subs(lin), 'numpy')
W_fn = sp.lambdify(args, W.subs(lin), 'numpy')
Lam_fn = sp.lambdify(args, ((f * uu - 1) / E).subs(lin), 'numpy')

def integra(mu_n, v_arrivo):
    # trasversalita': p_v(r1) = 0 fissa v'(r1)
    fl = 1 - 2 * (m0v + mu_n * v_arrivo) / r1
    lo = vp_radiale(fl, E_n) * (1 + 1e-9)
    hi = lo * 2
    while pv_fn(r1, v_arrivo, hi, mu_n) < 0:
        hi *= 2
    up1 = brentq(lambda z: pv_fn(r1, v_arrivo, z, mu_n), lo, hi)

    def rhs(rr, y):
        v_, u_ = y[0], y[1]
        return [u_, vpp_fn(rr, v_, u_, mu_n),
                np.sqrt(max(W_fn(rr, v_, u_, mu_n), 0.0)),
                Lam_fn(rr, v_, u_, mu_n)]

    sol = solve_ivp(rhs, [r1, r0], [v_arrivo, up1, 0.0, 0.0],
                    rtol=1e-10, atol=1e-12)
    y0 = sol.y[:, -1]
    return dict(v_part=y0[0], dphi=-y0[2], T_tau=-y0[3],
                T_v=v_arrivo - y0[0],
                pv_part=pv_fn(r0, y0[0], y0[1], mu_n))

v_arr = 40.0
res_v = integra(0.01, v_arr)       # accrescimento mu = 0.01
res_s = integra(0.0, v_arr)        # Schwarzschild m = m0

print(f"\n  salita r0={r0} -> r1={r1}, arrivo a v={v_arr}, E={E_n}, J={J_n}")
print(f"  {'':14s}{'Schwarzschild':>16s}{'Vaidya mu=0.01':>16s}")
for k, lab in (("v_part", "v partenza"), ("dphi", "Delta_phi"),
               ("T_v", "T_v = Dv"), ("T_tau", "T_tau"),
               ("pv_part", "p_v partenza")):
    print(f"  {lab:14s}{res_s[k]:16.6f}{res_v[k]:16.6f}")
print("\n  Vaidya: massa al passaggio m(v) in [%.3f, %.3f]" %
      (m0v + 0.01 * res_v['v_part'], m0v + 0.01 * v_arr))
print("  p_v(partenza) != 0 nel caso Vaidya (=0 esatto in Schwarzschild):")
print("  'memoria' non autonoma dell'accrescimento lungo il percorso.")
print("\nFATTO.")
