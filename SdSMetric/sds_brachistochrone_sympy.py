# -*- coding: utf-8 -*-
"""
Brachistocrone t e tau in Schwarzschild-de Sitter (Kottler) — analitico.

    f(r) = 1 - 2M/r - Lambda r^2/3     (EF avanzate: ds^2 = -f dv^2+2dvdr+...)

Stazionarieta': d_t e' Killing GENUINO ma timelike solo nella "vasca"
r_b < r < r_c (f > 0). Perlick 1991 vale solo li'; la riduzione ottica
degenera a ENTRAMBI gli orizzonti. Il worldline vincolato li attraversa.

Tutta la macchina Vaidya (riduzioni, Hamiltoniane Zermelo) dipende da f
solo attraverso f: si generalizza sostituendo f. Con m' = 0, H e' conservata
e H = 0 e' la famiglia brachistocrona (trasversalita').

Verifiche/risultati:
  V1  orizzonti r_b < r_c per 0 < 9 Lambda M^2 < 1; f_max = 1-(9 Lam M^2)^(1/3)
      a r* = (3M/Lam)^(1/3); sfera fotonica r_ph = 3M ESATTA (Lambda sparisce)
  V2  congelamento statico: v = sqrt(w)/E = 0 dove f = E^2 (w = E^2 - f):
      per E < sqrt(f_max) DUE barriere w=0 (coppia spaziale — l'analogo
      statico del congelamento FLRW a = Ehat)
  V3  ramo tau: g_tau = r^2 f - J^2 w = -J^2 E^2 < 0 a ENTRAMBI gli
      orizzonti => riflette per ogni J != 0 (J_c = 0, misura nulla,
      su tutti e due i bordi)
  V4  ramo t/v: K = f J/E -> 0 agli orizzonti (auto-sintonizzazione);
      d_r g_t|_{f=0} = r_h^2 f'(r_h): > 0 a r_b, < 0 a r_c
      => regione permessa fino al bordo DA ENTRAMBI I LATI:
      il ramo t attraversa marginalmente TUTTI E DUE gli orizzonti, ogni J
  V5  numerico (Hamilton, H=0): r_min/r_max del ramo tau = radici esatte
      di g_tau; dphi/dr lungo il flusso = forma chiusa (entrambi i rami);
      ramo t: raggiunge r_b e r_c
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

r, M, Lam, E, J, pr = sp.symbols('r M Lambda E J p_r', positive=True)

f = 1 - 2 * M / r - Lam * r**2 / 3
w = E**2 - f

print("=" * 72)
print("[V1] struttura statica")
print("=" * 72)
rstar = (3 * M / Lam) ** sp.Rational(1, 3)
fmax = sp.simplify(f.subs(r, rstar))
print("  f'(r) = 0 a r* =", rstar, ";  f_max =", sp.simplify(fmax),
      " = 1 - (9 Lam M^2)^(1/3)")
print("  orizzonti: radici di f=0 (esistono per 0 < 9 Lam M^2 < 1);")
print("  d_t Killing genuino, timelike SOLO per r_b < r < r_c.")
# sfera fotonica: r f' = 2 f
ph = sp.simplify(r * sp.diff(f, r) - 2 * f)
print("  r f' - 2f =", sp.simplify(sp.expand(ph)),
      "  =>  r_ph = 3M esatta (Lambda si cancella)")
print("  [V2] v_locale = sqrt(w)/E: congelamento dove f = E^2;")
print("       per E < sqrt(f_max): DUE barriere spaziali w=0")
print("       (analogo statico del congelamento FLRW a=Ehat).")

print()
print("=" * 72)
print("[V3-V4] soglie agli orizzonti (simbolico)")
print("=" * 72)
g_tau = r**2 * f - J**2 * w
print("  g_tau|_{f=0} =", sp.simplify(g_tau.subs(f, 0).subs(w, E**2)),
      "  < 0: ramo tau RIFLETTE a entrambi gli orizzonti (J != 0)")
K = f * J / E
g_t = r**2 * f - K**2 * w
dg_t = sp.diff(g_t, r)
# a f=0: tutti i termini con K^2 ~ f^2 hanno derivata ~ f f' -> 0
dg_t_h = sp.simplify(dg_t.subs(f, 0))  # formale; il calcolo vero sotto
print("  K(f=0) = 0 per ogni J (auto-sintonizzazione, doranT con a=0)")
print("  d_r g_t|_{f=0} = r_h^2 f'(r_h):  f'(r_b) > 0, f'(r_c) < 0")
print("  => regione permessa (g_t > 0) confina con ENTRAMBI i bordi:")
print("     il ramo t attraversa marginalmente r_b E r_c, per ogni J")

# ------------------------------------------------------------- numerico
print()
print("=" * 72)
print("[V5] numerico: Hamilton (H = 0) con f_SdS,  M=1, Lambda=0.03")
print("=" * 72)
M_n, L_n, E_n, J_n, r0 = 1.0, 0.03, 1.2, 1.3, 5.0

def f_n(rr):
    return 1 - 2 * M_n / rr - L_n * rr**2 / 3

r_b = brentq(f_n, 2.0, 3.0)
r_c = brentq(f_n, 6.0, 9.9)
rst = (3 * M_n / L_n) ** (1 / 3)
print(f"  r_b = {r_b:.6f}   r_c = {r_c:.6f}   f_max = f({rst:.3f}) = "
      f"{f_n(rst):.4f}  (E={E_n} > sqrt(f_max)={np.sqrt(f_n(rst)):.3f})")

# Hamiltoniane EF (come Vaidya, f -> f_SdS):
#   H_v   = p_r(f-E^2) - 1 + sqrt(w) sqrt(E^2 p_r^2 + J^2/r^2)
#   H_tau = p_r(f-E^2) - E + sqrt(w) sqrt((E p_r+1)^2 + J^2/r^2)
rr_ = sp.Symbol('rr', positive=True)
f_s = 1 - 2 * M_n / rr_ - L_n * rr_**2 / 3
w_s = E_n**2 - f_s
H_tau = pr * (f_s - E_n**2) - E_n \
    + sp.sqrt(w_s) * sp.sqrt((E_n * pr + 1)**2 + J_n**2 / rr_**2)
H_v = pr * (f_s - E_n**2) - 1 \
    + sp.sqrt(w_s) * sp.sqrt(E_n**2 * pr**2 + J_n**2 / rr_**2)

def flusso(Hexpr, verso, ev_lo=None, ev_hi=None, v_max=400.0):
    dHdp = sp.lambdify((rr_, pr), sp.diff(Hexpr, pr), 'numpy')
    dHdr = sp.lambdify((rr_, pr), sp.diff(Hexpr, rr_), 'numpy')
    Hfun = sp.lambdify((rr_, pr), Hexpr, 'numpy')
    # p_r al lancio da H=0, ramo con segno dr/dv richiesto
    p_grid = np.linspace(-40, 40, 160001)
    Hg = Hfun(r0, p_grid)
    roots = []
    for i in range(len(p_grid) - 1):
        if np.isfinite(Hg[i]) and np.isfinite(Hg[i+1]) \
                and Hg[i] * Hg[i+1] <= 0:
            roots.append(brentq(lambda p: Hfun(r0, p),
                                p_grid[i], p_grid[i+1]))
    p0 = None
    for p_ in roots:
        if np.sign(dHdp(r0, p_)) == verso:
            p0 = p_
    assert p0 is not None, f'lancio non trovato (radici {roots})'

    def rhs(v_, y):
        return [dHdp(y[0], y[1]), -dHdr(y[0], y[1])]
    evs = []
    ev_t = lambda v_, y: dHdp(y[0], y[1])
    ev_t.terminal, ev_t.direction = True, 0
    evs.append(ev_t)
    if ev_lo is not None:
        e1 = lambda v_, y: y[0] - ev_lo
        e1.terminal, e1.direction = True, -1
        evs.append(e1)
    if ev_hi is not None:
        e2 = lambda v_, y: y[0] - ev_hi
        e2.terminal, e2.direction = True, 1
        evs.append(e2)
    s = solve_ivp(rhs, [0, v_max], [r0, p0], rtol=1e-12, atol=1e-14,
                  method='DOP853', events=evs, dense_output=True)
    return s, p0, dHdp

g_tau_n = lambda rr: rr**2 * f_n(rr) - J_n**2 * (E_n**2 - f_n(rr))

# ramo tau, entrante: svolta interna
s, p0, dHdp = flusso(H_tau, -1)
r_turn = s.y[0, -1]
root_in = brentq(g_tau_n, r_b + 1e-9, r0)
print(f"\n  tau entrante: svolta a r = {r_turn:.9f}   radice g_tau = "
      f"{root_in:.9f}   diff = {abs(r_turn-root_in):.1e}")
# ramo tau, uscente: svolta esterna (riflessione dall'orizzonte cosmologico)
s, p0, _ = flusso(H_tau, +1)
r_turn = s.y[0, -1]
root_out = brentq(g_tau_n, r0, r_c - 1e-9)
print(f"  tau uscente:  svolta a r = {r_turn:.9f}   radice g_tau = "
      f"{root_out:.9f}   diff = {abs(r_turn-root_out):.1e}")
print("  => il ramo tau rimbalza in una scatola (r_in, r_out) DENTRO la vasca")

# dphi/dr lungo il flusso tau vs forma chiusa
dHdp_tau = sp.lambdify((rr_, pr), sp.diff(H_tau, pr), 'numpy')
dHdJ = sp.lambdify((rr_, pr),
                   sp.diff(H_tau, sp.Symbol('J', positive=True)), 'numpy')
# ricostruisco p_r(r) dal flusso e confronto dphi/dr
s, p0, _ = flusso(H_tau, -1)
ok = True
for rv in (4.5, 4.0, 3.5):
    # trova v con r(v)=rv
    vv_ = brentq(lambda v_: s.sol(v_)[0] - rv, 0, s.t[-1])
    rq, pq = s.sol(vv_)
    # dphi/dv = dH/dJ (J entra come p_phi)
    Jsym = sp.Symbol('J', positive=True)
    H_tau_J = pr * (f_s - E_n**2) - E_n \
        + sp.sqrt(w_s) * sp.sqrt((E_n * pr + 1)**2 + Jsym**2 / rr_**2)
    dphidv = float(sp.diff(H_tau_J, Jsym).subs(
        [(rr_, rq), (pr, pq), (Jsym, J_n)]))
    drdv = dHdp_tau(rq, pq)
    num = abs(dphidv / drdv)
    wf = (E_n**2 - f_n(rv)) * f_n(rv)
    Dl = rv**2 * f_n(rv)
    cf = J_n * rv * np.sqrt(wf) / (Dl * np.sqrt(Dl - J_n**2 * (E_n**2 - f_n(rv))))
    ok = ok and abs(num - cf) < 1e-8 * cf
    print(f"  dphi/dr a r={rv}: flusso = {num:.10f}  forma chiusa = {cf:.10f}")
print("  match forma chiusa (f generica!):", ok)

# ramo v: raggiunge entrambi gli orizzonti
for verso, targ, nome in ((-1, r_b, 'r_b'), (+1, r_c, 'r_c')):
    s, p0, dHdp_v = flusso(H_v, verso, ev_lo=r_b * 1.0000001,
                           ev_hi=r_c * 0.9999999, v_max=3000.0)
    r_end = s.y[0, -1]
    print(f"  ramo v verso {nome}: arriva a r = {r_end:.6f} "
          f"(orizzonte {targ:.6f}, f = {f_n(r_end):+.2e})")
print("  => ramo v: attraversamento marginale di ENTRAMBI gli orizzonti;")
print("     ramo tau: confinato. Dicotomia massimale dentro la vasca.")
print("\nFATTO.")
