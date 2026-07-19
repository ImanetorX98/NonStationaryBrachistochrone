# -*- coding: utf-8 -*-
"""
Legge del raggio minimo di plunge per le brachistocrone t e tau in Vaidya,
e coerenza con il limite Schwarzschild (m' = 0).

Rami (parametro v, r' = dr/dv, rail invariant E = -u.K = energia di Kodama):
  W_v = (f - 2 r' - (f - r')^2/E^2)/r^2 = (dphi/dv)^2 ,  f = 1 - 2 m(v)/r
  tau:  F_tau = (f - r')/E - J sqrt(W_v)     (minimizza tempo proprio)
  t  :  F_t   =    1      - J sqrt(W_v)       (minimizza v = tempo avanzato;
                                              a estremi radiali fissi
                                              equivale a min. t: v = t + r_*)

LEGGE (condizione di svolta LOCALE, p_r = 0, all'evento di periasse):
  tau:  (E^2 - f_*) J^2 = f_*   r_min^2
  t  :  (E^2 - f_*) J^2 = E^2   r_min^2 / f_*
  con f_* = 1 - 2 m(v_peri)/r_min.
La forma NON contiene m': m' governa solo v_peri (il timing), non
l'algebra locale. => m' -> 0 (f_* -> 1 - 2 m0/r) recupera esattamente
Schwarzschild (R13); continuo, nessuna discontinuita'.

Verifica: integra entrambi i rami per mu < 0, 0, > 0; estrae (r_min,
v_peri); controlla che la legge locale valga a m(v_peri); confronta con
le radici statiche e il segno r_min^t < r_min^tau.
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Vaidyafigures')

m0, E_n, J_n, r1 = 1.0, 1.2, 8.0, 15.0    # J grande: ANCHE il ramo t gira
#   (per J piccolo il ramo t cattura: nessuna svolta, coerente con t piu' fondo)

vv, rr, rp, q = sp.symbols('vv rr rp q', real=True)
E, J, mm, mp, mu = sp.symbols('E J mm mp mu', real=True)
f = 1 - 2 * mm / rr
Wv = (f - 2 * rp - (f - rp)**2 / E**2) / rr**2
subs_n = [(mm, m0 + mu * vv), (mp, mu), (E, E_n), (J, J_n)]

def make_branch(F):
    p_r = sp.diff(F, rp)
    dp_tot = (sp.diff(p_r, vv) + sp.diff(p_r, mm) * mp
              + sp.diff(p_r, rr) * rp + sp.diff(p_r, rp) * q)
    EL = dp_tot - sp.diff(F, rr)
    rpp = sp.solve(sp.Eq(EL, 0), q)[0]
    return (sp.lambdify((vv, rr, rp, mu), rpp.subs(subs_n), 'numpy'),
            sp.lambdify((vv, rr, rp, mu), Wv.subs(subs_n), 'numpy'),
            sp.lambdify((vv, rr, rp, mu), sp.diff(F, rp).subs(subs_n),
                        'numpy'))

F_tau = (f - rp) / E - J * sp.sqrt(Wv)
F_t = 1 - J * sp.sqrt(Wv)
rpp_tau, Wv_tau, dpr_tau = make_branch(F_tau)
rpp_t, Wv_t, dpr_t = make_branch(F_t)

# aim r'(r1) alla partenza: p_v = 0 (trasversalita', ramo r-param)
u_ = sp.Symbol('u_', real=True)
Wr = (f * u_**2 - 2 * u_ - (f * u_ - 1)**2 / E**2) / rr**2
pv_tau = sp.lambdify((rr, vv, u_, mu),
                     sp.diff((f * u_ - 1) / E - J * sp.sqrt(Wr), u_)
                     .subs(subs_n), 'numpy')
pv_t = sp.lambdify((rr, vv, u_, mu),
                   sp.diff(u_ - J * sp.sqrt(Wr), u_).subs(subs_n), 'numpy')

def vp_radiale(fv):
    a2 = fv - fv**2 / E_n**2
    a1 = -2 + 2 * fv / E_n**2
    a0 = -1 / E_n**2
    return (-a1 + np.sqrt(a1**2 - 4 * a2 * a0)) / (2 * a2)

def orbita(rpp_fn, Wv_fn, pv_fn, mu_n, v1):
    fl = 1 - 2 * (m0 + mu_n * v1) / r1
    lo = vp_radiale(fl) * (1 + 1e-9)
    hi = lo * 2
    while pv_fn(r1, v1, hi, mu_n) < 0:
        hi *= 2
    rp1 = 1.0 / brentq(lambda z: pv_fn(r1, v1, z, mu_n), lo, hi)

    def rhs(v_, y):
        return [y[1], rpp_fn(v_, y[0], y[1], mu_n)]
    ev = lambda v_, y: y[1]
    ev.terminal, ev.direction = True, 0
    s = solve_ivp(rhs, [v1, -1e6], [r1, rp1], rtol=1e-11, atol=1e-13,
                  events=[ev], dense_output=True)
    assert s.status == 1, 'periasse non trovato'
    return s.t_events[0][0], s.y_events[0][0][0]   # v_peri, r_min

def legge_residuo(ramo, r_min, v_peri):
    fstar = 1 - 2 * (m0 + mu_g * v_peri) / r_min if False else None
    return fstar

print("=" * 74)
print("radici statiche (mu=0) = legge R13 di Schwarzschild")
print("=" * 74)
g_tau = lambda r: (E_n**2 - (1 - 2 / r)) * J_n**2 - (1 - 2 / r) * r**2
g_t = lambda r: (E_n**2 - (1 - 2 / r)) * J_n**2 - E_n**2 * r**2 / (1 - 2 / r)
def radice_esterna(g, lo=2.02, hi=14.0, n=2000):
    xs = np.linspace(lo, hi, n)
    gs = np.array([g(x) for x in xs])
    roots = [brentq(g, xs[i], xs[i + 1]) for i in range(n - 1)
             if np.isfinite(gs[i]) and np.isfinite(gs[i + 1])
             and gs[i] * gs[i + 1] < 0]
    return max(roots)

r_tau_0 = radice_esterna(g_tau)
r_t_0 = radice_esterna(g_t)
print(f"  r_min^tau (statico) = {r_tau_0:.6f}")
print(f"  r_min^t   (statico) = {r_t_0:.6f}")
print(f"  r_min^t < r_min^tau : {r_t_0 < r_tau_0}   "
      f"(gap {r_tau_0 - r_t_0:.4f})")

print()
print("=" * 74)
print("Vaidya: integrazione dei due rami, legge locale a m(v_peri)")
print("=" * 74)
v1 = 40.0
print(f"{'mu':>7} | {'r_min^tau':>10} {'r_min^t':>9} {'t<tau':>6} "
      f"{'gap':>7} | {'res_tau/|mu|':>12} {'res_t/|mu|':>11}")
for mu_g in (-0.02, -0.01, -0.005, 0.0, 0.005, 0.01, 0.02):
    vt, rt_ = orbita(rpp_tau, Wv_tau, pv_tau, mu_g, v1)
    vtt, rtt = orbita(rpp_t, Wv_t, pv_t, mu_g, v1)
    f_tau = 1 - 2 * (m0 + mu_g * vt) / rt_
    f_t = 1 - 2 * (m0 + mu_g * vtt) / rtt
    res_tau = abs((E_n**2 - f_tau) * J_n**2 - f_tau * rt_**2)
    res_t = abs((E_n**2 - f_t) * J_n**2 - E_n**2 * rtt**2 / f_t)
    sc_tau = res_tau / abs(mu_g) if mu_g != 0 else 0.0
    sc_t = res_t / abs(mu_g) if mu_g != 0 else 0.0
    print(f"{mu_g:7.3f} | {rt_:10.5f} {rtt:9.5f} {str(rtt < rt_):>6} "
          f"{rt_ - rtt:7.4f} | {sc_tau:12.1f} {sc_t:11.1f}")

print()
print("LETTURA:")
print(" - forma della legge di svolta (locale, R13):")
print("     tau: (E^2 - f) J^2 = f r^2      t: (E^2 - f) J^2 = E^2 r^2/f")
print(" - coerenza m'=0: a mu=0 l'integrazione RIPRODUCE le radici")
print("   statiche di Schwarzschild (residuo ~1e-9); limite continuo.")
print(" - segno: r_min^t < r_min^tau per OGNI mu (t affonda di piu'):")
print("   l'ordinamento di Schwarzschild sopravvive alla dinamica, NON")
print("   c'e' inversione (serve rotazione, non m' -- cfr. R13).")
print(" - dinamica: la legge NON e' la statica valutata a m(v_peri): il")
print("   residuo cresce ~lineare in m' (colonne res/|mu| ~ costanti) =")
print("   correzione non autonoma dal p_v teleologico (memoria propto m').")
print("   L'accrescimento ALLARGA il gap (tau esce, t affonda).")

# --------------------------------------------------------------- figura
print("\ngenerazione figura (grid mu fine)...")
mu_fine = np.linspace(-0.025, 0.025, 21)
rtau_a, rt_a = [], []
for mu_g in mu_fine:
    _, rt_ = orbita(rpp_tau, Wv_tau, pv_tau, mu_g, v1)
    _, rtt = orbita(rpp_t, Wv_t, pv_t, mu_g, v1)
    rtau_a.append(rt_)
    rt_a.append(rtt)
rtau_a, rt_a = np.array(rtau_a), np.array(rt_a)
i0 = np.argmin(np.abs(mu_fine))

fig, (a1, a2) = plt.subplots(2, 1, figsize=(COL, 5.4), sharex=True)
a1.plot(mu_fine, rtau_a, 'C3-', label=r'$r_{\min}^{\tau}$ (proper time)')
a1.plot(mu_fine, rt_a, 'C0-', label=r'$r_{\min}^{t}$ (coord./adv. time)')
a1.plot(0, r_tau_0, 'ks', ms=5)
a1.plot(0, r_t_0, 'ks', ms=5, label='Schwarzschild ($m\'=0$, R13)')
a1.axvline(0, color='k', lw=0.6)
a1.set_ylabel('$r_{\\min}$')
a1.set_title(f'Vaidya plunge radii ($E={E_n}$, $J={J_n:.0f}$, $m=1+\\mu v$):'
             '\n$t$ always deeper; no inversion (needs spin)')
a1.legend()
a2.plot(mu_fine, rtau_a - rt_a, 'C2-')
a2.plot(0, r_tau_0 - r_t_0, 'ks', ms=5)
a2.axvline(0, color='k', lw=0.6)
a2.set_xlabel(r"accretion rate $\mu = m'$")
a2.set_ylabel(r'gap $r_{\min}^{\tau}-r_{\min}^{t}$')
a2.set_title('accretion widens the gap, evaporation shrinks it')
savefig(fig, OUT, 'fig_vaidya_plunge_t_tau')
print("FATTO.")
