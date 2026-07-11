# -*- coding: utf-8 -*-
"""
Brachistocrona tau in Vaidya — parametro v (regolare al periasse).

La parametrizzazione in r degenera al periasse (dv/dr -> oo). Con v come
parametro (' = d/dv):

    Lam = dtau/dv = sqrt(f - 2 r' - r^2 phi'^2)
    vincolo -u_v = E  =>  Lam = (f - r')/E
    =>  phi'^2 = W_v(v, r, r') = [ f - 2 r' - (f - r')^2/E^2 ] / r^2

Funzionale ridotto (J = momento di Fermat, costante):

    F_tau = (f - r')/E - J sqrt(W_v)       T_tau = int (f - r')/E dv

Identico al funzionale in parametro r via cambio di variabile
((f v'-1) dr = (f - r') dv), quindi STESSE curve estremali — ma l'EL in
v(r) resta regolare a r'=0: si integra ATTRAVERSO il rimbalzo.

Trasversalita' all'arrivo (istante v1 libero): presa dalla forma r-param
gia' validata, p_v(r1) = 0, che fissa r'(v1) = 1/v'(r1).

Verifiche / studio:
  V1  EL in parametro v: r'' esplicito; forzante dF/dv = (dF/dm) m'
  V2  cross-check r-param: stessi (v, Delta_phi, T_tau) al passaggio
      r = 3.5 per Schwarzschild e Vaidya mu=0.01 (tabella VaidyaResults R4)
  V3  rimbalzo: integrazione attraverso il periasse; statico: r_min
      coincide con la radice esatta di  r^2 f - J^2 (E^2 - f) = 0
      (tricotomia doranTau con a=0); Vaidya: periasse ~ dove la r-param
      moriva (r ~ 3.0105)
  V4  scan di timing: r_min in funzione dell'istante d'arrivo v1 per
      mu = +0.01 (accresce), 0 (statico), -0.01 (evapora):
      il caso statico e' piatto, i dinamici NO — prima misura della
      "finestra temporale" genuinamente non stazionaria.
Figure: Vaidyafigures/fig_vaidya_bounce, fig_vaidya_timing.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig

set_style()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Vaidyafigures')
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------- V1: sympy
vv, rr, rp, q = sp.symbols('vv rr rp q', real=True)     # v, r, r', r''
E, J, mm, mp, mu = sp.symbols('E J mm mp mu', real=True)

f = 1 - 2 * mm / rr
Wv = (f - 2 * rp - (f - rp)**2 / E**2) / rr**2     # phi'^2 = W_v
F_tau = (f - rp) / E - J * sp.sqrt(Wv)

p_r = sp.diff(F_tau, rp)
dp_tot = (sp.diff(p_r, vv) + sp.diff(p_r, mm) * mp
          + sp.diff(p_r, rr) * rp + sp.diff(p_r, rp) * q)
EL = dp_tot - (sp.diff(F_tau, rr) + sp.diff(F_tau, mm) * 0)
# ATTENZIONE: dF/dr e' derivata esplicita in r; la dipendenza da v e' via
# m(v) e sta in dp_tot (termine mp). dF/dv = (dF/dm) m' entra con segno
# opposto? No: EL = d/dv(dF/dr') - dF/dr; dF/dr non contiene m'.
rpp_sol = sp.solve(sp.Eq(EL, 0), q)[0]

print("=" * 72)
print("[V1] EL in parametro v")
print("=" * 72)
prop = sp.simplify(sp.diff(F_tau, vv)) == 0
print("  F non dipende esplicitamente da v (solo via m):", prop)
print("  dF/dm =", sp.factor_terms(sp.radsimp(sp.diff(F_tau, mm))))
print("  r'' risolto (lineare in q):  ok, regolare dove W_v > 0")

m0v, E_n, J_n = 1.0, 1.2, 1.3
r1 = 10.0
subs_n = [(mm, m0v + mu * vv), (mp, mu), (E, E_n), (J, J_n)]
args = (vv, rr, rp, mu)
rpp_fn = sp.lambdify(args, rpp_sol.subs(subs_n), 'numpy')
Wv_fn = sp.lambdify(args, Wv.subs(subs_n), 'numpy')
Lam_fn = sp.lambdify(args, ((f - rp) / E).subs(subs_n), 'numpy')

# trasversalita' dalla forma r-param (validata): p_v(r1) = 0
u_ = sp.Symbol('u_', real=True)
f_r = 1 - 2 * mm / rr
W_r = (f_r * u_**2 - 2 * u_ - (f_r * u_ - 1)**2 / E**2) / rr**2
F_r = (f_r * u_ - 1) / E - J * sp.sqrt(W_r)
pv_fn = sp.lambdify((rr, vv, u_, mu),
                    sp.diff(F_r, u_).subs(subs_n), 'numpy')

def vp_radiale(f_val):
    a2 = f_val - f_val**2 / E_n**2
    a1 = -2 + 2 * f_val / E_n**2
    a0 = -1 / E_n**2
    return (-a1 + np.sqrt(a1**2 - 4 * a2 * a0)) / (2 * a2)

def rp_arrivo(mu_n, v1):
    fl = 1 - 2 * (m0v + mu_n * v1) / r1
    lo = vp_radiale(fl) * (1 + 1e-9)
    hi = lo * 2
    while pv_fn(r1, v1, hi, mu_n) < 0:
        hi *= 2
    return 1.0 / brentq(lambda z: pv_fn(r1, v1, z, mu_n), lo, hi)

def orbita(mu_n, v1, v_min=-1e6):
    """Backward da (v1, r1) attraverso il periasse fino a r = r1."""
    rp1 = rp_arrivo(mu_n, v1)

    def rhs(v_, y):
        return [y[1], rpp_fn(v_, y[0], y[1], mu_n),
                np.sqrt(max(Wv_fn(v_, y[0], y[1], mu_n), 0.0)),
                Lam_fn(v_, y[0], y[1], mu_n)]

    ev_peri = lambda v_, y: y[1]
    ev_peri.terminal, ev_peri.direction = True, 0
    ev_r35 = lambda v_, y: y[0] - 3.5
    ev_r35.terminal = False
    s1 = solve_ivp(rhs, [v1, v_min], [r1, rp1, 0.0, 0.0],
                   rtol=1e-11, atol=1e-13, dense_output=True,
                   events=[ev_peri, ev_r35])
    assert s1.status == 1, 'periasse non trovato'
    v_p, y_p = s1.t_events[0][0], s1.y_events[0][0]
    ev_r10 = lambda v_, y: y[0] - r1
    ev_r10.terminal, ev_r10.direction = True, 0
    s2 = solve_ivp(rhs, [v_p, v_min], y_p, rtol=1e-11, atol=1e-13,
                   dense_output=True, events=[ev_r10])
    out = dict(v_peri=v_p, r_min=y_p[0], sol1=s1, sol2=s2)
    if len(s1.t_events[1]):        # passaggio a r=3.5 (ramo entrante? no:
        ve = s1.t_events[1][0]     # backward: e' il ramo USCENTE finale)
        ye = s1.y_events[1][0]
        out.update(v_35=ve, dphi_35=-ye[2], Ttau_35=-ye[3])
    if s2.status == 1:
        v0, y0 = s2.t_events[0][0], s2.y_events[0][0]
        out.update(v_start=v0, dphi_tot=-y0[2], Ttau_tot=-y0[3])
    return out

print()
print("=" * 72)
print("[V2] cross-check con la parametrizzazione in r (tabella R4)")
print("=" * 72)
rif = {0.0: (16.103201, 0.401351, 7.757186),
       0.01: (13.043582, 0.529307, 7.479964)}
for mu_n in (0.0, 0.01):
    o = orbita(mu_n, 40.0)
    v35, d35, T35 = o['v_35'], o['dphi_35'], o['Ttau_35']
    vr, dr_, Tr = rif[mu_n]
    print(f"  mu={mu_n:5.2f}:  v(3.5) = {v35:.6f}  (r-param {vr:.6f}, "
          f"diff {abs(v35-vr):.2e})")
    print(f"            Dphi   = {d35:.6f}  (r-param {dr_:.6f}, "
          f"diff {abs(d35-dr_):.2e})")
    print(f"            T_tau  = {T35:.6f}  (r-param {Tr:.6f}, "
          f"diff {abs(T35-Tr):.2e})")

print()
print("=" * 72)
print("[V3] rimbalzo attraverso il periasse")
print("=" * 72)
# statico: radice esatta di r^2 f - J^2 (E^2 - f) = 0
g_static = lambda r_: r_**2 * (1 - 2 / r_) - J_n**2 * (E_n**2 - 1 + 2 / r_)
r_p_exact = brentq(g_static, 2.01, 6.0)
o_s = orbita(0.0, 40.0)
print(f"  statico:  r_min integrato = {o_s['r_min']:.8f}   "
      f"radice esatta = {r_p_exact:.8f}   diff = "
      f"{abs(o_s['r_min']-r_p_exact):.2e}")
o_v = orbita(0.01, 40.0)
print(f"  Vaidya mu=0.01: r_min = {o_v['r_min']:.6f} a v = "
      f"{o_v['v_peri']:.4f}  (r-param moriva a r=3.0105)")
print(f"  orbita completa (r1 -> r_min -> r1): statico  v in "
      f"[{o_s['v_start']:.3f}, 40], Dphi = {o_s['dphi_tot']:.5f}, "
      f"T_tau = {o_s['Ttau_tot']:.5f}")
print(f"                                       mu=0.01  v in "
      f"[{o_v['v_start']:.3f}, 40], Dphi = {o_v['dphi_tot']:.5f}, "
      f"T_tau = {o_v['Ttau_tot']:.5f}")

print()
print("=" * 72)
print("[V4] scan di timing: r_min vs istante d'arrivo v1")
print("=" * 72)
mus = [(-0.01, 'C0', 'evapora $\\mu=-0.01$'),
       (0.0, 'k', 'statico $\\mu=0$'),
       (0.01, 'C3', 'accresce $\\mu=+0.01$')]
v1_grid = np.linspace(25.0, 70.0, 10)
curve = {}
for mu_n, _, _ in mus:
    rmins, vperis = [], []
    for v1 in v1_grid:
        o = orbita(mu_n, v1)
        rmins.append(o['r_min'])
        vperis.append(o['v_peri'])
    curve[mu_n] = (np.array(rmins), np.array(vperis))
    print(f"  mu={mu_n:+.2f}:  r_min: {curve[mu_n][0][0]:.4f} (v1={v1_grid[0]:.0f})"
          f"  ->  {curve[mu_n][0][-1]:.4f} (v1={v1_grid[-1]:.0f})")

# --------------------------------------------------------------- figure
def salva(fig, nome):
    savefig(fig, OUT, nome)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COL, 5.8))
for mu_n, col, lab in mus:
    o = orbita(mu_n, 40.0)
    vg1 = np.linspace(o['sol1'].t[0], o['sol1'].t[-1], 500)
    vg2 = np.linspace(o['sol2'].t[0], o['sol2'].t[-1], 500)
    vg = np.concatenate([vg2[::-1], vg1[::-1]])
    rg = np.concatenate([o['sol2'].sol(vg2)[0][::-1],
                         o['sol1'].sol(vg1)[0][::-1]])
    phg = np.concatenate([o['sol2'].sol(vg2)[2][::-1],
                          o['sol1'].sol(vg1)[2][::-1]])
    ax1.plot(vg, rg, color=col, label=lab)
    ax1.plot(o['v_peri'], o['r_min'], 'o', color=col, ms=4)
    m_of_v = m0v + mu_n * vg
    ax1.plot(vg, 2 * m_of_v, color=col, ls=':', lw=0.8)
    ph0 = phg[0]
    ax2.plot(rg * np.cos(phg - ph0), rg * np.sin(phg - ph0), color=col,
             label=lab)
ax1.set_xlabel('$v$')
ax1.set_ylabel('$r(v)$')
ax1.set_title('bounce through periapsis ($v$ parameter)\n'
              'dotted: apparent horizon $2m(v)$')
ax1.legend()
ax2.set_aspect('equal')
ax2.set_xlabel('$x$')
ax2.set_ylabel('$y$')
ax2.set_title('orbits ($E=1.2$, $J=1.3$, arrival $v=40$, $r_1=10$)')
ax2.legend(loc='lower right')
salva(fig, 'fig_vaidya_bounce')

fig, (axa, axb) = plt.subplots(2, 1, figsize=(COL, 5.6))
for mu_n, col, lab in mus:
    rmins, vperis = curve[mu_n]
    axa.plot(v1_grid, rmins, 'o-', color=col, ms=4, label=lab)
    axb.plot(v1_grid, rmins / (2 * (m0v + mu_n * vperis)), 'o-', color=col,
             ms=4, label=lab)
axa.axhline(r_p_exact, color='k', lw=0.6, ls='--')
axa.set_xlabel('arrival instant $v_1$')
axa.set_ylabel('$r_{min}$')
axa.set_title('timing window: periapsis depends on WHEN you arrive\n'
              '(static: flat = exact root, dashed)')
axa.legend()
axb.set_xlabel('arrival instant $v_1$')
axb.set_ylabel('$r_{min} / 2m(v_{peri})$')
axb.set_title('penetration relative to the apparent horizon\n'
              'at periapsis')
axb.legend()
salva(fig, 'fig_vaidya_timing')
print('\nFATTO.')
