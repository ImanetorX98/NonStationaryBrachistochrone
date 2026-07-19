# -*- coding: utf-8 -*-
"""
Mappa di penetrazione della brachistocrona tau in Vaidya lineare.

Setup (parametro v, ' = d/dv, E fissa, m(v) = m0 + mu*v):
  - lancio ENTRANTE da r0 = 10 all'istante v0, con la mira localmente
    ottima p_v = 0 (famiglia brachistocrona; in statico e' esatta e
    riproduce l'orbita chiusa doranTau a=0);
  - integrazione forward; esiti:
      RIFLETTE:  r' = 0 con r > 2m(v)  (periasse esterno)
      CATTURATA: attraversa r = 2m(v) verso l'interno
  - soglia J_c(v0): bisezione in J tra i due esiti.

Riferimento statico (Kerr a=0, tricotomia doranTau): il ramo tau penetra
SOLO a J = J_c = a/E = 0: soglia nulla, penetrazione a misura zero.
Atteso non stazionario: J_c FINITO (accrescimento: l'orizzonte cresce
verso l'orbita e inghiotte periassi marginali; evaporazione: l'orizzonte
si ritira e riflette anche mire piu' strette... o viceversa — misuriamo).

Stima analitica small-J (statico): il periasse sta sopra l'orizzonte di
  delta = r_p - 2m ~ J^2 E^2 / (2m)
=> in accrescimento la cattura richiede che 2m cresca di ~delta durante
l'avvicinamento: J_c ~ sqrt(2m * Delta(2m))/E, finito.

Output: fig_vaidya_penetration_map (Jc vs v0; traiettorie di esempio
vicino alla soglia con inseguimento dell'orizzonte).
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

# ------------------------------------------------------------- sympy EOM
vv, rr, rp, q = sp.symbols('vv rr rp q', real=True)
E, J, mm, mp, mu = sp.symbols('E J mm mp mu', real=True)

f = 1 - 2 * mm / rr
Wv = (f - 2 * rp - (f - rp)**2 / E**2) / rr**2
F_tau = (f - rp) / E - J * sp.sqrt(Wv)
p_r = sp.diff(F_tau, rp)
dp_tot = (sp.diff(p_r, vv) + sp.diff(p_r, mm) * mp
          + sp.diff(p_r, rr) * rp + sp.diff(p_r, rp) * q)
EL = dp_tot - sp.diff(F_tau, rr)
rpp_sol = sp.solve(sp.Eq(EL, 0), q)[0]

m0v, E_n = 1.0, 1.2
r0 = 10.0
subs_n = [(mm, m0v + mu * vv), (mp, mu), (E, E_n)]
args = (vv, rr, rp, mu, J)
rpp_fn = sp.lambdify(args, rpp_sol.subs(subs_n), 'numpy')
Wv_fn = sp.lambdify(args, Wv.subs(subs_n), 'numpy')

# p_v della forma r-param per la mira al lancio (ramo entrante u < u_-)
u_ = sp.Symbol('u_', real=True)
W_r = (f * u_**2 - 2 * u_ - (f * u_ - 1)**2 / E**2) / rr**2
F_r = (f * u_ - 1) / E - J * sp.sqrt(W_r)
pv_fn = sp.lambdify((rr, vv, u_, mu, J),
                    sp.diff(F_r, u_).subs(subs_n), 'numpy')

def rp_lancio(mu_n, v0, J_n):
    """r'(v0) entrante: mira p_v=0 (m congelata) + riflessione temporale.

    p_v = 0  <=>  sqrt(W) = J E W'/(2f)  <=>  N(u) = K (2 a2 u + a1)^2
    con N = a2 u^2 + a1 u + a0 (numeratore di W_r) e K = J^2E^2/(4f^2r^2):
    quadratica esatta, robusta anche a J piccolissimo. Poi riflessione
    temporale u_in = 2/f - u_out (vertice di W in u = 1/f).
    """
    fl = 1 - 2 * (m0v + mu_n * v0) / r0
    a2 = fl - fl**2 / E_n**2
    a1 = -2 + 2 * fl / E_n**2
    a0 = -1 / E_n**2
    K = J_n**2 * E_n**2 / (4 * fl**2 * r0**2)
    A = a2 * (1 - 4 * K * a2)
    B = a1 * (1 - 4 * K * a2)
    C = a0 - K * a1**2
    u_out = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)
    u_in = 2.0 / fl - u_out
    return 1.0 / u_in

def lancio(mu_n, v0, J_n, v_max_extra=400.0):
    """Integra forward; ritorna esito e diagnostica."""
    rp0 = rp_lancio(mu_n, v0, J_n)
    v_end = v0 + v_max_extra
    if mu_n < 0:                # m(v) > 0.03
        v_end = min(v_end, (m0v - 0.03) / (-mu_n))

    def rhs(v_, y):
        return [y[1], rpp_fn(v_, y[0], y[1], mu_n, J_n)]

    ev_turn = lambda v_, y: y[1]
    ev_turn.terminal, ev_turn.direction = True, 1        # r': - -> +
    ev_hor = lambda v_, y: y[0] - 2 * (m0v + mu_n * v_)
    ev_hor.terminal, ev_hor.direction = True, -1         # entra
    ev_W = lambda v_, y: Wv_fn(v_, y[0], y[1], mu_n, J_n) - 1e-12
    ev_W.terminal, ev_W.direction = True, -1
    s = solve_ivp(rhs, [v0, v_end], [r0, rp0], rtol=1e-12, atol=1e-14,
                  method='DOP853', events=[ev_turn, ev_hor, ev_W],
                  dense_output=True)
    if len(s.t_events[1]):
        return dict(esito='cattura', v_ev=s.t_events[1][0],
                    r_ev=s.y_events[1][0][0], sol=s)
    if len(s.t_events[0]):
        y = s.y_events[0][0]
        return dict(esito='riflette', v_ev=s.t_events[0][0], r_ev=y[0],
                    r_min=y[0], sol=s)
    if len(s.t_events[2]):
        # radializzazione W->0: dentro l'orizzonte = plunge catturato
        v_e, r_e = s.t_events[2][0], s.y_events[2][0][0]
        dentro = r_e < 2 * (m0v + mu_n * v_e)
        return dict(esito='cattura' if dentro else 'W=0-fuori',
                    v_ev=v_e, r_ev=r_e, sol=s)
    return dict(esito='indeterminato', v_ev=s.t[-1], r_ev=s.y[0, -1], sol=s)

def J_critico(mu_n, v0, J_hi=1.6, tol=2e-5):
    """J_c = sup{J : cattura}. Scala di sonde + bisezione.

    Nota: a J -> 0 la riduzione degenera (d2F/dr'^2 ∝ J): sotto J ~ 0.01
    l'EL e' stiff e gli esiti numerici non sono affidabili; le sonde
    partono da 0.01. J_c = 0 significa 'nessuna cattura sopra 0.01'.
    """
    probes = [0.01, 0.03, 0.1, 0.3]
    lo = None
    for p in probes:
        if lancio(mu_n, v0, p)['esito'] == 'cattura':
            lo = p
            break
    if lo is None:
        return 0.0
    e_hi = lancio(mu_n, v0, J_hi)['esito']
    assert e_hi == 'riflette', f'J_hi={J_hi} non riflette ({e_hi})'
    hi = J_hi
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if lancio(mu_n, v0, mid)['esito'] == 'cattura':
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

# --------------------------------------------------------------- studio
print("=" * 72)
print("[V1] validazione statica: J piccolo riflette (soglia J_c -> 0)")
print("=" * 72)
for J_t in (0.05, 0.02, 0.005):
    o = lancio(0.0, 0.0, J_t)
    extra = f", r_min = {o['r_min']:.6f}" if o['esito'] == 'riflette' else ''
    print(f"  statico J={J_t:6.3f}: {o['esito']}{extra}"
          f"   [periasse teorico 2m + J^2E^2/(2m) ~ "
          f"{2 + J_t**2 * E_n**2 / 2:.6f}]")
Jc_st = J_critico(0.0, 0.0)
print(f"  J_c statico (bisezione) = {Jc_st:.2e}")
print("  [floor numerico: il margine r_p - 2m = J^2 E^2/(2m) scende sotto")
print("   la risoluzione dell'integratore; teoricamente J_c = 0, misura nulla]")

print()
print("=" * 72)
print("[V2] mappa J_c(v0) per mu = +0.01 (accresce), -0.01 (evapora)")
print("=" * 72)
v0_grid = np.linspace(0.0, 60.0, 9)
mappe = {}
for mu_n, lab in ((0.01, 'accresce'), (-0.01, 'evapora')):
    Jcs = []
    for v0 in v0_grid:
        try:
            Jcs.append(J_critico(mu_n, v0))
        except AssertionError:
            Jcs.append(np.nan)
        print(f"  mu={mu_n:+.2f}  v0={v0:5.1f}:  J_c = {Jcs[-1]:.5f}")
    mappe[mu_n] = np.array(Jcs)

print()
print("=" * 72)
print("[V3] scala small-J (accrescimento): J_c vs sqrt(mu)")
print("=" * 72)
mu_scan = [0.0025, 0.005, 0.01, 0.02]
Jc_mu = [J_critico(m_, 0.0) for m_ in mu_scan]
for m_, jc in zip(mu_scan, Jc_mu):
    print(f"  mu={m_:7.4f}:  J_c = {jc:.5f}   J_c/sqrt(mu) = "
          f"{jc/np.sqrt(m_):.4f}")

# --------------------------------------------------------------- figure
def salva(fig, nome):
    savefig(fig, OUT, nome)

fig, (axa, axb) = plt.subplots(2, 1, figsize=(COL, 5.6))
axa.plot(v0_grid, mappe[0.01], 'C3o-', ms=4, label='accreting $\\mu=+0.01$')
axa.plot(v0_grid, mappe[-0.01], 'C0o-', ms=4, label='evaporating $\\mu=-0.01$')
axa.axhline(0, color='k', lw=1.0, label='static: $J_c=0$ (zero measure)')
axa.set_xlabel('launch instant $v_0$')
axa.set_ylabel('$J_c$')
axa.set_title('penetration threshold $J_c(v_0)$: below $J_c$ the\n'
              '$\\tau$-brachistochrone crosses $r=2m(v)$ ($E=1.2$, $r_0=10$)')
axa.legend()

# esempi vicino alla soglia (accrescimento, v0=0)
Jc0 = mappe[0.01][0]
v_ev_max = 0.0
for J_t, col, lab in ((Jc0 * 0.97, 'C3', 'below threshold: captured'),
                      (Jc0 * 1.03, 'C0', 'above threshold: reflects')):
    o = lancio(0.01, 0.0, J_t, v_max_extra=60.0)
    vg = np.linspace(o['sol'].t[0], o['v_ev'], 500)
    rg = o['sol'].sol(vg)[0]
    axb.plot(vg, rg, color=col, label=f'$J={J_t:.4f}$: {lab}')
    v_ev_max = max(v_ev_max, o['v_ev'])
vg = np.linspace(0, v_ev_max + 5, 200)
axb.plot(vg, 2 * (m0v + 0.01 * vg), 'k:', lw=1.2,
         label='apparent horizon $2m(v)$')
axb.set_xlabel('$v$')
axb.set_ylabel('$r(v)$')
axb.set_ylim(1.8, 5.0)
axb.set_xlim(0, v_ev_max + 5)
axb.set_title(f'horizon chase near the threshold\n'
              f'($\\mu=+0.01$, $v_0=0$, $J_c={Jc0:.4f}$)')
axb.legend()
salva(fig, 'fig_vaidya_penetration_map')
print('\nFATTO.')
