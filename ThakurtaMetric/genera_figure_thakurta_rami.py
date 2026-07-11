# -*- coding: utf-8 -*-
"""
Figura R8a — Thakurta-Kerr: una ellisse, tre costi, due fisiche.

  fig_thakurta_kerr_rami (Thakurtafigures/):
   (a) orbite di rimbalzo con J = ±1.3 (A=1, base Kerr):
       ramo tau  -> coppia SPECULARE (Riemanniana pura, niente vento)
       ramo eta  -> coppia ASIMMETRICA (Randers, vento gravitomagnetico)
       test implicito dell'intera catena: la simmetria non e' imposta,
       emerge dal flusso hamiltoniano con p_phi -> -p_phi.
   (b) i due scalari conformi che respirano sulla geometria rigida:
       n(eta,r) = Ehat/sqrt(Ehat^2 - A^2 f)   (indice del ramo eta/t)
       k_tau(eta,r) = A^2 f/(Ehat vbar)       (fattore del ramo tau)
       per A = 1, 1.6, 2.2 con asintoto al congelamento A^2 f = Ehat^2.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from paper_style import COL, set_style, savefig

set_style()
OUT = os.path.join(HERE, 'Thakurtafigures')

M_n, s_n, E_n, Jabs, r1 = 1.0, 0.9, 1.2, 1.3, 8.0

r, pr, pphi = sp.symbols('r p_r p_phi', real=True)
f = 1 - 2 * M_n / r
Dl = r**2 - 2 * M_n * r + s_n**2
P = r**2 + s_n**2 + 2 * M_n * s_n**2 / r
vb2 = 1 - f / E_n**2                      # A = 1
Pb = P + (2 * M_n * s_n / r)**2 / E_n**2
php0 = (2 * M_n * s_n / r) * vb2 / Pb
R_ = sp.sqrt(Dl * vb2 / Pb)

H_eta = pphi * php0 + R_ * sp.sqrt((Dl / r**2) * pr**2 + pphi**2 / Pb) - 1
ptphi = pphi - 2 * M_n * s_n / (r * E_n)
H_tau = ptphi * php0 + R_ * sp.sqrt((Dl / r**2) * pr**2
                                    + ptphi**2 / Pb) - f / E_n

r_plus = M_n + np.sqrt(M_n**2 - s_n**2)

def rimbalzo(Hx, Jv):
    """Backward dall'arrivo: rimbalzo (periasse) o spirale sull'orizzonte."""
    Hn = Hx.subs(pphi, Jv)
    dHdp = sp.lambdify((r, pr), sp.diff(Hn, pr), 'numpy')
    dHdr = sp.lambdify((r, pr), sp.diff(Hn, r), 'numpy')
    dHdJ = sp.lambdify((r, pr), sp.diff(Hx, pphi).subs(pphi, Jv), 'numpy')
    Hf = sp.lambdify((r, pr), Hn, 'numpy')
    pg = np.linspace(-120, 120, 480001)
    with np.errstate(invalid='ignore'):
        Hg = Hf(r1, pg)
    roots = [brentq(lambda p: Hf(r1, p), pg[i], pg[i + 1])
             for i in range(len(pg) - 1)
             if np.isfinite(Hg[i]) and np.isfinite(Hg[i + 1])
             and Hg[i] * Hg[i + 1] <= 0]
    p0 = [p_ for p_ in roots if dHdp(r1, p_) > 0][0]   # arrivo uscente

    def rhs(e_, y):
        return [dHdp(y[0], y[1]), -dHdr(y[0], y[1]), dHdJ(y[0], y[1])]

    ev_p = lambda e_, y: dHdp(y[0], y[1])
    ev_p.terminal, ev_p.direction = True, 0
    ev_h = lambda e_, y: y[0] - r_plus * 1.02
    ev_h.terminal, ev_h.direction = True, -1
    s1 = solve_ivp(rhs, [0, -400], [r1, p0, 0.0], rtol=1e-11, atol=1e-13,
                   method='DOP853', events=[ev_p, ev_h], dense_output=True)
    tt1 = np.linspace(0, s1.t[-1], 500)
    if len(s1.t_events[1]):          # spirale sull'orizzonte (ramo eta)
        rr_ = s1.sol(tt1)[0][::-1]
        ff_ = s1.sol(tt1)[2][::-1]
        dphi_tot = -s1.y_events[1][0][2]
        return rr_, ff_ - ff_[0], dphi_tot, s1.y_events[1][0][0], 'spirale'
    y_p = s1.y_events[0][0]
    ev_r = lambda e_, y: y[0] - r1
    ev_r.terminal, ev_r.direction = True, 1
    s2 = solve_ivp(rhs, [s1.t_events[0][0], -800], y_p, rtol=1e-11,
                   atol=1e-13, method='DOP853', events=[ev_r],
                   dense_output=True)
    tt2 = np.linspace(s2.t[0], s2.t[-1], 500)
    rr_ = np.concatenate([s2.sol(tt2)[0][::-1], s1.sol(tt1)[0][::-1]])
    ff_ = np.concatenate([s2.sol(tt2)[2][::-1], s1.sol(tt1)[2][::-1]])
    dphi_tot = -s2.y_events[0][0][2]
    return rr_, ff_ - ff_[0], dphi_tot, y_p[0], 'rimbalzo'

print('pannello (a): rimbalzi J = +-1.3, rami tau ed eta...')
fig, (axa, axb) = plt.subplots(2, 1, figsize=(COL, 6.2))
res = {}
for Hx, nome, cols in ((H_tau, 'tau', ('C3', 'C1')),
                       (H_eta, 'eta', ('C0', 'C9'))):
    for Jv, col, ls in ((+Jabs, cols[0], '-'), (-Jabs, cols[1], '--')):
        rr_, ff_, dphi, rmin, tipo = rimbalzo(Hx, Jv)
        res[(nome, Jv)] = (dphi, rmin, tipo)
        lab = (f'$\\{nome}$' if nome == 'tau' else '$\\eta$') \
            + f', $J={Jv:+.1f}$: ' \
            + (f'$\\Delta\\phi={dphi:+.3f}$' if tipo == 'rimbalzo'
               else f'spiral, $\\Delta\\phi={dphi:+.3f}$')
        axa.plot(rr_ * np.cos(ff_), rr_ * np.sin(ff_), color=col, ls=ls,
                 lw=1.3, label=lab)
th = np.linspace(0, 2 * np.pi, 200)
axa.plot(2 * np.cos(th), 2 * np.sin(th), 'k:', lw=0.8)
axa.plot((1 + np.sqrt(1 - s_n**2)) * np.cos(th),
         (1 + np.sqrt(1 - s_n**2)) * np.sin(th), 'k-', lw=0.8)
axa.set_aspect('equal')
axa.set_xlabel('$x$')
axa.set_ylabel('$y$')
axa.set_title(f'Kerr ($A=1$, $s={s_n}$, $E={E_n}$, arrival $r={r1}$, '
              f'$|J|={Jabs}$):\n$\\tau$ BOUNCES (mirror pair, no '
              'wind); $\\eta$ SPIRALS onto the horizon')
axa.legend(fontsize=6, loc='lower left')

asym_tau = abs(res[('tau', Jabs)][1] - res[('tau', -Jabs)][1])
asym_eta = res[('eta', Jabs)][0] + res[('eta', -Jabs)][0]
print(f"  specularita' tau: |r_min(+J)-r_min(-J)| = {asym_tau:.2e}")
print(f"  asimmetria eta (vento): Dphi(+J)+Dphi(-J) = {asym_eta:+.4f}")

print('pannello (b): scalari conformi...')
rg = np.logspace(np.log10(2.05), np.log10(60), 600)
fg = 1 - 2 * M_n / rg
for Av, col in ((1.0, 'C0'), (1.6, 'C2'), (2.2, 'C3')):
    vb = 1 - Av**2 * fg / E_n**2
    okm = vb > 1e-12
    n_ = np.where(okm, E_n / np.sqrt(np.abs(E_n**2 - Av**2 * fg)), np.nan)
    kt = np.where(okm, Av**2 * fg / (E_n * np.sqrt(np.abs(vb))), np.nan)
    axb.plot(rg[okm], n_[okm], color=col, lw=1.4, label=f'$n$, $A={Av}$')
    axb.plot(rg[okm], kt[okm], color=col, lw=1.4, ls='--',
             label=f'$k_\\tau$, $A={Av}$')
    if Av > E_n:
        rfz = 2 * M_n / (1 - E_n**2 / Av**2)
        axb.axvline(rfz, color=col, lw=0.7, ls=':')
axb.set_xscale('log')
axb.set_yscale('log')
axb.set_xlabel('$r$')
axb.set_ylabel('conformal scalars')
axb.set_title('the two scalars breathing on the rigid geometry '
              '$\\alpha_K$:\n$n=\\hat E/\\sqrt{\\hat E^2-A^2f}$ '
              '($\\eta,t$ branches), $k_\\tau=A^2f/(\\hat E\\bar v)$ '
              '($\\tau$ branch)')
axb.legend(fontsize=6, ncol=2)
savefig(fig, OUT, 'fig_thakurta_kerr_rami')
print('FATTO.')
