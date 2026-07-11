# -*- coding: utf-8 -*-
"""
fig_brachistocrone_confronto (PaperFigures/): le brachistocrone nei tre
spaziotempi del paper, una accanto all'altra.

  (a) FLRW de Sitter: worldline comoventi, congelamento sotto l'orizzonte
  (b) Vaidya: rimbalzi tau (statico vs accresce vs evapora), stesso arrivo
  (c) Thakurta-Kerr (A=1): tau speculare (+-J) vs eta che spirala
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))

fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))

# ------------------------------------------------------------- (a) FLRW
ax = axs[0]
H = 1.0
for Eh, c in zip((1.5, 2.0, 3.0, 5.0),
                 plt.cm.viridis(np.linspace(0.15, 0.85, 4))):
    eta = np.linspace(-1.0, -1.0 / Eh, 3000)
    aa = -1.0 / eta
    v = np.sqrt(np.clip(1 - aa**2 / Eh**2, 0, None))
    x = np.concatenate([[0], np.cumsum((v[1:] + v[:-1]) / 2
                                       * np.diff(eta))])
    ax.plot(x, eta, color=c, label=rf'$\hat E={Eh}$')
    ax.plot(x[-1], eta[-1], 'o', color=c, ms=4)
eta = np.linspace(-1.0, -1e-3, 100)
ax.plot(-eta, eta, 'k--', lw=1.0, label='orizzonte Hubble')
ax.plot(eta + 1.0, eta, 'r-', lw=0.9, alpha=0.7, label='raggio nullo')
ax.set_xlim(0, 1.02)
ax.set_ylim(-1.01, 0)
ax.set_xlabel('$x$ comovente')
ax.set_ylabel(r'$\eta$')
ax.set_title('FLRW (de Sitter): rette comoventi,\ncongelamento a $a=\\hat E$')
ax.legend(fontsize=7)

# ----------------------------------------------------------- (b) Vaidya
ax = axs[1]
E_n, J_n, m0 = 1.2, 1.3, 1.0
r1, v_arr = 10.0, 40.0
vv, rr, pr = sp.symbols('vv rr p_r', real=True)
mu_s = sp.Symbol('mu', real=True)
f_v = 1 - 2 * (m0 + mu_s * vv) / rr
w_v = E_n**2 - f_v
H_tv = pr * (f_v - E_n**2) - E_n \
    + sp.sqrt(w_v) * sp.sqrt((E_n * pr + 1)**2 + J_n**2 / rr**2)
dHp = sp.lambdify((vv, rr, pr, mu_s), sp.diff(H_tv, pr), 'numpy')
dHr = sp.lambdify((vv, rr, pr, mu_s), sp.diff(H_tv, rr), 'numpy')
dHJ_expr = sp.diff(H_tv, sp.Symbol('J')) if False else None
Jsym = sp.Symbol('Js', positive=True)
H_tvJ = H_tv.subs(J_n, Jsym) if False else None
# dphi/dv esplicito: dH/dJ con J numerico
dHJ = sp.lambdify((vv, rr, pr, mu_s),
                  (sp.sqrt(w_v) * (J_n / rr**2)
                   / sp.sqrt((E_n * pr + 1)**2 + J_n**2 / rr**2)), 'numpy')
Hf = sp.lambdify((vv, rr, pr, mu_s), H_tv, 'numpy')

for mu_n, col, lab in ((0.0, 'k', 'statico'),
                       (0.01, 'C3', 'accresce $\\mu=+0.01$'),
                       (-0.01, 'C0', 'evapora $\\mu=-0.01$')):
    pg = np.linspace(-30, 30, 120001)
    with np.errstate(invalid='ignore'):
        Hg = Hf(v_arr, r1, pg, mu_n)
    roots = [brentq(lambda p: Hf(v_arr, r1, p, mu_n), pg[i], pg[i + 1])
             for i in range(len(pg) - 1)
             if np.isfinite(Hg[i]) and np.isfinite(Hg[i + 1])
             and Hg[i] * Hg[i + 1] <= 0]
    p0 = [p_ for p_ in roots if dHp(v_arr, r1, p_, mu_n) > 0][0]

    def rhs(v_, y):
        return [dHp(v_, y[0], y[1], mu_n), -dHr(v_, y[0], y[1], mu_n),
                dHJ(v_, y[0], y[1], mu_n)]
    ev_p = lambda v_, y: dHp(v_, y[0], y[1], mu_n)
    ev_p.terminal, ev_p.direction = True, 0
    s1 = solve_ivp(rhs, [v_arr, -200], [r1, p0, 0.0], rtol=1e-11,
                   atol=1e-13, method='DOP853', events=[ev_p],
                   dense_output=True)
    ev_r = lambda v_, y: y[0] - r1
    ev_r.terminal, ev_r.direction = True, 1
    s2 = solve_ivp(rhs, [s1.t_events[0][0], -400], s1.y_events[0][0],
                   rtol=1e-11, atol=1e-13, method='DOP853', events=[ev_r],
                   dense_output=True)
    tt1 = np.linspace(v_arr, s1.t[-1], 400)
    tt2 = np.linspace(s2.t[0], s2.t[-1], 400)
    rr_ = np.concatenate([s2.sol(tt2)[0][::-1], s1.sol(tt1)[0][::-1]])
    ff_ = np.concatenate([s2.sol(tt2)[2][::-1], s1.sol(tt1)[2][::-1]])
    ff_ = ff_ - ff_[0]
    ax.plot(rr_ * np.cos(ff_), rr_ * np.sin(ff_), color=col, lw=1.2,
            label=lab)
th = np.linspace(0, 2 * np.pi, 100)
ax.plot(2 * np.cos(th), 2 * np.sin(th), 'k:', lw=0.8)
ax.set_aspect('equal')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title('Vaidya: rimbalzo $\\tau$, stesso arrivo\n'
             '($E=1.2$, $J=1.3$): il timing conta')
ax.legend(fontsize=7, loc='lower right')

# ---------------------------------------------------- (c) Thakurta-Kerr
ax = axs[2]
s_n = 0.9
r_, pr_ = sp.symbols('r_th p_r_th', real=True)
f_k = 1 - 2 / r_
Dl = r_**2 - 2 * r_ + s_n**2
P_k = r_**2 + s_n**2 + 2 * s_n**2 / r_
vb2 = 1 - f_k / E_n**2
Pb = P_k + (2 * s_n / r_)**2 / E_n**2
php0 = (2 * s_n / r_) * vb2 / Pb
R_ = sp.sqrt(Dl * vb2 / Pb)
pphi = sp.Symbol('p_phi_c', real=True)
H_e = pphi * php0 + R_ * sp.sqrt((Dl / r_**2) * pr_**2 + pphi**2 / Pb) - 1
ptp = pphi - 2 * s_n / (r_ * E_n)
H_t = ptp * php0 + R_ * sp.sqrt((Dl / r_**2) * pr_**2 + ptp**2 / Pb) \
    - f_k / E_n
r1k = 8.0
r_plus = 1 + np.sqrt(1 - s_n**2)

def traccia(Hx, Jv):
    Hn = Hx.subs(pphi, Jv)
    dHp_ = sp.lambdify((r_, pr_), sp.diff(Hn, pr_), 'numpy')
    dHr_ = sp.lambdify((r_, pr_), sp.diff(Hn, r_), 'numpy')
    dHJ_ = sp.lambdify((r_, pr_), sp.diff(Hx, pphi).subs(pphi, Jv), 'numpy')
    Hf_ = sp.lambdify((r_, pr_), Hn, 'numpy')
    pg = np.linspace(-120, 120, 480001)
    with np.errstate(invalid='ignore'):
        Hg = Hf_(r1k, pg)
    roots = [brentq(lambda p: Hf_(r1k, p), pg[i], pg[i + 1])
             for i in range(len(pg) - 1)
             if np.isfinite(Hg[i]) and np.isfinite(Hg[i + 1])
             and Hg[i] * Hg[i + 1] <= 0]
    p0 = [p_ for p_ in roots if dHp_(r1k, p_) > 0][0]

    def rhs(e_, y):
        return [dHp_(y[0], y[1]), -dHr_(y[0], y[1]), dHJ_(y[0], y[1])]
    ev_p = lambda e_, y: dHp_(y[0], y[1])
    ev_p.terminal, ev_p.direction = True, 0
    ev_h = lambda e_, y: y[0] - r_plus * 1.02
    ev_h.terminal, ev_h.direction = True, -1
    s1 = solve_ivp(rhs, [0, -400], [r1k, p0, 0.0], rtol=1e-11, atol=1e-13,
                   method='DOP853', events=[ev_p, ev_h], dense_output=True)
    tt1 = np.linspace(0, s1.t[-1], 500)
    if len(s1.t_events[1]):
        rr_ = s1.sol(tt1)[0][::-1]
        ff_ = s1.sol(tt1)[2][::-1]
        return rr_, ff_ - ff_[0]
    ev_r = lambda e_, y: y[0] - r1k
    ev_r.terminal, ev_r.direction = True, 1
    s2 = solve_ivp(rhs, [s1.t_events[0][0], -800], s1.y_events[0][0],
                   rtol=1e-11, atol=1e-13, method='DOP853', events=[ev_r],
                   dense_output=True)
    tt2 = np.linspace(s2.t[0], s2.t[-1], 500)
    rr_ = np.concatenate([s2.sol(tt2)[0][::-1], s1.sol(tt1)[0][::-1]])
    ff_ = np.concatenate([s2.sol(tt2)[2][::-1], s1.sol(tt1)[2][::-1]])
    return rr_, ff_ - ff_[0]

for Hx, Jv, col, ls, lab in ((H_t, +1.3, 'C3', '-', '$\\tau$, $J=+1.3$'),
                             (H_t, -1.3, 'C1', '--', '$\\tau$, $J=-1.3$'),
                             (H_e, +1.3, 'C0', '-', '$\\eta$: spirala')):
    rr_, ff_ = traccia(Hx, Jv)
    ax.plot(rr_ * np.cos(ff_), rr_ * np.sin(ff_), color=col, ls=ls,
            lw=1.2, label=lab)
th = np.linspace(0, 2 * np.pi, 100)
ax.plot(2 * np.cos(th), 2 * np.sin(th), 'k:', lw=0.8)
ax.plot(r_plus * np.cos(th), r_plus * np.sin(th), 'k-', lw=0.8)
ax.set_aspect('equal')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title('Thakurta-Kerr ($A=1$): $\\tau$ speculare,\n'
             '$\\eta$ spirala (dragging + auto-sintonizzazione)')
ax.legend(fontsize=7, loc='lower right')

fig.suptitle('Brachistocrone a rotaia nei tre spaziotempi '
             '($-u\\cdot W = \\hat E$)', y=1.03)
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(HERE, f'fig_brachistocrone_confronto.{ext}'),
                dpi=200, bbox_inches='tight')
plt.close(fig)
print('salvata fig_brachistocrone_confronto')
