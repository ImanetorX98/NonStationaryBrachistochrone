# -*- coding: utf-8 -*-
"""
Testimonianza: esiste un J di penetrazione (cattura all'ergosfera r_e=2M)
anche per J NEGATIVO, sia per il ramo t che per il ramo tau di
Thakurta-Kerr equatoriale.

(a) r_min raggiunto vs J (entrambi i rami): plateau a r_e nella banda di
    cattura, con soglie negative J_neg^tau = -sA^2/Ehat e J_neg^t ~ -8.
(b) traiettorie equatoriali reali a J<0: catturate (toccano r_e) vs
    diffuse, per ciascun ramo.
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
OUT = os.path.join(HERE, 'Thakurtafigures')
M, a, Eh, A_val, r0 = 1.0, 0.9, 1.2, 1.0, 8.0
r_e = 2 * M
r_plus = M + np.sqrt(M**2 - a**2)

r, pr, J = sp.symbols('r p_r J', real=True)
f = 1 - 2 * M / r
Dl = r**2 - 2 * M * r + a**2
P = r**2 + a**2 + 2 * M * a**2 / r
b = 2 * M * a / r
vb2 = 1 - A_val**2 * f / Eh**2
Pb = P + A_val**2 * b**2 / Eh**2
php0 = b * vb2 / Pb
R2 = vb2 * Dl / Pb

def build(branch):
    if branch == 'tau':
        ptp = J - A_val**2 * b / Eh
        H = ptp * php0 + sp.sqrt(R2) * sp.sqrt((Dl / r**2) * pr**2
                                               + ptp**2 / Pb) - A_val**2 * f / Eh
    else:
        H = J * php0 + sp.sqrt(R2) * sp.sqrt((Dl / r**2) * pr**2
                                             + J**2 / Pb) - 1
    return (sp.lambdify((r, pr, J), H, 'numpy'),
            sp.lambdify((r, pr, J), sp.diff(H, pr), 'numpy'),
            sp.lambdify((r, pr, J), sp.diff(H, r), 'numpy'),
            sp.lambdify((r, pr, J), sp.diff(H, J), 'numpy'))

def traiettoria(fns, J_val):
    """Ritorna (r_arr, phi_arr, r_min, catturata)."""
    Hf, dHp, dHr, dHJ = fns
    pg = np.linspace(-400, 400, 200001)
    with np.errstate(invalid='ignore'):
        Hv = Hf(r0, pg, J_val)
    roots = [brentq(lambda p: Hf(r0, p, J_val), pg[i], pg[i + 1])
             for i in range(len(pg) - 1)
             if np.isfinite(Hv[i]) and np.isfinite(Hv[i + 1])
             and Hv[i] * Hv[i + 1] < 0]
    ing = [p for p in roots if dHp(r0, p, J_val) < 0]
    if not ing:
        return None
    p0 = min(ing)
    ev_turn = lambda t_, y: y[1]
    ev_turn.terminal, ev_turn.direction = True, 1
    ev_erg = lambda t_, y: y[0] - (r_e - 1e-4)
    ev_erg.terminal, ev_erg.direction = True, -1
    s = solve_ivp(lambda t_, y: [dHp(y[0], y[1], J_val),
                                 -dHr(y[0], y[1], J_val),
                                 dHJ(y[0], y[1], J_val)],
                  [0, 5000], [r0, p0, 0.0], rtol=1e-11, atol=1e-13,
                  method='DOP853', events=[ev_turn, ev_erg],
                  dense_output=True, max_step=0.3)
    tt = np.linspace(0, s.t[-1], 600)
    Y = s.sol(tt)
    rmin = Y[0].min()
    return Y[0], Y[2], rmin, (rmin <= r_e + 1e-3)

fns = {'tau': build('tau'), 't': build('t')}

# (a) r_min vs J
print("r_min(J)...")
Jg_tau = np.linspace(-1.6, 1.6, 121)
Jg_t = np.linspace(-11.0, 3.0, 141)
rmin_tau = [(traiettoria(fns['tau'], Jv) or (None, None, np.nan, None))[2]
            for Jv in Jg_tau]
rmin_t = [(traiettoria(fns['t'], Jv) or (None, None, np.nan, None))[2]
          for Jv in Jg_t]
Jc = a * A_val**2 / Eh

fig, (a1, a2) = plt.subplots(2, 1, figsize=(COL, 6.2))
a1.plot(Jg_t, rmin_t, 'C0-', label='$t$ (η) branch')
a1.plot(Jg_tau, rmin_tau, 'C3-', label=r'$\tau$ branch')
a1.axhline(r_e, color='k', ls='--', lw=0.9, label='ergosphere $r_e=2M$')
a1.axvline(-Jc, color='C3', ls=':', lw=0.9)
a1.axvline(-8.05, color='C0', ls=':', lw=0.9)
a1.annotate(f'$J_{{neg}}^{{\\tau}}=-sA^2/\\hat E={-Jc:.2f}$', (-Jc, 4.2),
            fontsize=6, color='C3',
            xytext=(-4.5, 5.0),
            arrowprops=dict(arrowstyle='->', lw=0.5, color='C3'))
a1.annotate(r'$J_{neg}^{t}\approx-8.05$', (-8.05, 2.0), fontsize=6,
            color='C0', xytext=(-10.5, 4.0),
            arrowprops=dict(arrowstyle='->', lw=0.5, color='C0'))
a1.set_xlabel('$J$ (particle angular momentum; $J<0$ = retrograde)')
a1.set_ylabel('$r_{\\min}$ reached')
a1.set_ylim(1.8, 8.2)
a1.set_title('Capture reaches $r_e$ for a BAND of $J$, incl. negative:\n'
             r'$\tau$ symmetric $[-sA^2/\hat E,+sA^2/\hat E]$; $t$ wide')
a1.legend(fontsize=6, loc='upper center', ncol=2)

# (b) traiettorie a J<0
print("traiettorie J<0...")
casi = [('tau', -0.5, 'C3', '-', r'$\tau$, $J=-0.5$: capture'),
        ('tau', -1.3, 'C3', ':', r'$\tau$, $J=-1.3$: scatter'),
        ('t', -5.0, 'C0', '-', r'$t$, $J=-5$: capture'),
        ('t', -9.5, 'C0', ':', r'$t$, $J=-9.5$: scatter')]
for branch, Jv, col, ls, lab in casi:
    out = traiettoria(fns[branch], Jv)
    if out is None:
        continue
    rA, phiA, rmin, cap = out
    a2.plot(rA * np.cos(phiA), rA * np.sin(phiA), color=col, ls=ls,
            lw=1.6 if cap else 1.1, label=lab)
th = np.linspace(0, 2 * np.pi, 200)
a2.plot(r_e * np.cos(th), r_e * np.sin(th), 'k--', lw=1.0,
        label='ergosphere')
a2.plot(r_plus * np.cos(th), r_plus * np.sin(th), 'k-', lw=0.8)
a2.plot([r0], [0], 'k^', ms=6)
a2.set_aspect('equal')
a2.set_xlabel('$x$')
a2.set_ylabel('$y$')
a2.set_title('retrograde ($J<0$) orbits: some REACH $r_e$ (capture)')
a2.legend(fontsize=5.5, loc='lower left')
savefig(fig, OUT, 'fig_thakurta_tricotomia_Jneg')
print("FATTO.")
