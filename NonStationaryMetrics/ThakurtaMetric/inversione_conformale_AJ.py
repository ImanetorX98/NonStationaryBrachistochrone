# -*- coding: utf-8 -*-
"""
Inversione di plunge conformale in Thakurta-Kerr equatoriale, piano (J, A).
Stavolta col RANGE DI J GIUSTO: l'inversione vive a J ~ 16-24, fuori dal
vecchio scan [3.9, 5.5]. Mostra il rosso (tau piu' profondo).

r_min di ciascun ramo = radice esterna della condizione di svolta
(p_r = 0, H = 0), calcolata direttamente (nessuna integrazione):
  tau:  H_tau(r,0,J) = ptp*php0 + sqrt(R2)|ptp|/sqrt(Pb) - A^2 f/Ehat = 0
  t  :  H_t  (r,0,J) = J*php0   + sqrt(R2)|J|/sqrt(Pb)   - 1        = 0
Curva analitica di inversione (r_min^tau = r_min^t):
  cond(r,A) = A^2 b Q - Pb (Ehat - A^2 f) = 0 ,  Q = b vb2 + sqrt(vb2 Dl)
con  cond(A=1) < 0,  cond(A_freeze) = Ehat r Dl (Ehat-1)/(r-2M) > 0
=> A_inv in (1, A_freeze): l'inversione esiste sempre ed e' raggiungibile.
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
OUT = os.path.join(HERE, 'Thakurtafigures')
M, a, Eh = 1.0, 0.9, 1.2
r_e = 2 * M

# ---- condizioni di svolta (numeriche) --------------------------------
def geom(r, A):
    f = 1 - 2 * M / r
    Dl = r**2 - 2 * M * r + a**2
    P = r**2 + a**2 + 2 * M * a**2 / r
    b = 2 * M * a / r
    vb2 = 1 - A**2 * f / Eh**2
    Pb = P + A**2 * b**2 / Eh**2
    return f, Dl, P, b, vb2, Pb

def Htau0(r, A, J):
    f, Dl, P, b, vb2, Pb = geom(r, A)
    if vb2 <= 0:
        return np.nan
    ptp = J - A**2 * b / Eh
    php0 = b * vb2 / Pb
    S = np.sqrt(vb2 * Dl) / Pb
    return ptp * php0 + S * abs(ptp) - A**2 * f / Eh

def Ht0(r, A, J):
    f, Dl, P, b, vb2, Pb = geom(r, A)
    if vb2 <= 0:
        return np.nan
    php0 = b * vb2 / Pb
    S = np.sqrt(vb2 * Dl) / Pb
    return J * php0 + S * abs(J) - 1.0

def r_outer(func, A, J):
    rs = np.linspace(20.0, r_e + 1e-4, 3000)
    vals = np.array([func(rr, A, J) for rr in rs])
    for i in range(len(rs) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) \
                and vals[i] * vals[i + 1] < 0:
            return brentq(lambda rr: func(rr, A, J), rs[i + 1], rs[i])
    return np.nan

# ---- griglia (J, A) --------------------------------------------------
J_g = np.linspace(5.0, 30.0, 60)
A_g = np.linspace(1.0, 1.5, 55)
Dr = np.full((len(A_g), len(J_g)), np.nan)
for i, Av in enumerate(A_g):
    for j, Jv in enumerate(J_g):
        rtau = r_outer(Htau0, Av, Jv)
        rt = r_outer(Ht0, Av, Jv)
        if np.isfinite(rtau) and np.isfinite(rt):
            Dr[i, j] = rt - rtau
n_red = int(np.nansum(Dr > 0))
print(f"celle con INVERSIONE (Dr>0, tau piu' profondo): {n_red}/"
      f"{np.sum(np.isfinite(Dr))}")
print(f"  Dr range: min {np.nanmin(Dr):+.3f}, max {np.nanmax(Dr):+.3f}")

# ---- curva analitica cond=0: J_inv(A) --------------------------------
r_s, A_s = sp.symbols('r A', positive=True)
f_s = 1 - 2 * M / r_s
Dl_s = r_s**2 - 2 * M * r_s + a**2
P_s = r_s**2 + a**2 + 2 * M * a**2 / r_s
b_s = 2 * M * a / r_s
vb2_s = 1 - A_s**2 * f_s / Eh**2
Pb_s = P_s + A_s**2 * b_s**2 / Eh**2
Q_s = b_s * vb2_s + sp.sqrt(vb2_s * Dl_s)
cond = A_s**2 * b_s * Q_s - Pb_s * (Eh - A_s**2 * f_s)
condf = sp.lambdify((r_s, A_s), cond, 'numpy')
Jinvf = sp.lambdify((r_s, A_s), Pb_s / Q_s, 'numpy')
A_curve, J_curve = [], []
for Av in np.linspace(A_g[0] + 0.02, A_g[-1], 120):
    Afr = Eh / np.sqrt(1 - 2 * M / 3.0)   # dummy large; trova radice r*>2M
    rs = np.linspace(2.05, 40.0, 8000)
    with np.errstate(invalid='ignore'):
        cv = np.array([condf(rr, Av) for rr in rs])
    rstar = None
    for k in range(len(rs) - 1):
        if np.isfinite(cv[k]) and np.isfinite(cv[k + 1]) \
                and cv[k] * cv[k + 1] < 0:
            rstar = brentq(lambda rr: condf(rr, Av), rs[k], rs[k + 1])
            break
    if rstar is not None:
        A_curve.append(Av)
        J_curve.append(float(Jinvf(rstar, Av)))
A_curve, J_curve = np.array(A_curve), np.array(J_curve)

# --------------------------------------------------------------- figura
fig, ax = plt.subplots(figsize=(COL, COL * 0.92))
# scala colore clippata perche' il blu profondo (Dr~-5) slaverebbe il rosso
vmax = min(1.5, np.nanmax(Dr) * 4) if np.nanmax(Dr) > 0 else 1.5
vmax = max(vmax, 1.0)
im = ax.pcolormesh(J_g, A_g, Dr, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                   shading='auto')
cb = plt.colorbar(im, ax=ax, extend='min',
                  label='$\\Delta r = r_{\\min}^{t}-r_{\\min}^{\\tau}$')
CS = ax.contour(J_g, A_g, Dr, levels=[0.0], colors='k', linewidths=1.6)
msk = J_curve <= J_g[-1]
ax.plot(J_curve[msk], A_curve[msk], 'w--', lw=1.6,
        label='analytic $\\mathrm{cond}(r,A)=0$')
ax.text(9, 1.05, 'blue: $t$ deeper', fontsize=6.5, color='magenta',
        fontweight='bold')
ax.text(24, 1.42, 'red:\n$\\tau$ deeper', fontsize=6.5, color='darkred',
        ha='center', fontweight='bold')
ax.set_xlabel('$J$ (particle angular momentum)')
ax.set_ylabel('$A$ (conformal factor)')
ax.set_xlim(J_g[0], J_g[-1])
ax.set_ylim(A_g[0], A_g[-1])
ax.set_title('Conformal plunge inversion (Thakurta-Kerr, $s=0.9$,\n'
             f'$\\hat E=1.2$): $\\Delta r$ changes sign; inversion curve\n'
             'analytic $=$ numeric, REACHABLE ($A_{inv}<A_{freeze}$)')
ax.legend(loc='lower right', fontsize=6)
savefig(fig, OUT, 'fig_thakurta_kerr_inversione_AJ')
print('FATTO.')
