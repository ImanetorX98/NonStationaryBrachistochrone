# -*- coding: utf-8 -*-
"""
Validation figures — FLRW brachistochrones (de Sitter, H=1, a=-1/eta).
Single-column layout, English labels. Output in FLRWfigures/ (pdf+png):

  fig_flrw_cinematica   v(t) and E_phys(t)=Ehat/a for several Ehat; freezing a=Ehat
  fig_flrw_worldlines   x(eta) vs comoving Hubble horizon; analytic Dx_max
                        vs numerical quadrature (validation R5)
  fig_flrw_variazionale T_t(eps) and T_tau(eps) over a family of perturbed paths:
                        minimum at eps=0 for BOTH branches (validation R2)
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'FLRWfigures')

H = 1.0
EHATS = [1.5, 2.0, 3.0, 5.0]
COLS = plt.cm.viridis(np.linspace(0.15, 0.85, len(EHATS)))

def dx_max_analitica(Eh):
    return (np.sqrt(Eh**2 - 1) + np.arcsin(1 / Eh) - np.pi / 2) / (Eh * H)

# ---------------------------------------------------------------- fig 1
# two stacked panels (column width)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COL, 4.4))
t = np.linspace(0, np.log(max(EHATS)) + 0.3, 600)
a = np.exp(H * t)
for Eh, c in zip(EHATS, COLS):
    v = np.sqrt(np.clip(1 - a**2 / Eh**2, 0, None))
    v[a > Eh] = np.nan
    ax1.plot(t, v, color=c, label=rf'$\hat E={Eh}$')
    ax1.axvline(np.log(Eh) / H, color=c, ls=':', lw=0.8)
    Eph = Eh / a
    Eph[a > Eh] = np.nan
    ax2.plot(t, Eph, color=c)
    ax2.axvline(np.log(Eh) / H, color=c, ls=':', lw=0.8)
ax1.set_ylabel(r'local speed $v$')
ax1.set_title(r'freezing: $v\to 0$ at $a=\hat E$')
ax1.legend(ncol=2)
ax2.axhline(1, color='k', lw=0.6, ls='--')
ax2.set_xlabel(r'cosmic time $t\,H$')
ax2.set_ylabel(r'$E_{\rm phys}=\gamma=\hat E/a$')
ax2.set_title(r'rail redshift ($E_{\rm phys}\propto 1/a$)')
fig.suptitle(r'de Sitter FLRW: conformal-rail kinematics $-u_\eta=\hat E$')
savefig(fig, OUT, 'fig_flrw_cinematica')

# ---------------------------------------------------------------- fig 2
fig, ax = plt.subplots(figsize=(COL, COL * 0.85))
eta0 = -1.0 / H
err_max = 0.0
for Eh, c in zip(EHATS, COLS):
    eta_f = -1.0 / (Eh * H)
    eta = np.linspace(eta0, eta_f, 4000)
    aa = -1.0 / (H * eta)
    v = np.sqrt(np.clip(1 - aa**2 / Eh**2, 0, None))
    x = np.concatenate([[0], np.cumsum((v[1:] + v[:-1]) / 2 * np.diff(eta))])
    ax.plot(x, eta, color=c, label=rf'$\hat E={Eh}$')
    dx_num, dx_an = x[-1], dx_max_analitica(Eh)
    err_max = max(err_max, abs(dx_num - dx_an))
    ax.plot(dx_an, eta_f, 'o', color=c, ms=5, mfc='none')
    ax.plot(dx_num, eta_f, '.', color=c, ms=4)
eta = np.linspace(eta0, -1e-3, 200)
ax.plot(-eta, eta, 'k--', lw=1.2, label=r'Hubble horizon $x=-\eta$')
ax.plot(eta - eta0, eta, 'r-', lw=1.0, alpha=0.7,
        label=r'null ray ($\hat E\to\infty$)')
ax.set_xlabel(r'comoving distance $x$')
ax.set_ylabel(r'conformal time $\eta$')
ax.set_title(r'comoving worldlines: freezing below the horizon'
             '\n'
             rf'(circles = analytic $\Delta x_{{\max}}$; max quad. err.'
             rf' $= {err_max:.1e}$)')
ax.legend(loc='upper right', ncol=2)
ax.set_xlim(0, 1.05)
ax.set_ylim(eta0 * 1.02, 0)
savefig(fig, OUT, 'fig_flrw_worldlines')
print(f'  [validation] |Dx_num - Dx_analytic| max = {err_max:.3e}')

# ---------------------------------------------------------------- fig 3
Eh, L = 3.0, 0.3
eta_f = -1.0 / (Eh * H)
eta = np.linspace(eta0, eta_f - 1e-12, 60000)
aa = -1.0 / (H * eta)
v = np.sqrt(np.clip(1 - aa**2 / Eh**2, 0, None))
cum = np.concatenate([[0], np.cumsum((v[1:] + v[:-1]) / 2 * np.diff(eta))])

def arrivo(ell):
    return np.interp(ell, cum, eta)

sg = np.linspace(0, 1, 2000)
eps_grid = np.linspace(-0.16, 0.16, 81)
Tt, Ttau = [], []
for eps in eps_grid:
    dl = np.sqrt(L**2 + (eps * np.pi * np.cos(np.pi * sg))**2)
    ell = np.trapezoid(dl, sg)
    eta1 = arrivo(ell)
    a1 = -1.0 / (H * eta1)
    Tt.append(np.log(a1) / H)
    Ttau.append((a1 - 1) / (Eh * H))
Tt, Ttau = np.array(Tt), np.array(Ttau)

fig, ax = plt.subplots(figsize=(COL, COL * 0.78))
ax.plot(eps_grid, Tt / Tt[len(Tt) // 2], 'C0-',
        label=r'$T_t(\varepsilon)/T_t(0)$')
ax.plot(eps_grid, Ttau / Ttau[len(Ttau) // 2], 'C3--',
        label=r'$T_\tau(\varepsilon)/T_\tau(0)$')
ax.axvline(0, color='k', lw=0.6)
ax.set_xlabel(r'perturbation amplitude $\varepsilon$ '
              r'($y=\varepsilon\sin\pi s$)')
ax.set_ylabel(r'travel time (normalized)')
ax.set_title(r'variational check R2: the comoving straight line'
             '\n'
             rf'minimizes BOTH functionals ($\hat E={Eh}$, $L={L}$)')
ax.legend()
savefig(fig, OUT, 'fig_flrw_variazionale')
print(f'  [validation] argmin T_t   at eps = {eps_grid[np.argmin(Tt)]:.4f}')
print(f'  [validation] argmin T_tau at eps = {eps_grid[np.argmin(Ttau)]:.4f}')
print('DONE.')
