# -*- coding: utf-8 -*-
"""
Colormap Delta_r(A, J) per brachistocrone t/tau a ESTREMI FISSI simmetrici
(Schwarzschild, fattore conforme via E_eff = Ehat/A). Asse J = momento del
ramo tau (fissa l'apertura angolare Phi degli estremi (r0, +/-Phi)).

Per ogni (E_eff, J_tau):
  - ramo tau con J_tau: r_min^tau e Phi (semi-apertura);
  - ramo t agli STESSI estremi (r0, +/-Phi): spara J_t -> r_min^t;
  - Delta_r = r_min^t - r_min^tau.

Teorema: n_t/n_tau = E_eff/f > 1 per ogni r,E => Delta_r > 0 SEMPRE
(nessuna inversione conforme a estremi fissi). La colormap conferma.
"""

import os
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x, **k):
        return x

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, r0, Eh = 1.0, 8.0, 1.2

def f(r):
    return 1 - 2 * M / r

def n_tau(r, E):
    return np.sqrt(f(r) / (E**2 - f(r)))

def n_t(r, E):
    return E / np.sqrt(f(r) * (E**2 - f(r)))

def rmin_of_J(N, E, J):
    rs = np.linspace(r0 * (1 - 1e-7), 2.0001 * M, 4000)
    g = N(rs, E) * rs - J
    for i in range(len(rs) - 1):
        if g[i] * g[i + 1] < 0:
            return brentq(lambda r: N(r, E) * r - J, rs[i + 1], rs[i])
    return None

def dphi_half(N, E, J, rm):
    integ = lambda r: J / (r * np.sqrt(f(r))
                           * np.sqrt(N(r, E)**2 * r**2 - J**2))
    return quad(integ, rm + 1e-10, r0, limit=150)[0]

def t_rmin_at_Phi(E, Phi):
    """spara J_t: dphi_half(t) = Phi. Ritorna r_min^t (shallow)."""
    Jmax = n_t(r0, E) * r0 * (1 - 1e-7)
    Js = np.linspace(0.05 * Jmax, Jmax, 200)
    sols, prev = [], None
    for J in Js:
        rm = rmin_of_J(n_t, E, J)
        if rm is None:
            prev = None
            continue
        d = dphi_half(n_t, E, J, rm) - Phi
        if prev is not None and prev[1] * d < 0:
            try:
                Js = brentq(lambda x: dphi_half(n_t, E, x,
                            rmin_of_J(n_t, E, x)) - Phi, prev[0], J)
                sols.append(rmin_of_J(n_t, E, Js))
            except (ValueError, TypeError):
                pass
        prev = (J, d)
    return max(sols) if sols else np.nan

def Dr(E, J_tau):
    rm_tau = rmin_of_J(n_tau, E, J_tau)
    if rm_tau is None:
        return np.nan
    Phi = dphi_half(n_tau, E, J_tau, rm_tau)
    rm_t = t_rmin_at_Phi(E, Phi)
    return rm_t - rm_tau if np.isfinite(rm_t) else np.nan

# griglia (J_tau, E_eff)  ->  A = Ehat/E_eff
Jg = np.linspace(0.5, 4.5, 44)
Eg = np.linspace(1.05, 2.4, 40)
Z = np.full((len(Eg), len(Jg)), np.nan)
for i, E in enumerate(tqdm(Eg, desc='E_eff (=Eh/A)', ncols=70)):
    for j, J in enumerate(Jg):
        Z[i, j] = Dr(E, J)

print(f"Delta_r in [{np.nanmin(Z):+.3f}, {np.nanmax(Z):+.3f}]", flush=True)
print(">>> " + ("INVERSIONE (cambia segno)" if np.nanmin(Z) * np.nanmax(Z) < 0
                else "NESSUNA inversione: Delta_r > 0 ovunque (teorema)"),
      flush=True)

# --------------------------------------------------------------- figura
Ag = Eh / Eg                                   # asse A
fig, ax = plt.subplots(figsize=(COL, COL * 0.9))
im = ax.pcolormesh(Jg, Ag, Z, cmap='viridis', shading='auto')
plt.colorbar(im, ax=ax, label=r'$\Delta r = r_{\min}^{t}-r_{\min}^{\tau}$')
ax.set_xlabel(r'$J$ ($\tau$-branch momentum $\to$ endpoint angle)')
ax.set_ylabel(r'$A$ (conformal factor, $E_{eff}=\hat E/A$)')
ax.set_title('Fixed symmetric endpoints: $\\Delta r$ in $(J,A)$\n'
             '$>0$ everywhere (no conformal inversion; $n_t/n_\\tau=E/f>1$)')
savefig(fig, HERE, 'fig_colormap_conforme_estremi_fissi')
print("FATTO.", flush=True)
