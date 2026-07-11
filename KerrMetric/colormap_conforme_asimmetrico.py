# -*- coding: utf-8 -*-
"""
Colormap Delta_r(A, J) per brachistocrone t/tau a ESTREMI FISSI
ASIMMETRICI (Schwarzschild, fattore conforme via E_eff = Ehat/A).
Estremi A=(rA,0), B=(rB, Delta_phi) con rA != rB. L'asse J = momento del
ramo tau, che (dati rA,rB) fissa l'apertura angolare Delta_phi.

Per ogni (E_eff, J_tau):
  - ramo tau con J_tau: periasse r_min^tau e Delta_phi = ang(rmin,rA)
    + ang(rmin,rB) (arco che scende da rA al periasse e risale a rB);
  - ramo t agli STESSI estremi: spara J_t -> Delta_phi_t = Delta_phi;
  - Delta_r = r_min^t - r_min^tau.

Teorema (indipendente dagli estremi): n_t/n_tau = E_eff/f > 1 => il ramo
t svolta piu' in alto ovunque => Delta_r > 0 SEMPRE. La colormap conferma
che il verdetto non dipende dalla simmetria degli estremi.
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
M, rA, rB, Eh = 1.0, 10.0, 6.0, 1.2        # estremi asimmetrici

def f(r):
    return 1 - 2 * M / r

def n_tau(r, E):
    return np.sqrt(f(r) / (E**2 - f(r)))

def n_t(r, E):
    return E / np.sqrt(f(r) * (E**2 - f(r)))

def rmin_of_J(N, E, J):
    """periasse: radice di N(r) r = J con r < rB (deve raggiungere rB)."""
    rs = np.linspace(rB * (1 - 1e-7), 2.0001 * M, 5000)
    g = N(rs, E) * rs - J
    for i in range(len(rs) - 1):
        if g[i] * g[i + 1] < 0:
            return brentq(lambda r: N(r, E) * r - J, rs[i + 1], rs[i])
    return None

def ang(N, E, J, rm, r_hi):
    integ = lambda r: J / (r * np.sqrt(f(r))
                           * np.sqrt(N(r, E)**2 * r**2 - J**2))
    return quad(integ, rm + 1e-10, r_hi, limit=150)[0]

def dphi_tot(N, E, J, rm):
    return ang(N, E, J, rm, rA) + ang(N, E, J, rm, rB)

def t_rmin_at_Phi(E, Phi):
    """spara J_t: dphi_tot(t) = Phi. Ritorna r_min^t (shallow)."""
    Jmax = n_t(rB, E) * rB * (1 - 1e-7)        # periasse < rB
    Js = np.linspace(0.05 * Jmax, Jmax, 220)
    sols, prev = [], None
    for J in Js:
        rm = rmin_of_J(n_t, E, J)
        if rm is None:
            prev = None
            continue
        d = dphi_tot(n_t, E, J, rm) - Phi
        if prev is not None and prev[1] * d < 0:
            try:
                Js = brentq(lambda x: dphi_tot(n_t, E, x,
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
    Phi = dphi_tot(n_tau, E, J_tau, rm_tau)
    rm_t = t_rmin_at_Phi(E, Phi)
    return rm_t - rm_tau if np.isfinite(rm_t) else np.nan

# griglia (J_tau, E_eff)  ->  A = Ehat/E_eff
Jg = np.linspace(0.5, 5.0, 44)
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
Ag = Eh / Eg
fig, ax = plt.subplots(figsize=(COL, COL * 0.9))
im = ax.pcolormesh(Jg, Ag, Z, cmap='viridis', shading='auto')
plt.colorbar(im, ax=ax, label=r'$\Delta r = r_{\min}^{t}-r_{\min}^{\tau}$')
if np.nanmin(Z) * np.nanmax(Z) < 0:
    try:
        ax.contour(Jg, Ag, Z, levels=[0], colors='r', linewidths=1.4)
    except Exception:
        pass
ax.set_xlabel(r'$J$ ($\tau$-branch momentum $\to$ endpoint angle)')
ax.set_ylabel(r'$A$ (conformal factor, $E_{eff}=\hat E/A$)')
ax.set_title('Fixed ASYMMETRIC endpoints '
             f'($r_A={rA:.0f},r_B={rB:.0f}$): $\\Delta r$ in $(J,A)$\n'
             '$>0$ everywhere (no conformal inversion; $n_t/n_\\tau=E/f>1$)')
savefig(fig, HERE, 'fig_colormap_conforme_asimmetrico')
print("FATTO.", flush=True)
