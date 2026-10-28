# -*- coding: utf-8 -*-
"""
Colormap Delta_r(A, J) for t/tau brachistochrones with FIXED symmetric ENDPOINTS
(Schwarzschild, conformal factor through E_eff = Ehat/A). J axis = momentum of the
tau branch (it fixes the opening angle Phi of the endpoints (r0, +/-Phi)).

For each (E_eff, J_tau):
  - tau branch with J_tau: r_min^tau and Phi (half-opening);
  - t branch with the SAME endpoints (r0, +/-Phi): shoot J_t -> r_min^t;
  - Delta_r = r_min^t - r_min^tau.

Theorem: n_t/n_tau = E_eff/f > 1 for every r, E => Delta_r > 0 ALWAYS
(no conformal inversion with fixed endpoints). The colormap confirms it.
"""

import os
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
# paper_style.py lives in NonStationaryMetrics/; a fresh clone must find it
# without relying on an unarchived sibling directory
sys.path.insert(0, os.path.join(_ROOT, "NonStationaryMetrics"))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x, **k):
        return x

import warnings
# NaN-producing branches of the lambdified radicals are expected and filtered by
# np.isfinite below; silence only these messages so they do not break the bar.
warnings.filterwarnings('ignore', message='invalid value encountered', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='divide by zero encountered', category=RuntimeWarning)

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
    """shoot J_t: dphi_half(t) = Phi. Returns r_min^t (shallow solution)."""
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

# grid (J_tau, E_eff)  ->  A = Ehat/E_eff
Jg = np.linspace(0.5, 4.5, 44)
Eg = np.linspace(1.05, 2.4, 40)
Z = np.full((len(Eg), len(Jg)), np.nan)
# one progress step per cell (same row-major order as before), with an ETA
cells = [(i, j) for i in range(len(Eg)) for j in range(len(Jg))]
for i, j in tqdm(cells, desc='cells', ncols=70, mininterval=1.0):
    Z[i, j] = Dr(Eg[i], Jg[j])

print(f"Delta_r in [{np.nanmin(Z):+.3f}, {np.nanmax(Z):+.3f}]", flush=True)
print(">>> " + ("INVERSION (sign change)" if np.nanmin(Z) * np.nanmax(Z) < 0
                else "NO inversion: Delta_r > 0 everywhere (theorem)"),
      flush=True)

# --------------------------------------------------------------- figure
Ag = Eh / Eg                                   # A axis
fig, ax = plt.subplots(figsize=(COL, COL * 0.9))
im = ax.pcolormesh(Jg, Ag, Z, cmap='viridis', shading='auto')
plt.colorbar(im, ax=ax, label=r'$\Delta r = r_{\min}^{t}-r_{\min}^{\tau}$')
ax.set_xlabel(r'$J$ ($\tau$-branch momentum $\to$ endpoint angle)')
ax.set_ylabel(r'$A$ (conformal factor, $E_{eff}=\hat E/A$)')
ax.set_title('Fixed symmetric endpoints: $\\Delta r$ in $(J,A)$\n'
             '$\\Delta r>0$ on the sampled grid (no conformal inversion there;\n'
             '$n_t/n_\\tau=E/f>1$)')
savefig(fig, HERE, 'fig_colormap_conformal_fixed_endpoints')
print("DONE.", flush=True)
