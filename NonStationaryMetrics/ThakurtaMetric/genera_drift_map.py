# -*- coding: utf-8 -*-
"""
Drift map (theta, r) della quasi-costante dinamica K_NLO(E_eff) in
Thakurta-Kerr — stile drift-colormap del progetto Kerr.

Drift puntuale totale lungo il flusso PMP 3D dinamico:

    D(r, theta) = | {K, H_3D} + eps a^2 S1 | / |K|

con momenti di riferimento fissati (p_r, p_th, J), A e eps = A'/A dati.
Colore: verde (poco drift) -> rosso (molto drift), r LOGARITMICO,
ergosfera tratteggiata r_e(th) = M + sqrt(M^2 - a^2 cos^2 th),
orizzonte pieno, regione oltre il muro di congelamento mascherata.
Output: Thakurtafigures/fig_thakurta_kerr_drift_map
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from paper_style import COL, set_style, savefig

set_style()
OUT = os.path.join(HERE, 'Thakurtafigures')

M_n, Eh_n = 1.0, 1.2
a_v, A_v, Hc_v = 0.9, 1.25, 0.02          # spin, fattore conforme, eps=Hc*A
J_v, pr_v, pth_v = 1.2, -0.3, 1.0
eps_v = Hc_v * A_v
E_eff = Eh_n / A_v
r_plus = M_n + np.sqrt(M_n**2 - a_v**2)

# ------------------------------------------------ H 3D e K_NLO simbolici
r, th = sp.symbols('r theta', positive=True)
pr, pth = sp.symbols('p_r p_theta', real=True)

Sig = r**2 + a_v**2 * sp.cos(th)**2
Dl = r**2 - 2 * M_n * r + a_v**2
f_S = 1 - 2 * M_n * r / Sig
b = 2 * M_n * a_v * r * sp.sin(th)**2 / Sig
G = (r**2 + a_v**2 + 2 * M_n * a_v**2 * r * sp.sin(th)**2 / Sig) \
    * sp.sin(th)**2
vbS2 = 1 - A_v**2 * f_S / Eh_n**2
Gb = G + A_v**2 * b**2 / Eh_n**2
php0 = b * vbS2 / Gb
R2 = vbS2 * Dl * sp.sin(th)**2 / Gb
ptphi = J_v - A_v**2 * b / Eh_n
H3 = ptphi * php0 + sp.sqrt(R2) * sp.sqrt((Dl / Sig) * pr**2
                                          + pth**2 / Sig + ptphi**2 / Gb) \
    - A_v**2 * f_S / Eh_n

D0 = r**2 - 2 * M_n * r
DE0 = (E_eff**2 - 1) * r**2 + 2 * M_n * r
N2 = (E_eff**2 - 1) * r**2 + 4 * M_n * r - 4 * M_n**2
f2 = N2 * sp.cos(2 * th) / (2 * r**2 * D0 * DE0)
K = (pth**2 + J_v**2 / sp.tan(th)**2 - a_v**2 * E_eff**2 * sp.cos(th)**2
     + a_v**2 * f2 * pth**2)

PB = (sp.diff(K, r) * sp.diff(H3, pr) + sp.diff(K, th) * sp.diff(H3, pth)
      - sp.diff(K, pth) * sp.diff(H3, th))
S1 = 2 * E_eff**2 * (sp.cos(th)**2
                     + M_n * sp.cos(2 * th) * pth**2 / (r * DE0**2))
drift = sp.Abs(PB + eps_v * a_v**2 * S1) / sp.Abs(K)

print("lambdify (cse)...")
Dr_f = sp.lambdify((r, th, pr, pth), drift, 'numpy', cse=True)
vb_f = sp.lambdify((r, th), vbS2, 'numpy')

# ------------------------------------------------------------- griglia
rg = np.logspace(np.log10(r_plus * 1.02), np.log10(30.0), 320)
tg = np.linspace(np.deg2rad(4), np.deg2rad(176), 320)
RR, TT = np.meshgrid(rg, tg)
with np.errstate(invalid='ignore', divide='ignore'):
    DD = Dr_f(RR, TT, pr_v, pth_v)
    VV = vb_f(RR, TT)
DD = np.where(VV > 0, DD, np.nan)          # oltre il muro: inaccessibile
DD = np.where(np.isfinite(DD), DD, np.nan)

print(f"drift: min {np.nanmin(DD):.2e}  mediana {np.nanmedian(DD):.2e}"
      f"  max {np.nanmax(DD):.2e}")

# --------------------------------------------------------------- figura
fig, ax = plt.subplots(figsize=(COL, COL * 0.92))
vmin, vmax = np.nanpercentile(DD, 2), np.nanpercentile(DD, 98)
im = ax.pcolormesh(np.rad2deg(TT), RR, np.log10(DD),
                   cmap='RdYlGn_r',
                   vmin=np.log10(vmin), vmax=np.log10(vmax),
                   shading='auto')
cb = plt.colorbar(im, ax=ax, label='$\\log_{10}$ relative drift '
                  '$|dK/d\\eta|/|K|$')
th_line = np.linspace(np.deg2rad(4), np.deg2rad(176), 400)
r_ergo = M_n + np.sqrt(M_n**2 - a_v**2 * np.cos(th_line)**2)
ax.plot(np.rad2deg(th_line), r_ergo, 'k--', lw=1.6, label='ergosphere')
ax.axhline(r_plus, color='k', lw=1.8, label='horizon $r_+$')
ax.set_yscale('log')
ax.set_ylim(rg[0], rg[-1])
ax.set_xlim(4, 176)
ax.set_xlabel('$\\theta$ [deg]')
ax.set_ylabel('$r$  (log)')
ax.set_title('drift of the quasi-constant $K_{NLO}(E_{eff})$\n'
             f'$s={a_v}$, $\\hat E={Eh_n}$, $A={A_v}$ '
             f'($E_{{eff}}={E_eff:.3f}$), $\\varepsilon={eps_v}$\n'
             'grey = beyond the freezing wall')
ax.set_facecolor('#cccccc')
ax.legend(loc='upper right')
savefig(fig, OUT, 'fig_thakurta_kerr_drift_map')
