# -*- coding: utf-8 -*-
"""
Drift map nel piano (A, r): il fattore conforme come asse.

Drift puntuale totale D(r; A) = |{K_NLO(E_eff), H_3D(A)} + eps a^2 S1|/|K|
a theta e momenti fissati, con eps = Hc*A. r LOGARITMICO in ordinata,
A in ascissa. Sovrapposte: orizzonte (linea), polo r = 2M (punteggiata),
muro di congelamento r_w(A) = 2M/(1 - Ehat^2/A^2) per A > Ehat
(tratteggiata; = polo di DE0). Verde poco drift -> rosso molto.
Output: Thakurtafigures/fig_thakurta_kerr_drift_map_A
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
a_v, Hc_v = 0.9, 0.02
J_v, pr_v, pth_v = 1.2, -0.3, 1.0
th_v = np.deg2rad(60.0)
r_plus = M_n + np.sqrt(M_n**2 - a_v**2)

r, th = sp.symbols('r theta', positive=True)
pr, pth = sp.symbols('p_r p_theta', real=True)
A = sp.Symbol('A', positive=True)

Sig = r**2 + a_v**2 * sp.cos(th)**2
Dl = r**2 - 2 * M_n * r + a_v**2
f_S = 1 - 2 * M_n * r / Sig
b = 2 * M_n * a_v * r * sp.sin(th)**2 / Sig
G = (r**2 + a_v**2 + 2 * M_n * a_v**2 * r * sp.sin(th)**2 / Sig) \
    * sp.sin(th)**2
vbS2 = 1 - A**2 * f_S / Eh_n**2
Gb = G + A**2 * b**2 / Eh_n**2
php0 = b * vbS2 / Gb
R2 = vbS2 * Dl * sp.sin(th)**2 / Gb
ptphi = J_v - A**2 * b / Eh_n
H3 = ptphi * php0 + sp.sqrt(R2) * sp.sqrt((Dl / Sig) * pr**2
                                          + pth**2 / Sig + ptphi**2 / Gb) \
    - A**2 * f_S / Eh_n

E_eff = Eh_n / A
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
eps = Hc_v * A
drift = sp.Abs(PB + eps * a_v**2 * S1) / sp.Abs(K)

print("lambdify (cse)...")
Dr_f = sp.lambdify((r, th, pr, pth, A), drift, 'numpy', cse=True)
vb_f = sp.lambdify((r, th, A), vbS2, 'numpy')

Ag = np.linspace(0.7, 1.6, 340)
rg = np.logspace(np.log10(r_plus * 1.02), np.log10(40.0), 340)
AA, RR = np.meshgrid(Ag, rg)
with np.errstate(invalid='ignore', divide='ignore'):
    DD = Dr_f(RR, th_v, pr_v, pth_v, AA)
    VV = vb_f(RR, th_v, AA)
DD = np.where(VV > 0, DD, np.nan)
DD = np.where(np.isfinite(DD), DD, np.nan)
print(f"drift: min {np.nanmin(DD):.2e}  mediana {np.nanmedian(DD):.2e}")

fig, ax = plt.subplots(figsize=(COL, COL * 0.92))
vmin, vmax = np.nanpercentile(DD, 2), np.nanpercentile(DD, 98)
im = ax.pcolormesh(AA, RR, np.log10(DD), cmap='RdYlGn_r',
                   vmin=np.log10(vmin), vmax=np.log10(vmax),
                   shading='auto')
plt.colorbar(im, ax=ax,
             label='$\\log_{10}$ relative drift $|dK/d\\eta|/|K|$')
ax.axhline(r_plus, color='k', lw=1.8, label='horizon $r_+$')
ax.axhline(2 * M_n, color='k', lw=1.1, ls=':', label='pole $r=2M$ ($D_0$)')
A_w = np.linspace(Eh_n * 1.001, Ag[-1], 200)
r_w = 2 * M_n / (1 - Eh_n**2 / A_w**2)
ax.plot(A_w, r_w, 'k--', lw=1.6,
        label='wall $r_w(A)$ (pole of $DE_0$)')
ax.axvline(Eh_n, color='b', lw=0.8, ls=':')
ax.text(Eh_n + 0.01, 25, '$A=\\hat E$', fontsize=7, color='b')
ax.set_yscale('log')
ax.set_xlim(Ag[0], Ag[-1])
ax.set_ylim(rg[0], rg[-1])
ax.set_xlabel('$A$ (conformal factor)')
ax.set_ylabel('$r$  (log)')
ax.set_title('drift of $K_{NLO}(E_{eff})$ in the $(A, r)$ plane, '
             f'$\\theta={np.rad2deg(th_v):.0f}^\\circ$\n'
             f'$s={a_v}$, $\\hat E={Eh_n}$, $\\varepsilon=H_cA$ '
             f'($H_c={Hc_v}$); grey = beyond the wall')
ax.set_facecolor('#cccccc')
ax.legend(loc='upper left')
savefig(fig, OUT, 'fig_thakurta_kerr_drift_map_A')
