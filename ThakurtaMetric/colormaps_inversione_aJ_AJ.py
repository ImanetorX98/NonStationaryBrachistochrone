# -*- coding: utf-8 -*-
"""
Due colormap di Delta_r = r_min^t - r_min^tau (stesso-lancio), colore =
Delta_r della cella:
  (a) piano (J, a): fattore conforme A=1 (Kerr), inversione ROTAZIONALE;
  (b) piano (J, A): spin a fisso, inversione CONFORME.
Contorno nero Delta_r=0 = luogo di inversione. Blu: t piu' fondo;
rosso: tau piu' fondo.

r_min di ciascun ramo = radice esterna della condizione di svolta
(p_r=0, H=0), equatoriale Thakurta-Kerr.
"""

import os
import sys
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x, **k):
        return x

set_style()
OUT = os.path.join(HERE, 'Thakurtafigures')
M, Eh = 1.0, 1.2
r_e = 2 * M

def geom(r, a, A):
    f = 1 - 2 * M / r
    Dl = r**2 - 2 * M * r + a**2
    P = r**2 + a**2 + 2 * M * a**2 / r
    b = 2 * M * a / r
    vb2 = 1 - A**2 * f / Eh**2
    Pb = P + A**2 * b**2 / Eh**2
    return f, Dl, P, b, vb2, Pb

def Htau0(r, a, A, J):
    f, Dl, P, b, vb2, Pb = geom(r, a, A)
    if vb2 <= 0:
        return np.nan
    ptp = J - A**2 * b / Eh
    return ptp * b * vb2 / Pb + np.sqrt(vb2 * Dl) / Pb * abs(ptp) - A**2 * f / Eh

def Ht0(r, a, A, J):
    f, Dl, P, b, vb2, Pb = geom(r, a, A)
    if vb2 <= 0:
        return np.nan
    return J * b * vb2 / Pb + np.sqrt(vb2 * Dl) / Pb * abs(J) - 1.0

def r_outer(func, a, A, J):
    r_plus = M + np.sqrt(max(M**2 - a**2, 0))   # scan fin dentro l'ergosfera
    rs = np.linspace(30.0, r_plus * 1.002, 3000)
    v = np.array([func(rr, a, A, J) for rr in rs])
    for i in range(len(rs) - 1):
        if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and v[i] * v[i + 1] < 0:
            return brentq(lambda rr: func(rr, a, A, J), rs[i + 1], rs[i])
    return np.nan

def has_orbit(func, a, A, J):
    """orbita esiste a r0=12 se min_pr H = H(r0,0,J) < 0 (regione ammessa)."""
    with np.errstate(invalid='ignore'):
        return func(12.0, a, A, J) < 0

def r_min_branch(func, a, A, J):
    """r_min: punto di svolta se scatter, r_e se catturata."""
    rt = r_outer(func, a, A, J)
    if np.isfinite(rt):
        return rt
    # nessuna svolta sopra r_e: catturata (grazing a r_e) se l'orbita esiste
    if has_orbit(func, a, A, J):
        return r_e
    return np.nan

def Dr(a, A, J):
    rt = r_min_branch(Htau0, a, A, J)
    rr = r_min_branch(Ht0, a, A, J)
    if np.isfinite(rt) and np.isfinite(rr):
        return rr - rt
    return np.nan

def mappa(xs, ys, get):
    Z = np.full((len(ys), len(xs)), np.nan)
    for i, yv in enumerate(tqdm(ys, desc='  righe', ncols=70, leave=False)):
        for j, xv in enumerate(xs):
            Z[i, j] = get(xv, yv)
    return Z

def plot(ax, X, Y, Z, xlabel, ylabel, title):
    vmax = min(1.5, np.nanmax(np.abs(Z)))
    im = ax.pcolormesh(X, Y, Z, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                       shading='auto')
    cb = plt.colorbar(im, ax=ax, extend='both',
                      label=r'$\Delta r = r_{\min}^{t}-r_{\min}^{\tau}$')
    try:
        ax.contour(X, Y, Z, levels=[0.0], colors='k', linewidths=1.4)
    except Exception:
        pass
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

# ---- (a) piano (J, a): A=1, inversione rotazionale ------------------
print("colormap (J, a), A=1 ...", flush=True)
J_a = np.linspace(0.3, 5.0, 60)
a_a = np.linspace(0.05, 0.95, 55)
Za = mappa(J_a, a_a, lambda Jv, av: Dr(av, 1.0, Jv))

# ---- (b) piano (J, A): a=0.9, inversione conforme -------------------
print("colormap (J, A), a=0.9 ...", flush=True)
a_fix = 0.9
J_A = np.linspace(5.0, 30.0, 60)
A_A = np.linspace(1.0, 1.5, 55)
Zb = mappa(J_A, A_A, lambda Jv, Av: Dr(a_fix, Av, Jv))

print(f"(J,a): Dr in [{np.nanmin(Za):+.2f},{np.nanmax(Za):+.2f}], "
      f"inversione: {'SI' if np.nanmin(Za)*np.nanmax(Za)<0 else 'no'}",
      flush=True)
print(f"(J,A): Dr in [{np.nanmin(Zb):+.2f},{np.nanmax(Zb):+.2f}], "
      f"inversione: {'SI' if np.nanmin(Zb)*np.nanmax(Zb)<0 else 'no'}",
      flush=True)

# --------------------------------------------------------------- figura
fig, (a1, a2) = plt.subplots(2, 1, figsize=(COL, 6.4))
plot(a1, J_a, a_a, Za, '$J$ (particle ang. mom.)', '$a$ (spin)',
     'Rotational: $\\Delta r$ in the $(J,a)$ plane ($A=1$)')
plot(a2, J_A, A_A, Zb, '$J$ (particle ang. mom.)', '$A$ (conformal factor)',
     f'Conformal: $\\Delta r$ in the $(J,A)$ plane ($a={a_fix}$)')
savefig(fig, OUT, 'fig_colormaps_inversione')
print("FATTO.", flush=True)
