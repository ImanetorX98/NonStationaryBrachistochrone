# -*- coding: utf-8 -*-
"""
Colormap Delta_r(a, J) per brachistocrone t/tau a ESTREMI FISSI (Kerr
equatoriale, flusso di Hamilton -- pesante). Asse J = momento del ramo
tau, che fissa l'apertura angolare Delta_phi degli estremi (entrambi a
r0). Per ogni (a, J_tau):
  - ramo tau con J_tau: (Delta_phi, r_min^tau) dall'integrazione;
  - ramo t agli STESSI estremi: spara J_t -> Delta_phi_t = Delta_phi;
    prende la soluzione SHALLOW (r_min max = brachistocrona minima);
  - Delta_r = r_min^t - r_min^tau.
Atteso (dal caso conforme e dallo studio a Phi fisso): Delta_r > 0
ovunque -> nessuna inversione a estremi fissi anche con spin.
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
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
M, E, r0 = 1.0, 1.2, 6.0

r, pr, J, a = sp.symbols('r p_r J a', real=True)
f = 1 - 2 * M / r
Dl = r**2 - 2 * M * r + a**2
P = r**2 + a**2 + 2 * M * a**2 / r
b = 2 * M * a / r
vb2 = 1 - f / E**2
Pb = P + b**2 / E**2
php0 = b * vb2 / Pb
R2 = vb2 * Dl / Pb

def build(branch):
    if branch == 'tau':
        ptp = J - b / E
        H = ptp * php0 + sp.sqrt(R2) * sp.sqrt((Dl / r**2) * pr**2
                                               + ptp**2 / Pb) - f / E
    else:
        H = J * php0 + sp.sqrt(R2) * sp.sqrt((Dl / r**2) * pr**2
                                             + J**2 / Pb) - 1
    return (sp.lambdify((a, r, pr, J), H, 'numpy'),
            sp.lambdify((a, r, pr, J), sp.diff(H, pr), 'numpy'),
            sp.lambdify((a, r, pr, J), sp.diff(H, r), 'numpy'),
            sp.lambdify((a, r, pr, J), sp.diff(H, J), 'numpy'))

FN = {'tau': build('tau'), 't': build('t')}

def orbita(branch, av, Jv):
    Hf, dHp, dHr, dHJ = FN[branch]
    pg = np.linspace(-50, 50, 2001)
    with np.errstate(invalid='ignore'):
        Hv = Hf(av, r0, pg, Jv)
    roots = [brentq(lambda p: Hf(av, r0, p, Jv), pg[i], pg[i + 1])
             for i in range(len(pg) - 1)
             if np.isfinite(Hv[i]) and np.isfinite(Hv[i + 1])
             and Hv[i] * Hv[i + 1] < 0]
    ing = [p for p in roots if dHp(av, r0, p, Jv) < 0]
    if not ing:
        return None
    p0 = min(ing)
    ev = lambda t_, y: y[0] - r0
    ev.terminal, ev.direction = True, 1
    s = solve_ivp(lambda t_, y: [dHp(av, y[0], y[1], Jv),
                                 -dHr(av, y[0], y[1], Jv),
                                 dHJ(av, y[0], y[1], Jv)],
                  [0, 5000], [r0, p0, 0.0], rtol=1e-9, atol=1e-11,
                  method='DOP853', events=[ev], max_step=2.0)
    if not len(s.t_events[0]):
        return None
    return s.y_events[0][0][2], s.y[0].min()      # (Dphi, r_min)

def t_rmin(av, Dphi_target):
    """spara J_t: Dphi_t = target; shallow (r_min max)."""
    Js = np.linspace(0.5, 16.0, 45)
    dat = []
    for Jv in Js:
        o = orbita('t', av, Jv)
        if o is not None:
            dat.append((Jv, o[0], o[1]))
    sols = []
    for i in range(len(dat) - 1):
        if (dat[i][1] - Dphi_target) * (dat[i + 1][1] - Dphi_target) < 0:
            try:
                Js2 = brentq(lambda x: orbita('t', av, x)[0] - Dphi_target,
                             dat[i][0], dat[i + 1][0])
                o = orbita('t', av, Js2)
                if o is not None:
                    sols.append(o[1])
            except (ValueError, TypeError):
                pass
    return max(sols) if sols else np.nan

def Dr(av, Jtau):
    o = orbita('tau', av, Jtau)
    if o is None:
        return np.nan
    Dphi, rm_tau = o
    rm_t = t_rmin(av, Dphi)
    return rm_t - rm_tau if np.isfinite(rm_t) else np.nan

a_g = np.linspace(0.05, 0.95, 20)
J_g = np.linspace(1.2, 4.0, 18)          # momento del ramo tau
Z = np.full((len(a_g), len(J_g)), np.nan)
for i, av in enumerate(tqdm(a_g, desc='spin a', ncols=70)):
    for j, Jt in enumerate(J_g):
        Z[i, j] = Dr(av, Jt)

print(f"Delta_r in [{np.nanmin(Z):+.3f}, {np.nanmax(Z):+.3f}]", flush=True)
print(">>> " + ("INVERSIONE (cambia segno)" if np.nanmin(Z)*np.nanmax(Z) < 0
                else "NESSUNA inversione: Delta_r > 0 ovunque"), flush=True)

# --------------------------------------------------------------- figura
fig, ax = plt.subplots(figsize=(COL, COL * 0.9))
im = ax.pcolormesh(J_g, a_g, Z, cmap='viridis', shading='auto')
plt.colorbar(im, ax=ax, label=r'$\Delta r = r_{\min}^{t}-r_{\min}^{\tau}$')
if np.nanmin(Z) * np.nanmax(Z) < 0:
    try:
        ax.contour(J_g, a_g, Z, levels=[0], colors='r', linewidths=1.4)
    except Exception:
        pass
ax.set_xlabel(r'$J$ ($\tau$-branch momentum $\to$ endpoint angle)')
ax.set_ylabel('$a$ (black-hole spin)')
ax.set_title('Fixed symmetric endpoints, Kerr: $\\Delta r$ in $(J,a)$\n'
             'no inversion ($\\Delta r>0$): spin does not invert either')
savefig(fig, HERE, 'fig_colormap_spin_estremi_fissi')
print("FATTO.", flush=True)
