# -*- coding: utf-8 -*-
"""
Same-launch depth comparison in equatorial Thakurta-Kerr, (J, A) plane, in TWO normalisations
of the costate: (1) the protocol declared in the paper, t vs tau at J_t = J_tau (H_t = h - A);
(2) eta vs tau at J_eta = J_tau, that is J_t = A J_tau (H_eta = h - 1).  Only (2) inverts the ordering:
the "conformal" inversion is an effect of the costate normalisation, not of the rail.
(An earlier version used H_eta and called it H_t.)
For (1) there is no inversion in the prograde sector: Paper II, prop:same-launch.

EARLIER DESCRIPTION (OBSOLETE, kept for the record; it applies to normalisation (2), not to H_t):
This time with the RIGHT RANGE OF J: the inversion lives at J ~ 16-24, outside the
old scan [3.9, 5.5]. It shows the red region (tau deeper).

r_min of each branch = outer root of the turning condition
(p_r = 0, H = 0), computed directly (no integration):
  tau:  H_tau(r,0,J) = ptp*php0 + sqrt(R2)|ptp|/sqrt(Pb) - A^2 f/Ehat = 0
  t  :  H_t  (r,0,J) = J*php0   + sqrt(R2)|J|/sqrt(Pb)   - 1        = 0
Analytic inversion curve (r_min^tau = r_min^t):
  cond(r,A) = A^2 b Q - Pb (Ehat - A^2 f) = 0 ,  Q = b vb2 + sqrt(vb2 Dl)
with  cond(A=1) < 0,  cond(A_freeze) = Ehat r Dl (Ehat-1)/(r-2M) > 0
=> A_inv in (1, A_freeze): the inversion always exists and is reachable.
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

# ---- turning conditions (numerical) --------------------------------
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

def Heta0(r, A, J):
    """conformal-time branch H_eta = h - 1 at p_r = 0."""
    f, Dl, P, b, vb2, Pb = geom(r, A)
    if vb2 <= 0:
        return np.nan
    php0 = b * vb2 / Pb
    S = np.sqrt(vb2 * Dl) / Pb
    return J * php0 + S * abs(J) - 1.0

def Ht0(r, A, J):
    """coordinate-time branch H_t = h - A at p_r = 0 (paper eq:Ht-tk); = H_eta at A = 1."""
    return Heta0(r, A, J) + 1.0 - A

def r_outer(func, A, J):
    rs = np.linspace(20.0, r_e + 1e-4, 3000)
    vals = np.array([func(rr, A, J) for rr in rs])
    for i in range(len(rs) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) \
                and vals[i] * vals[i + 1] < 0:
            return brentq(lambda rr: func(rr, A, J), rs[i + 1], rs[i])
    return np.nan

# ---- griglia (J, A): declared protocol t/tau at J_t = J_tau, and eta/tau at J_eta = J_tau ----
J_g = np.linspace(5.0, 30.0, 60)
A_g = np.linspace(1.0, 1.5, 55)
Dr = {'t': np.full((len(A_g), len(J_g)), np.nan), 'eta': np.full((len(A_g), len(J_g)), np.nan)}
for i, Av in enumerate(A_g):
    for j, Jv in enumerate(J_g):
        rtau = r_outer(Htau0, Av, Jv)
        for k, fn in (('t', Ht0), ('eta', Heta0)):
            rx = r_outer(fn, Av, Jv)
            if np.isfinite(rtau) and np.isfinite(rx):
                Dr[k][i, j] = rx - rtau
for k in ('t', 'eta'):
    print(f"{k}/tau: cells with tau deeper (Dr>0): {int(np.nansum(Dr[k] > 0))}/"
          f"{np.sum(np.isfinite(Dr[k]))};  Dr range [{np.nanmin(Dr[k]):+.3f}, {np.nanmax(Dr[k]):+.3f}]")

# ---- algebraic coincidence conditions of the two outer turning points (prograde sector) ----
#   eta/tau:  A^2 b Q - Pb (Ehat - A^2 f) = 0      t/tau:  A b Q - Pb (Ehat - A f) = 0
r_s, A_s = sp.symbols('r A', positive=True)
f_s = 1 - 2 * M / r_s
Dl_s = r_s**2 - 2 * M * r_s + a**2
P_s = r_s**2 + a**2 + 2 * M * a**2 / r_s
b_s = 2 * M * a / r_s
vb2_s = 1 - A_s**2 * f_s / Eh**2
Pb_s = P_s + A_s**2 * b_s**2 / Eh**2
Q_s = b_s * vb2_s + sp.sqrt(vb2_s * Dl_s)
cond = A_s**2 * b_s * Q_s - Pb_s * (Eh - A_s**2 * f_s)
cond_t = A_s * b_s * Q_s - Pb_s * (Eh - A_s * f_s)
condf = sp.lambdify((r_s, A_s), cond, 'numpy')
condtf = sp.lambdify((r_s, A_s), cond_t, 'numpy')
Jinvf = sp.lambdify((r_s, A_s), Pb_s / Q_s, 'numpy')
A_curve, J_curve = [], []
for Av in np.linspace(A_g[0] + 0.02, A_g[-1], 120):
    rs = np.linspace(2.05, 40.0, 8000)
    with np.errstate(invalid='ignore'):
        cv = np.array([condf(rr, Av) for rr in rs])
    for k in range(len(rs) - 1):
        if np.isfinite(cv[k]) and np.isfinite(cv[k + 1]) and cv[k] * cv[k + 1] < 0:
            rstar = brentq(lambda rr: condf(rr, Av), rs[k], rs[k + 1])
            A_curve.append(Av); J_curve.append(float(Jinvf(rstar, Av)))
            break
A_curve, J_curve = np.array(A_curve), np.array(J_curve)
with np.errstate(invalid='ignore'):
    ct_max = np.nanmax([condtf(rr, Av) for Av in np.linspace(1.0, 1.5, 120)
                        for rr in np.linspace(2.05, 40.0, 800)])
print(f"max of the t/tau coincidence condition on the grid: {ct_max:+.3f} (negative: no coincidence)")

# --------------------------------------------------------------- figure
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COL, 2 * COL * 0.85))
vmax = 1.5
for ax, k, ttl in ((ax1, 't', 'declared protocol: $t$ vs $\\tau$ at $J_t=J_\\tau$'),
                   (ax2, 'eta', 'conformal normalisation: $\\eta$ vs $\\tau$ at $J_\\eta=J_\\tau$ ($J_t=AJ_\\tau$)')):
    im = ax.pcolormesh(J_g, A_g, Dr[k], cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
    plt.colorbar(im, ax=ax, extend='min', label='$\\Delta r$ (branch $-$ $\\tau$)')
    if np.nanmax(Dr[k]) > 0:
        ax.contour(J_g, A_g, Dr[k], levels=[0.0], colors='k', linewidths=1.4)
    ax.set_xlabel('$J$'); ax.set_ylabel('$A$')
    ax.set_title(ttl + f"\n$\\tau$ deeper in {int(np.nansum(Dr[k] > 0))} of {int(np.sum(np.isfinite(Dr[k])))} cells",
                 fontsize=6.6)
msk = J_curve <= J_g[-1]
ax2.plot(J_curve[msk], A_curve[msk], 'w--', lw=1.4, label='analytic coincidence condition')
ax2.legend(loc='lower right', fontsize=6)
fig.suptitle(f'Thakurta--Kerr, $a={a}$, $\\hat E={Eh}$: the same-launch ordering depends on the costate normalisation',
             fontsize=6.6)
savefig(fig, OUT, 'fig_thakurta_kerr_inversion_AJ')
print('DONE.')
