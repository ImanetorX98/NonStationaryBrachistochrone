# -*- coding: utf-8 -*-
"""
Figure di validazione — Thakurta (Schwarzschild e Kerr conforme).

Output in Thakurtafigures/:
  fig_thakurta_cattura        cattura da espansione: orbita vs superficie
                              di congelamento + regione accessibile (eta,r)
  fig_thakurta_kerr_superfici (A, r): orizzonte rigido vs congelamento
                              che respira; shell accessibile
  fig_thakurta_kerr_residui   flusso H vs forme chiuse Kerr (A=1),
                              residui entrambi i rami
e in ../PaperFigures/:
  fig_indicatrici             figura concettuale: le indicatrici ellittiche
                              nei tre casi (FLRW, Vaidya, Thakurta-Kerr)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from paper_style import COL, DCOL, set_style, savefig

set_style()
OUT = os.path.join(HERE, 'Thakurtafigures')
OUTP = os.path.join(os.path.dirname(HERE), 'PaperFigures')

def salva(fig, nome, cart=OUT):
    savefig(fig, cart, nome)

M_n, E_n, J_n = 1.0, 1.2, 1.3

# ------------------------------------------------------- fig 1: cattura
print('fig 1: cattura da espansione (Thakurta-Schwarzschild)...')
H_c, r_l = 0.02, 4.0
eta0 = -1.0 / H_c

ee, rr, pr = sp.symbols('ee rr p_r', real=True)
a_s = -1 / (H_c * ee)
f_s = 1 - 2 * M_n / rr
v_s = sp.sqrt(1 - a_s**2 * f_s / E_n**2)
H_tau = v_s * sp.sqrt(f_s) * sp.sqrt(f_s * pr**2 + J_n**2 / rr**2) \
    - a_s**2 * f_s / E_n
dHdp = sp.lambdify((ee, rr, pr), sp.diff(H_tau, pr), 'numpy')
dHdr = sp.lambdify((ee, rr, pr), sp.diff(H_tau, rr), 'numpy')
Hfun = sp.lambdify((ee, rr, pr), H_tau, 'numpy')

pg = np.linspace(-200, 200, 400001)
Hg = Hfun(eta0, r_l, pg)
roots = [brentq(lambda p: Hfun(eta0, r_l, p), pg[i], pg[i + 1])
         for i in range(len(pg) - 1)
         if np.isfinite(Hg[i]) and np.isfinite(Hg[i + 1])
         and Hg[i] * Hg[i + 1] <= 0]
p0 = [p_ for p_ in roots if dHdp(eta0, r_l, p_) > 0][-1]
s = solve_ivp(lambda e_, y: [dHdp(e_, y[0], y[1]), -dHdr(e_, y[0], y[1])],
              [eta0, -1e-3], [r_l, p0], rtol=1e-11, atol=1e-13,
              method='DOP853', dense_output=True)
eta_c = s.t[-1]

def a_n(e_):
    return -1.0 / (H_c * e_)

def r_freeze(e_):
    t_ = E_n**2 / a_n(e_)**2
    return 2 * M_n / (1 - t_) if t_ < 1 else np.inf

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(COL, 5.8))
eg = np.linspace(eta0, eta_c, 400)
ax1.plot(eg, s.sol(eg)[0], 'C0-', lw=1.5, label='orbit (launched outward)')
eg2 = np.linspace(eta0, -0.5, 2000)
rf = np.array([r_freeze(e_) for e_ in eg2])
ax1.plot(eg2, rf, 'C3--', lw=1.2,
         label='freezing $a^2f=\\hat E^2$')
eg3 = np.linspace(eta_c, -0.5, 300)
ax1.plot(eg3, [r_freeze(e_) for e_ in eg3], 'C0:', lw=2.0,
         label='after contact: rides the surface')
ax1.axhline(2 * M_n, color='k', lw=0.8)
ax1.text(-48, 2.1, 'horizon $2M$', fontsize=6.5)
ax1.plot(eta_c, s.y[0, -1], 'ko', ms=5)
ax1.set_ylim(1.5, 16)
ax1.set_xlabel('$\\eta$')
ax1.set_ylabel('$r$')
ax1.set_title('capture by expansion (Thakurta-Schw., $a=-1/H\\eta$,\n'
              f'$\\hat E={E_n}$, $J={J_n}$): contact at '
              f'$\\eta={eta_c:.1f}$, $r={s.y[0,-1]:.2f}$')
ax1.legend()

rgrid = np.linspace(1.6, 16, 400)
egrid = np.linspace(eta0, -0.5, 400)
EE, RR = np.meshgrid(egrid, rgrid)
acc = (1 - 2 * M_n / RR) * a_n(EE)**2 < E_n**2
acc &= RR > 2 * M_n
ax2.contourf(EE, RR, acc, levels=[0.5, 1.5], colors=['#cce5ff'])
ax2.plot(eg, s.sol(eg)[0], 'C0-', lw=1.5)
ax2.plot(eg2, rf, 'C3--', lw=1.2)
ax2.axhline(2 * M_n, color='k', lw=0.8)
ax2.set_ylim(1.5, 16)
ax2.set_xlabel('$\\eta$')
ax2.set_ylabel('$r$')
ax2.set_title('rail-accessible region (blue):\n'
              'the horizon-freezing shell closes')
salva(fig, 'fig_thakurta_cattura')

# ------------------------------------- fig 2: superfici Thakurta-Kerr
print('fig 2: superfici critiche Thakurta-Kerr...')
s_n = 0.9
r_p = M_n + np.sqrt(M_n**2 - s_n**2)
fig, ax = plt.subplots(figsize=(COL, COL * 0.85))
Ag = np.linspace(1.0, 5.0, 400)
rf_k = np.array([2 * M_n / (1 - E_n**2 / Av**2) if Av > E_n else np.inf
                 for Av in Ag])
ax.fill_between(Ag, r_p, np.minimum(rf_k, 40), color='#cce5ff',
                label='accessible shell')
ax.axhline(r_p, color='k', lw=1.6,
           label='horizon $r_+$ (RIGID: conformally invariant)')
ax.axhline(2 * M_n, color='C2', lw=1.0, ls=':',
           label='ergosphere $r_e=2M$ (crossed regularly)')
ax.plot(Ag, rf_k, 'C3--', lw=1.4,
        label='freezing $A^2f=\\hat E^2$ (BREATHES)')
ax.axvline(E_n, color='C3', lw=0.6, ls=':')
ax.text(E_n + 0.04, 20, '$A=\\hat E$', fontsize=7, color='C3')
ax.set_yscale('log')
ax.set_ylim(1.2, 40)
ax.set_xlabel('scale factor $A(\\eta)$')
ax.set_ylabel('$r$')
ax.set_title(f'Thakurta-Kerr ($s={s_n}$, $\\hat E={E_n}$): '
             '$R^2=\\Delta\\bar v^2/\\bar P$\n'
             'two critical surfaces: one rigid, one breathing')
ax.legend(loc='upper right')
salva(fig, 'fig_thakurta_kerr_superfici')

# ------------------------------------- fig 3: residui Kerr statico
print('fig 3: residui forme chiuse Kerr (A=1)...')
r_, pr_ = sp.symbols('r_ p_r_', positive=True), sp.Symbol('p_r_', real=True)
r_ = sp.Symbol('r_', positive=True)
f_k = 1 - 2 * M_n / r_
Dl = r_**2 - 2 * M_n * r_ + s_n**2
P_k = r_**2 + s_n**2 + 2 * M_n * s_n**2 / r_
vb2 = 1 - f_k / E_n**2                    # A=1
Pb = P_k + (2 * M_n * s_n / r_)**2 / E_n**2
php0 = (2 * M_n * s_n / r_) * vb2 / Pb
R2 = Dl * vb2 / Pb
pphi = sp.Symbol('p_phi', real=True)
R_ = sp.sqrt(R2)
H_eta = pphi * php0 + R_ * sp.sqrt((Dl / r_**2) * pr_**2 + pphi**2 / Pb) - 1
ptphi = pphi - 2 * M_n * s_n / (r_ * E_n)
H_tauK = ptphi * php0 + R_ * sp.sqrt((Dl / r_**2) * pr_**2
                                     + ptphi**2 / Pb) - f_k / E_n

fig, ax = plt.subplots(figsize=(COL, COL * 0.82))
rgrid = np.linspace(3.0, 14.0, 34)
for Hx, Kfun, st, lab in (
        (H_tauK, lambda rv, fv: J_n, 'o-', '$\\tau$ branch'),
        (H_eta, lambda rv, fv: (fv * J_n + 2 * M_n * s_n / rv) / E_n,
         's--', '$t$ branch')):
    Hn = Hx.subs(pphi, J_n)
    dHdp_ = sp.lambdify((r_, pr_), sp.diff(Hn, pr_), 'numpy')
    dHdJ_ = sp.lambdify((r_, pr_), sp.diff(Hx, pphi).subs(pphi, J_n),
                        'numpy')
    Hf = sp.lambdify((r_, pr_), Hn, 'numpy')
    res = []
    for rv in rgrid:
        pg = np.linspace(-80, 80, 320001)
        Hg = Hf(rv, pg)
        roots = [brentq(lambda p: Hf(rv, p), pg[i], pg[i + 1])
                 for i in range(len(pg) - 1)
                 if np.isfinite(Hg[i]) and np.isfinite(Hg[i + 1])
                 and Hg[i] * Hg[i + 1] <= 0]
        p_in = [p_ for p_ in roots if dHdp_(rv, p_) < 0][0]
        num = abs(dHdJ_(rv, p_in) / dHdp_(rv, p_in))
        fv = 1 - 2 * M_n / rv
        wv = E_n**2 - fv
        Dv = rv**2 - 2 * M_n * rv + s_n**2
        Kv = Kfun(rv, fv)
        cf = Kv * rv * np.sqrt(wv * fv) / (Dv * np.sqrt(Dv - Kv**2 * wv))
        res.append(abs(num - cf) / cf)
    ax.semilogy(rgrid, np.maximum(res, 1e-17), st, ms=3, lw=0.9, label=lab)
ax.set_xlabel('$r/M$')
ax.set_ylabel('relative residual $|\\Delta(d\\phi/dr)|$')
ax.set_title('Thakurta-Kerr, static limit $A=1$: Hamilton flow\n'
             f'vs closed-form equatorial Kerr ($s={s_n}$, $E={E_n}$, '
             f'$J={J_n}$)')
ax.legend()
salva(fig, 'fig_thakurta_kerr_residui')

# ------------------------------------- fig 4: indicatrici (concettuale)
print('fig 4: indicatrici (figura concettuale, PaperFigures)...')
fig, axs = plt.subplots(3, 1, figsize=(COL, 8.4))

# (a) FLRW: shrinking circles v(eta)
ax = axs[0]
for aa, c in zip((1.0, 1.5, 2.0, 2.75), plt.cm.viridis([0.2, 0.45, 0.7, 0.9])):
    v = np.sqrt(max(1 - aa**2 / 3.0**2, 0))
    ax.add_patch(plt.Circle((0, 0), v, fill=False, color=c, lw=1.6,
                            label=f'$a={aa}$'))
ax.plot(0, 0, 'k+', ms=8)
ax.set_title('FLRW ($\\hat E=3$): circles,\nno wind; freezing $a\\to\\hat E$')
ax.legend(fontsize=6, loc='upper right')
ax.set_xlim(-1.15, 1.15)
ax.set_ylim(-1.15, 1.15)

# (b) Vaidya EF: ingoing radial wind f - E^2
ax = axs[1]
for rv, c in zip((10.0, 4.0, 2.5), plt.cm.plasma([0.15, 0.5, 0.8])):
    fv = 1 - 2.0 / rv
    wv = E_n**2 - fv
    cx = fv - E_n**2
    ax.add_patch(Ellipse((cx, 0), 2 * E_n * np.sqrt(wv), 2 * np.sqrt(wv),
                         fill=False, color=c, lw=1.6, label=f'$r={rv}$'))
    ax.plot(cx, 0, '.', color=c, ms=5)
ax.plot(0, 0, 'k+', ms=8)
ax.axvline(0, color='k', lw=0.5)
ax.set_title('Vaidya (EF, $E=1.2$): ingoing radial\nwind $f-E^2<0$; '
             '$m(v)$ moves everything')
ax.legend(fontsize=6, loc='upper right')
ax.set_xlim(-3.6, 1.2)
ax.set_ylim(-2.0, 2.0)

# (c) Thakurta-Kerr: angular wind (dragging), ellipse shrinks with A
ax = axs[2]
rv = 2.2   # dentro quasi-ergosfera per vento visibile
fv = 1 - 2.0 / rv
Dv = rv**2 - 2 * rv + s_n**2
Pv = rv**2 + s_n**2 + 2 * s_n**2 / rv
for Av, c in zip((1.0, 2.0, 3.0), plt.cm.viridis([0.25, 0.55, 0.85])):
    vb = 1 - Av**2 * fv / E_n**2
    Pbv = Pv + Av**2 * (2 * s_n / rv)**2 / E_n**2
    ph0 = (2 * s_n / rv) * vb / Pbv
    R2v = Dv * vb / Pbv
    ax.add_patch(Ellipse((0, ph0 * rv), 2 * np.sqrt(R2v * Dv) / rv,
                         2 * np.sqrt(R2v / Pbv) * rv, fill=False, color=c,
                         lw=1.6, label=f'$A={Av}$'))
    ax.plot(0, ph0 * rv, '.', color=c, ms=5)
ax.plot(0, 0, 'k+', ms=8)
ax.axhline(0, color='k', lw=0.5)
ax.set_title(f'Thakurta-Kerr ($r={rv}$, $s={s_n}$): angular\nwind '
             '(dragging); $A(\\eta)$ shrinks it')
ax.legend(fontsize=6, loc='upper right')
ax.set_xlim(-1.3, 1.3)
ax.set_ylim(-0.6, 1.5)
for ax in axs:
    ax.set_xlabel("$r'$")
    ax.set_ylabel("$r\\,\\phi'$")
    ax.set_aspect('equal')
fig.suptitle('Rail indicatrices (linear constraint $\\cap$ normalization'
             '\n= ellipse): the unifying formalism')
salva(fig, 'fig_indicatrici', OUTP)
print('\nFATTO.')
