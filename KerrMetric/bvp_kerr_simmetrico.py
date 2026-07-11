# -*- coding: utf-8 -*-
"""
Inversione ROTAZIONALE a ESTREMI FISSI (simmetrici) in Kerr equatoriale.

Estremi A=(r0,-Phi), B=(r0,+Phi). Per ciascun ramo (t/eta, tau) si spara
il momento J (conservato) finche' la brachistocrona da A a B (dip verso
il periasse e ritorno a r0) spazza Delta_phi = 2 Phi; r_min = periasse.
Delta_r = r_min^t - r_min^tau vs spin a: si cerca l'inversione (Delta_r=0).

Hamiltoniane equatoriali (A=1, rotaia E):
  f=1-2M/r, Dl=r^2-2Mr+a^2, P=r^2+a^2+2Ma^2/r, b=2Ma/r
  vb2=1-f/E^2, Pb=P+b^2/E^2, php0=b vb2/Pb, R2=vb2 Dl/Pb
  tau: ptp=J-b/E; H=ptp php0 + sqrt(R2) sqrt((Dl/r^2)pr^2+ptp^2/Pb) - f/E
  t  : H=J php0   + sqrt(R2) sqrt((Dl/r^2)pr^2+J^2/Pb)   - 1
"""

import os
import sys
import time
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
try:
    from tqdm import tqdm
except ImportError:                          # fallback: no-op
    def tqdm(x, **k):
        return x

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, E, r0, Phi = 1.0, 1.2, 6.0, 0.9
_t0 = time.time()

def log(msg):
    print(f"[{time.time()-_t0:6.1f}s] {msg}", flush=True)

r, pr, J, a = sp.symbols('r p_r J a', real=True)
f = 1 - 2 * M / r
Dl = r**2 - 2 * M * r + a**2
P = r**2 + a**2 + 2 * M * a**2 / r
b = 2 * M * a / r
vb2 = 1 - f / E**2
Pb = P + b**2 / E**2
php0 = b * vb2 / Pb
R2 = vb2 * Dl / Pb
rad = sp.sqrt((Dl / r**2) * pr**2)

def build(branch):
    if branch == 'tau':
        ptp = J - b / E
        H = ptp * php0 + sp.sqrt(R2) * sp.sqrt((Dl / r**2) * pr**2
                                               + ptp**2 / Pb) - f / E
    else:
        H = J * php0 + sp.sqrt(R2) * sp.sqrt((Dl / r**2) * pr**2
                                             + J**2 / Pb) - 1
    Hf = sp.lambdify((a, r, pr, J), H, 'numpy')
    dHp = sp.lambdify((a, r, pr, J), sp.diff(H, pr), 'numpy')
    dHr = sp.lambdify((a, r, pr, J), sp.diff(H, r), 'numpy')
    dHJ = sp.lambdify((a, r, pr, J), sp.diff(H, J), 'numpy')
    return Hf, dHp, dHr, dHJ

FN = {'tau': build('tau'), 't': build('t')}

def orbita(branch, av, Jv):
    """da r0, ramo ENTRANTE, fino al ritorno a r0. Ritorna (Dphi, r_min)."""
    Hf, dHp, dHr, dHJ = FN[branch]
    pg = np.linspace(-60, 60, 4001)        # root-finding leggero
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
    ev_back = lambda t_, y: y[0] - r0       # ritorno a r0 uscendo
    ev_back.terminal, ev_back.direction = True, 1
    s = solve_ivp(lambda t_, y: [dHp(av, y[0], y[1], Jv),
                                 -dHr(av, y[0], y[1], Jv),
                                 dHJ(av, y[0], y[1], Jv)],
                  [0, 6000], [r0, p0, 0.0], rtol=1e-10, atol=1e-12,
                  method='DOP853', events=[ev_back], max_step=1.0)
    if not len(s.t_events[0]):
        return None
    return s.y_events[0][0][2], s.y[0].min()   # (Dphi, r_min)

def rmin_bvp(branch, av):
    """spara J: Dphi(J)=2 Phi. Ritorna r_min (e J)."""
    # t deflette molto piu' di tau -> serve J grande (J_opt^t ~ 8 gia' a a=0)
    Js = np.linspace(0.5, 14.0, 60)        # scan grezzo, poi brentq
    dat = [(Jv, orbita(branch, av, Jv))
           for Jv in tqdm(Js, desc=f'  scan J ({branch}, a={av:+.2f})',
                          leave=False, ncols=70)]
    dat = [(Jv, o[0], o[1]) for Jv, o in dat if o is not None]
    # la BVP puo' avere piu' soluzioni (deep/shallow): raccogli tutte
    # e prendi la SUPERFICIALE (r_min piu' grande = brachistocrona minima)
    sols = []
    for i in range(len(dat) - 1):
        g0, g1 = dat[i][1] - 2 * Phi, dat[i + 1][1] - 2 * Phi
        if g0 * g1 < 0:
            Jstar = brentq(lambda Jv: orbita(branch, av, Jv)[0] - 2 * Phi,
                           dat[i][0], dat[i + 1][0])
            o = orbita(branch, av, Jstar)
            if o is not None:
                sols.append((o[1], Jstar))
    if not sols:
        return None, None
    return max(sols, key=lambda s: s[0])       # r_min massimo (shallow)

log(f"Kerr equatoriale, estremi simmetrici (r0={r0}, +/-{Phi}), E={E}")
log("Delta_r = r_min^t - r_min^tau  vs spin a")
print(f"{'a':>6} {'r_min^tau':>10} {'r_min^t':>9} {'Delta_r':>9}", flush=True)
a_grid = np.linspace(-0.9, 0.9, 19)      # prograde E retrogrado (|a|<1)
res = []
bar = tqdm(a_grid, desc='spin a', ncols=70)
for av in bar:
    rt, _ = rmin_bvp('tau', av)
    rr, _ = rmin_bvp('t', av)
    if rt is None or rr is None:
        bar.write(f"  a={av:+.2f}: nessuna soluzione BVP")
        continue
    res.append((av, rt, rr, rr - rt))
    bar.set_postfix(Dr=f'{rr-rt:+.3f}', rtau=f'{rt:.2f}', rt=f'{rr:.2f}')
    print(f"{av:6.2f} {rt:10.4f} {rr:9.4f} {rr-rt:+9.4f}", flush=True)

res = np.array(res)
c = None
if len(res) > 2:
    a_g, dr = res[:, 0], res[:, 3]
    c = np.polyfit(a_g, dr, 2)
    print(f"\nfit Delta_r(a) = {c[2]:+.4f} {c[1]:+.4f} a {c[0]:+.4f} a^2",
          flush=True)
    print(f"  Delta_r_0 (Schw) = {c[2]:+.4f},  Delta_r_1 = {c[1]:+.4f}",
          flush=True)
    if c[1] != 0:
        print(f"  a_inv (lineare) ~ {-c[2]/c[1]:.4f}", flush=True)
    rr_ = np.roots(c)
    rr_ = rr_[(np.abs(rr_.imag) < 1e-9) & (rr_.real > 0)].real
    if len(rr_):
        print(f"  a_inv (quadratico) ~ {min(rr_):.4f}", flush=True)

# --------------------------------------------------------------- figura
if len(res) > 2:
    fig, ax = plt.subplots(figsize=(COL, COL * 0.75))
    ax.plot(res[:, 0], res[:, 3], 'C2o-', ms=4, label=r'$\Delta r$ (numeric)')
    ag = np.linspace(0, res[:, 0].max(), 100)
    ax.plot(ag, np.polyval(c, ag), 'k--', lw=0.8, label='quadratic fit')
    ax.axhline(0, color='k', lw=0.6)
    a_inv = -c[2] / c[1] if c[1] != 0 else np.nan
    if 0 < a_inv < res[:, 0].max():
        ax.axvline(a_inv, color='C3', ls=':', lw=1.0)
        ax.text(a_inv, 0.1, f'$a_{{inv}}={a_inv:.3f}$', fontsize=6, color='C3')
    ax.set_xlabel('$a$ (black-hole spin)')
    ax.set_ylabel(r'$\Delta r = r_{\min}^{t}-r_{\min}^{\tau}$')
    ax.set_title('Rotational inversion at fixed symmetric endpoints\n'
                 f'(Kerr, $r_0={r0}$, $\\pm{Phi}$, $E={E}$)')
    ax.legend(fontsize=6)
    savefig(fig, HERE, 'fig_bvp_kerr_inversione')
log("FATTO.")
