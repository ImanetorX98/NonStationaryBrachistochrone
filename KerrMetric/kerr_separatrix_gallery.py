# -*- coding: utf-8 -*-
"""
Galleria di separatrici: il test a 3 traiettorie (ODE, quadratura,
Weierstrass valutata) ripetuto su piu' parametri (a, E), per verificare
che la sovrapposizione non sia un caso speciale.

Ogni (a, E) ha il SUO reticolo: g2, g3, tau cambiano — e puo' cambiare
il TIPO (rombico Re tau = 1/2, discriminante < 0; rettangolare
Re tau = 0, discriminante > 0): la ricerca del reticolo prova entrambe
le famiglie e valida con l'uniformizzazione P(z(r)) = A/r + B.

Output: fig_separatrix_gallery + tabella deviazioni.
"""

import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import sympy as sp
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import DCOL, set_style, savefig
set_style()

mp.mp.dps = 25
M = 1.0
CASI = [(0.9, 1.2), (0.5, 1.2), (0.9, 1.05), (0.7, 1.5)]
r0 = 9.0

def caso(a, E):
    c = a / E
    Jc = a / E
    rp_ = M + np.sqrt(M**2 - a**2)
    rm_ = M - np.sqrt(M**2 - a**2)
    r_end = rp_ + 0.05
    Q4 = lambda r: r * (2 * M + (E**2 - 1) * r) * (r**2 + c**2)
    Dl = lambda r: r**2 - 2 * M * r + a**2
    dphidr = lambda r: (a / E) * np.sqrt(Q4(r)) / (Dl(r) * (r**2 + c**2))

    # T2: quadratura
    rg = np.linspace(r_end, r0, 240)
    phi2 = np.array([quad(dphidr, rv, r0, limit=300)[0] for rv in rg])

    # T1: ODE PMP
    r_s, pr_s, pphi_s = sp.symbols('r p_r p_phi', real=True)
    f_s = 1 - 2 * M / r_s
    Dl_s = r_s**2 - 2 * M * r_s + a**2
    P_s = r_s**2 + a**2 + 2 * M * a**2 / r_s
    vb2 = 1 - f_s / E**2
    Pb = P_s + (2 * M * a / r_s)**2 / E**2
    php0 = (2 * M * a / r_s) * vb2 / Pb
    R2s = Dl_s * vb2 / Pb
    ptp = pphi_s - 2 * M * a / (r_s * E)
    Hg_ = ptp * php0 + sp.sqrt(R2s) * sp.sqrt((Dl_s / r_s**2) * pr_s**2
                                              + ptp**2 / Pb) - f_s / E
    Hn = Hg_.subs(pphi_s, Jc)
    Hf = sp.lambdify((r_s, pr_s), Hn, 'numpy')
    dHp = sp.lambdify((r_s, pr_s), sp.diff(Hn, pr_s), 'numpy')
    dHr = sp.lambdify((r_s, pr_s), sp.diff(Hn, r_s), 'numpy')
    dHJ = sp.lambdify((r_s, pr_s),
                      sp.diff(Hg_, pphi_s).subs(pphi_s, Jc), 'numpy')
    with np.errstate(invalid='ignore'):
        pg = np.linspace(-200, 200, 200001)
        Hv = Hf(r0, pg)
    roots = [brentq(lambda p: Hf(r0, p), pg[i], pg[i + 1])
             for i in range(len(pg) - 1)
             if np.isfinite(Hv[i]) and np.isfinite(Hv[i + 1])
             and Hv[i] * Hv[i + 1] <= 0]
    p0 = [p_ for p_ in roots if dHp(r0, p_) < 0][0]
    ev = lambda t_, y: y[1]
    ev.terminal, ev.direction = True, 1

    def rhs(t_, y):
        return [dHp(y[0], y[1]), -dHr(y[0], y[1]), dHJ(y[0], y[1])]
    s = solve_ivp(rhs, [0, 6000], [r0, p0, 0.0], rtol=1e-12, atol=1e-14,
                  method='DOP853', events=[ev], dense_output=True,
                  max_step=2.0)
    tA = np.linspace(0, s.t[-1], 400)
    rA, phiA = s.sol(tA)[0], s.sol(tA)[2]

    # T3: Weierstrass sul reticolo vero
    A_ = M * a**2 / (2 * E**2)
    B_ = a**2 * (E**2 - 1) / (12 * E**2)
    g2 = a**4 * (E**2 - 1)**2 / (12 * E**4) - a**2 * M**2 / E**2
    g3 = a**4 * (36 * E**2 * M**2 * (1 - E**2)
                 - a**2 * (E**2 - 1)**3) / (216 * E**6)
    Q4m = lambda r: r * (2 * M + (E**2 - 1) * r) * (r**2 + c**2)
    z_of_r = lambda rv: mp.quad(lambda t: 1 / mp.sqrt(Q4m(t)), [0, rv])
    r_neg = -2 * M / (E**2 - 1)
    w1 = mp.quad(lambda t: 1 / mp.sqrt(Q4m(t)), [0, mp.inf]) \
        + mp.quad(lambda t: 1 / mp.sqrt(abs(Q4m(t))), [-mp.inf, r_neg])

    def make_wz(re_tau, y_):
        tau = mp.mpc(re_tau, y_)
        q = mp.exp(1j * mp.pi * tau)
        th1 = lambda v: mp.jtheta(1, v, q)
        th1p = lambda v: mp.jtheta(1, v, q, 1)
        eta1 = -(mp.pi**2 / (12 * w1)) * mp.jtheta(1, 0, q, 3) / th1p(0)

        def zeta_(z):
            v = mp.pi * z / (2 * w1)
            return eta1 * z / w1 + (mp.pi / (2 * w1)) * th1p(v) / th1(v)

        def sigma_(z):
            v = mp.pi * z / (2 * w1)
            return (2 * w1 / mp.pi) * mp.exp(eta1 * z**2 / (2 * w1)) \
                * th1(v) / th1p(0)
        wp_ = lambda z: -mp.diff(zeta_, z)
        return zeta_, sigma_, wp_

    z_t = mp.mpf(0.6) * w1

    def resid(re_tau, y_):
        _, _, wp_ = make_wz(re_tau, y_)
        p = wp_(z_t)
        return float(abs((mp.diff(wp_, z_t))**2
                         - (4 * p**3 - g2 * p - g3)))

    best = None
    for re_tau in (0.5, 0.0):
        ys = np.linspace(0.15, 2.2, 30)
        rs = [resid(re_tau, y_) for y_ in ys]
        y0 = ys[int(np.argmin(rs))]
        lo, hi = y0 - 0.08, y0 + 0.08
        for _ in range(55):
            m1, m2 = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
            if resid(re_tau, m1) < resid(re_tau, m2):
                hi = m2
            else:
                lo = m1
        y_o = 0.5 * (lo + hi)
        rr = resid(re_tau, y_o)
        if best is None or rr < best[2]:
            best = (re_tau, y_o, rr)
    re_t, y_o, rr = best
    zeta_w, sigma_w, wp_w = make_wz(re_t, y_o)
    err_p = float(abs(wp_w(z_of_r(5.0)) - (A_ / 5.0 + B_)))

    alpha, lamk, vk = {}, {}, {}
    for rk in (rp_, rm_):
        other = rm_ if rk == rp_ else rp_
        alpha[rk] = rk * ((E**2 - 1) * rk + 2 * M) / (rk - other)
        lamk[rk] = alpha[rk] / np.sqrt(Q4(rk))
        vk[rk] = z_of_r(rk)
    Lam0 = (E**2 - 1) - alpha[rp_] / rp_ - alpha[rm_] / rm_ \
        + 2 * float(mp.re(zeta_w(vk[rp_]))) * lamk[rp_] \
        + 2 * float(mp.re(zeta_w(vk[rm_]))) * lamk[rm_]

    def Phi(z):
        s_ = Lam0 * z
        for rk in (rp_, rm_):
            s_ += lamk[rk] * mp.log(sigma_w(z - vk[rk])
                                    / sigma_w(z + vk[rk]))
        return (a / E) * s_

    z0 = z_of_r(r0)
    phi3 = np.array([float(mp.re(Phi(z0) - Phi(z_of_r(rv))))
                     for rv in rg])

    devW = float(np.max(np.abs(phi3 - phi2)))
    devO = float(np.max(np.abs(phiA - np.interp(rA, rg, phi2))[rA > 2.2]))
    return dict(a=a, E=E, rg=rg, phi2=phi2, phi3=phi3, rA=rA, phiA=phiA,
                rp=rp_, r_graze=float(rA.min()), devW=devW, devO=devO,
                re_tau=re_t, y=y_o, res_lat=rr, err_p=err_p)

fig, axs = plt.subplots(2, 2, figsize=(DCOL, DCOL * 0.9))
print(f"{'a':>4} {'E':>5} | {'reticolo tau':>16} {'residuo':>9} "
      f"{'unif.':>9} | {'|W-quad|':>9} {'graze':>9} {'|ODE-q|':>9}")
for ax, (a, E) in zip(axs.flat, CASI):
    d = caso(a, E)
    print(f"{a:4.1f} {E:5.2f} | {d['re_tau']:.1f}+{d['y']:.4f}i "
          f"{d['res_lat']:9.1e} {d['err_p']:9.1e} | {d['devW']:9.1e} "
          f"{d['r_graze']:9.6f} {d['devO']:9.1e}", flush=True)
    ax.plot(d['rA'] * np.cos(d['phiA']), d['rA'] * np.sin(d['phiA']),
            'k-', lw=3.2, label='ODE')
    ax.plot(d['rg'] * np.cos(d['phi2']), d['rg'] * np.sin(d['phi2']),
            'C1--', lw=1.6, label='quadrature')
    ax.plot(d['rg'] * np.cos(d['phi3']), d['rg'] * np.sin(d['phi3']),
            'C2:', lw=1.8, label='Weierstrass')
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(2 * M * np.cos(th), 2 * M * np.sin(th), 'b--', lw=1.0)
    ax.plot(d['rp'] * np.cos(th), d['rp'] * np.sin(th), 'k-', lw=0.8)
    ax.set_aspect('equal')
    ax.set_title(f"$a={a}$, $E={E}$:  $|W-q|_{{max}}={d['devW']:.0e}$, "
                 f"graze $r={d['r_graze']:.4f}$", fontsize=9)
    ax.legend(fontsize=6.5, loc='lower left')
fig.suptitle('separatrices $J=J_c=a/E$: ODE vs quadrature vs evaluated '
             'Weierstrass — four parameters, four lattices')
HERE = os.path.dirname(os.path.abspath(__file__))
savefig(fig, HERE, 'fig_separatrix_gallery')
print('salvata fig_separatrix_gallery')
