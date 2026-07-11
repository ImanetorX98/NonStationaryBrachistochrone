# -*- coding: utf-8 -*-
"""
Figure di validazione — brachistocrone in Vaidya (accrescimento lineare).

Output in Vaidyafigures/ (png+pdf):
  fig_vaidya_orizzonti     (a) EH vs AH per accrescimento con stop
                           (b) raggi autosimilari x(v) e x_pm per m = mu*v
  fig_vaidya_kerr_a0       residui |dphi/dr - forma chiusa Kerr a=0| (m'=0)
  fig_vaidya_traiettorie   orbite tau Schwarzschild vs Vaidya + p_v(r)
  fig_vaidya_variazionale  T_tau su famiglia perturbata a Delta_phi fisso:
                           minimo sulla soluzione EL (validazione diretta)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig

set_style()
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Vaidyafigures')

def salva(fig, nome):
    savefig(fig, OUT, nome)

# ------------------------------------------------------------------ sympy
r, vv, uu, q = sp.symbols('r vv uu q', real=True)
E, J, mm, mp, mu = sp.symbols('E J mm mp mu', positive=True)
f = 1 - 2 * mm / r
W = (f * uu**2 - 2 * uu - (f * uu - 1)**2 / E**2) / r**2
phip = sp.sqrt(W)
F_tau = (f * uu - 1) / E - J * phip

p_tau = sp.diff(F_tau, uu)
dp_tot = (sp.diff(p_tau, r) + sp.diff(p_tau, vv) * uu
          + sp.diff(p_tau, mm) * mp * uu + sp.diff(p_tau, uu) * q)
EL = dp_tot - sp.diff(F_tau, mm) * mp
vpp_tau = sp.solve(sp.Eq(EL, 0), q)[0]

m0v, E_n, J_n = 1.0, 1.2, 1.3
r1, r0 = 10.0, 3.5     # periasse (E=1.2, J=1.3) a r~3.01: restare sopra
lin = [(mm, m0v + mu * vv), (mp, mu), (E, E_n), (J, J_n)]
args = (r, vv, uu, mu)
vpp_fn = sp.lambdify(args, vpp_tau.subs(lin), 'numpy')
pv_fn = sp.lambdify(args, p_tau.subs(lin), 'numpy')
W_fn = sp.lambdify(args, W.subs(lin), 'numpy')
Lam_fn = sp.lambdify(args, ((f * uu - 1) / E).subs(lin), 'numpy')

def vp_radiale(f_val, E_val):
    a2 = f_val - f_val**2 / E_val**2
    a1 = -2 + 2 * f_val / E_val**2
    a0 = -1 / E_val**2
    return (-a1 + np.sqrt(a1**2 - 4 * a2 * a0)) / (2 * a2)

def integra(mu_n, v_arrivo, dense=False):
    fl = 1 - 2 * (m0v + mu_n * v_arrivo) / r1
    lo = vp_radiale(fl, E_n) * (1 + 1e-9)
    hi = lo * 2
    while pv_fn(r1, v_arrivo, hi, mu_n) < 0:
        hi *= 2
    up1 = brentq(lambda z: pv_fn(r1, v_arrivo, z, mu_n), lo, hi)

    def rhs(rr, y):
        return [y[1], vpp_fn(rr, y[0], y[1], mu_n),
                np.sqrt(max(W_fn(rr, y[0], y[1], mu_n), 0.0)),
                Lam_fn(rr, y[0], y[1], mu_n)]

    sol = solve_ivp(rhs, [r1, r0], [v_arrivo, up1, 0.0, 0.0],
                    rtol=1e-11, atol=1e-13, dense_output=True)
    assert sol.status == 0, f'integrazione fallita a r={sol.t[-1]}'
    return sol

# ---------------------------------------------------------------- fig 1
print('fig 1: orizzonti...')
mu_a, v_stop = 0.02, 50.0

def mass(vq):
    return m0v + mu_a * min(vq, v_stop)

def ray(rq0, v_end=2000.0):
    s = solve_ivp(lambda vq, y: [0.5 * (1 - 2 * mass(vq) / y[0])],
                  [0, v_end], [rq0], rtol=1e-10, atol=1e-12,
                  dense_output=True, max_step=1.0)
    return s

lo, hi = 2.0, 6.0
for _ in range(55):
    mid = 0.5 * (lo + hi)
    if ray(mid).y[0, -1] > 10 * mass(1e9):
        hi = mid
    else:
        lo = mid
sol_eh = ray(0.5 * (lo + hi))
vg = np.linspace(0, 100, 400)
r_eh = sol_eh.sol(vg)[0]
r_ah = np.array([2 * mass(x) for x in vg])

fig, (axa, axb) = plt.subplots(2, 1, figsize=(COL, 5.4))
axa.plot(vg, r_eh, 'C0-', label='event horizon (critical generator)')
axa.plot(vg, r_ah, 'C3--', label='apparent horizon $r=2m(v)$')
axa.axvline(v_stop, color='k', ls=':', lw=0.8)
axa.annotate('accretion stops', (v_stop + 1, 2.3), fontsize=6.5)
axa.set_xlabel('$v$')
axa.set_ylabel('$r$')
axa.set_title(f'$m=1+{mu_a}v$ up to $v={v_stop:.0f}$: EH outside AH,\n'
              'then they merge (teleological EH)')
axa.legend()

mu_n = 0.02
disc = np.sqrt(1 - 16 * mu_n)
xm, xp = (1 - disc) / 4, (1 + disc) / 4
vgl = np.logspace(0, 5, 500)
for x0 in (0.030, 0.040, 0.0438, 0.045, 0.10, 0.30, 0.60):
    s = solve_ivp(lambda vq, y: [((1 - 2 * mu_n / y[0]) / 2 - y[0]) / vq],
                  [1, 1e5], [x0], rtol=1e-11, atol=1e-14, dense_output=True)
    xv = s.sol(np.clip(vgl, 1, s.t[-1]))[0]
    axb.plot(vgl, xv, lw=0.9,
             color='C3' if x0 < xm else 'C0')
axb.axhline(xm, color='k', lw=1.2)
axb.axhline(xp, color='k', lw=1.2, ls='--')
axb.axhline(2 * mu_n, color='C2', lw=1.0, ls=':')
axb.text(1e3, xm * 1.06, '$x_-$ (self-similar EH)', fontsize=6.5)
axb.text(1e3, xp * 1.05, '$x_+$', fontsize=6.5)
axb.text(1e3, 2 * mu_n * 0.86, '$x_{AH}=2\\mu$', fontsize=6.5, color='C2')
axb.set_xscale('log')
axb.set_yscale('log')
axb.set_ylim(0.015, 0.8)
axb.set_xlabel('$v$')
axb.set_ylabel('$x = r/v$')
axb.set_title(f'$m=\\mu v$, $\\mu={mu_n}$: outgoing rays;\n'
              '$x_\\pm=(1\\pm\\sqrt{1-16\\mu})/4$, critical $\\mu_c=1/16$')
salva(fig, 'fig_vaidya_orizzonti')

# ---------------------------------------------------------------- fig 2
print('fig 2: residui forma chiusa Kerr a=0...')
F_t = uu - J * phip
w_ = E**2 - f
Dl = r**2 * f
fig, ax = plt.subplots(figsize=(COL, COL * 0.82))
rgrid = np.linspace(3.0, 15.0, 40)
for (nome, F, K), st in ((('$\\tau$ branch', F_tau, J), 'o-'),
                         (('$t$ branch', F_t, f * J / E), 's--')):
    pv = sp.diff(F, uu)
    tgt = K * sp.sqrt(w_ * f) * r / (Dl * sp.sqrt(Dl - K**2 * w_))
    for Jv, col in ((0.8, 'C0'), (1.3, 'C3')):
        res = []
        for rv in rgrid:
            sub = {mm: m0v, E: E_n, J: Jv, r: rv}
            pv_l = sp.lambdify(uu, pv.subs(sub), 'numpy')
            lo = vp_radiale(float(f.subs(sub)), E_n) * (1 + 1e-9)
            hi = lo * 2
            while pv_l(hi) < 0:
                hi *= 2
            ustar = brentq(pv_l, lo, hi, xtol=1e-15, rtol=8.9e-16)
            num = float(phip.subs(sub).subs(uu, ustar))
            cf = float(tgt.subs(sub))
            res.append(abs(num - cf) / cf)
        ax.semilogy(rgrid, np.maximum(res, 1e-17), st, color=col, ms=3,
                    lw=0.9, label=f'{nome}, $J={Jv}$')
ax.set_xlabel('$r/M$')
ax.set_ylabel('relative residual $|\\Delta(d\\phi/dr)|$')
ax.set_title('stationary-limit check ($m\'=0$, $p_v=0$):\n'
             'numerical orbit vs closed-form Kerr $a=0$ ($E=1.2$)')
ax.legend()
salva(fig, 'fig_vaidya_kerr_a0')

# ---------------------------------------------------------------- fig 3
print('fig 3: traiettorie...')
v_arr = 40.0
sol_s = integra(0.0, v_arr)
sol_v = integra(0.01, v_arr)

fig, (axo, axp) = plt.subplots(2, 1, figsize=(COL, 5.8))
for sol, mu_c, col, lab in ((sol_s, 0.0, 'C0', 'Schwarzschild $m=1$'),
                            (sol_v, 0.01, 'C3', 'Vaidya $m=1+0.01v$')):
    rr = np.linspace(r1, r0, 800)
    Y = sol.sol(rr)
    phi = Y[2] - Y[2, -1]           # phi=0 alla partenza
    axo.plot(rr * np.cos(phi), rr * np.sin(phi), color=col, label=lab)
    m_dep = m0v + mu_c * Y[0, -1]
    m_arr = m0v + mu_c * v_arr
    th = np.linspace(0, 2 * np.pi, 200)
    axo.plot(2 * m_dep * np.cos(th), 2 * m_dep * np.sin(th), color=col,
             ls=':', lw=0.9)
    axo.plot(2 * m_arr * np.cos(th), 2 * m_arr * np.sin(th), color=col,
             ls='--', lw=0.9)
    pv_line = [pv_fn(x, sol.sol(x)[0], sol.sol(x)[1], mu_c) for x in rr]
    axp.plot(rr, pv_line, color=col, label=lab)
axo.plot([r0], [0], 'k^', ms=7, label=f'start $r_0={r0}$')
axo.set_aspect('equal')
axo.set_xlabel('$x$')
axo.set_ylabel('$y$')
axo.set_title('$\\tau$-brachistochrones ($E=1.2$, $J=1.3$, arrival $v=40$)\n'
              'AH at start (dotted) and arrival (dashed)')
axo.legend(loc='lower right')
axp.axhline(0, color='k', lw=0.6)
axp.set_xlabel('$r$')
axp.set_ylabel('$p_v(r)$')
axp.set_title('transversality $p_v(r_1)=0$; in Vaidya $p_v\\neq 0$ along\n'
              'the path: non-autonomous memory ($\\propto m\'$)')
axp.legend()
salva(fig, 'fig_vaidya_traiettorie')

# ---------------------------------------------------------------- fig 4
print('fig 4: validazione variazionale diretta...')
mu_c = 0.01
sol = sol_v
rr = np.linspace(r0, r1, 1500)
v_el = CubicSpline(rr, sol.sol(rr)[0])
dphi_target = float(sol.sol(r0)[2] * -1)   # Delta_phi della soluzione EL
T_el = float(-sol.sol(r0)[3])

psi1 = lambda x: np.sin(np.pi * (x - r0) / (r1 - r0))     # 0 agli estremi
psi2 = lambda x: (r1 - x) / (r1 - r0)                     # libera in partenza
dpsi1 = lambda x: np.pi / (r1 - r0) * np.cos(np.pi * (x - r0) / (r1 - r0))
dpsi2 = lambda x: -1.0 / (r1 - r0) + 0 * x

def funzionali(e1, e2):
    vq = v_el(rr) + e1 * psi1(rr) + e2 * psi2(rr)
    uq = v_el(rr, 1) + e1 * dpsi1(rr) + e2 * dpsi2(rr)
    Wq = W_fn(rr, vq, uq, mu_c)
    if np.any(Wq < 0):
        return None, None
    dphi = np.trapezoid(np.sqrt(Wq), rr)
    Tq = np.trapezoid(Lam_fn(rr, vq, uq, mu_c), rr)
    return Tq, dphi

def risolvi_e2(e1):
    """Trova e2 con Delta_phi = target: scansione + brentq locale."""
    def gap(e2):
        _, d = funzionali(e1, e2)
        return np.nan if d is None else d - dphi_target
    grid = np.linspace(-1.5, 1.5, 301)
    gv = np.array([gap(x) for x in grid])
    for i in range(len(grid) - 1):
        if np.isfinite(gv[i]) and np.isfinite(gv[i + 1]) \
                and gv[i] * gv[i + 1] <= 0:
            return brentq(gap, grid[i], grid[i + 1], xtol=1e-12)
    return None

eps_grid = np.linspace(-0.12, 0.05, 35)
Tt = []
for e1 in eps_grid:
    e2s = risolvi_e2(e1)
    Tq = np.nan if e2s is None else funzionali(e1, e2s)[0]
    Tt.append(Tq)
Tt = np.array(Tt)

fig, ax = plt.subplots(figsize=(COL, COL * 0.78))
ax.plot(eps_grid, Tt, 'C3.-', lw=0.9)
ax.axhline(T_el, color='k', lw=0.7, ls='--',
           label=f'EL solution: $T_\\tau={T_el:.6f}$')
ax.axvline(0, color='k', lw=0.6)
ax.set_xlabel(r'perturbation $\varepsilon_1$ on $v(r)$ '
              r'($\Delta\phi$ fixed via $\varepsilon_2$)')
ax.set_ylabel(r'$T_\tau$')
ax.set_title('variational check: the non-autonomous EL solution\n'
             f'minimizes $T_\\tau$ at fixed $\\Delta\\phi$ (Vaidya $\\mu={mu_c}$)')
ax.legend()
salva(fig, 'fig_vaidya_variazionale')
imin = np.nanargmin(Tt)
print(f'  [validazione] argmin T_tau a eps1 = {eps_grid[imin]:.4f} '
      f'(T_min={np.nanmin(Tt):.8f} vs T_EL={T_el:.8f})')
print('FATTO.')
