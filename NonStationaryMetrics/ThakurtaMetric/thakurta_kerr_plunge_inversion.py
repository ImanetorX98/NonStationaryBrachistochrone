# -*- coding: utf-8 -*-
"""
Inversione del ramo tau e tricotomia all'ergosfera in Thakurta-Kerr.

NOTA: la tricotomia classifica il comportamento all'ERGOSFERA r_e = 2M
(superficie di luce): J>J_c periasse liscio; J<J_c cuspide su r_e;
J=J_c attraversa e spirala su r_+ (cattura, misura nulla).

FORMULA ANALITICA (derivazione): con A costante, g = A^2 g_Kerr e' Kerr
riscalato: il worldline di rotaia con Ehat = -u_eta equivale al problema
di Kerr con  E_eff = Ehat/A  e momento di Fermat  J_eff = J/A
(T_taubar = A T_tau => p_phi si riscala di A). Quindi il raggio di
inversione (periasse) del ramo tau e' la radice esterna di

    g(r) = Dl(r) - (J/A)^2 [ (Ehat/A)^2 - f(r) ]        (*)

e la linea critica di attraversamento (tricotomia doranTau, J_c = s/E):

    J_c(A) = s A^2 / Ehat                                (**)

Verifiche:
  V1  (*) coincide con la condizione p_r = 0, H_tau = 0 (simbolico-numerico)
  V2  COLORMAP statica (J, A): r_min dal flusso di Hamilton vs radice
      di (*): residuo relativo (atteso ~1e-9) + contorni analitici
      sovrapposti; linea (**) sovrapposta
  V3  dinamico (a = -1/(H eta)): r_min numerico vs (*) valutata ad
      A(eta_periasse): validita' quasi-statica e sua rottura
Figura: Thakurtafigures/fig_thakurta_kerr_plunge_map
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from paper_style import COL, set_style, savefig

set_style()
OUT = os.path.join(HERE, 'Thakurtafigures')

M_n, s_n, E_n = 1.0, 0.9, 1.2
r0 = 8.0
r_plus = M_n + np.sqrt(M_n**2 - s_n**2)

# --------------------------------------------------- H_tau (A parametro)
r, pr, A = sp.symbols('r p_r A', positive=True)
pphi = sp.Symbol('p_phi', real=True)
f = 1 - 2 * M_n / r
Dl = r**2 - 2 * M_n * r + s_n**2
P = r**2 + s_n**2 + 2 * M_n * s_n**2 / r
vb2 = 1 - A**2 * f / E_n**2
Pb = P + A**2 * (2 * M_n * s_n / r)**2 / E_n**2
php0 = (2 * M_n * s_n / r) * vb2 / Pb
R_ = sp.sqrt(Dl * vb2 / Pb)
ptp = pphi - 2 * M_n * s_n * A**2 / (r * E_n)
H_tau = ptp * php0 + R_ * sp.sqrt((Dl / r**2) * pr**2 + ptp**2 / Pb) \
    - A**2 * f / E_n

Hf = sp.lambdify((r, pr, A, pphi), H_tau, 'numpy')
dHp = sp.lambdify((r, pr, A, pphi), sp.diff(H_tau, pr), 'numpy')
dHr = sp.lambdify((r, pr, A, pphi), sp.diff(H_tau, r), 'numpy')

def g_inv(rv, Jv, Av):
    """funzione di svolta analitica (*)"""
    fv = 1 - 2 * M_n / rv
    return (rv**2 - 2 * M_n * rv + s_n**2) \
        - (Jv / Av)**2 * ((E_n / Av)**2 - fv)

def r_inv_analitico(Jv, Av, con_cuspide=True):
    """r_min analitico: radice esterna di (*); tricotomia doranTau
    riscalata: J < J_c(A)=sA^2/Ehat => cuspide a r_e=2M; J=J_c attraversa."""
    Jc = s_n * Av**2 / E_n
    if con_cuspide and Jv < Jc:
        return 2 * M_n                    # riflessione a r_e (cuspide)
    rg = np.linspace(r0, 2 * M_n * 0.999, 4000)
    gv = g_inv(rg, Jv, Av)
    for i in range(len(rg) - 1):
        if gv[i] > 0 >= gv[i + 1]:
            rt = brentq(lambda x: g_inv(x, Jv, Av), rg[i + 1], rg[i])
            return max(rt, 2 * M_n) if con_cuspide else rt
    return np.nan          # nessuna svolta esterna: attraversa

print("=" * 72)
print("[V1] la condizione p_r=0, H_tau=0 e' la (*)")
print("=" * 72)
rng = np.random.default_rng(3)
ok = True
for _ in range(12):
    Av = rng.uniform(0.7, 1.5)
    Jv = rng.uniform(1.0, 3.0)
    if Jv <= s_n * Av**2 / E_n * 1.02:    # sotto/vicino a J_c: cuspide
        continue
    ra = r_inv_analitico(Jv, Av, con_cuspide=False)
    if not np.isfinite(ra):
        continue
    h0 = Hf(ra, 0.0, Av, Jv)
    ok = ok and abs(h0) < 1e-10
print("  H_tau(p_r=0, r=r_inv(*)) = 0  (J > J_c, svolta liscia):", ok)
print("  J < J_c(A) = sA^2/Ehat: cuspide a r_e = 2M (tricotomia riscalata)")

# ------------------------------------------------- V2: colormap statica
print()
print("=" * 72)
print("[V2] colormap statica (J, A): flusso numerico vs (*)")
print("=" * 72)

def r_min_numerico(Jv, Av):
    """lancio entrante da r0 con H=0 (A costante): svolta o plunge."""
    pg = np.linspace(-150, 150, 300001)
    with np.errstate(invalid='ignore'):
        Hg = Hf(r0, pg, Av, Jv)
    roots = [brentq(lambda p: Hf(r0, p, Av, Jv), pg[i], pg[i + 1])
             for i in range(len(pg) - 1)
             if np.isfinite(Hg[i]) and np.isfinite(Hg[i + 1])
             and Hg[i] * Hg[i + 1] <= 0]
    cand = [p_ for p_ in roots if dHp(r0, p_, Av, Jv) < 0]
    if not cand:
        return np.nan
    p0 = cand[-1]

    def rhs(e_, y):
        return [dHp(y[0], y[1], Av, Jv), -dHr(y[0], y[1], Av, Jv)]
    ev_t = lambda e_, y: dHp(y[0], y[1], Av, Jv)
    ev_t.terminal, ev_t.direction = True, 1
    ev_h = lambda e_, y: y[0] - r_plus * 1.005
    ev_h.terminal, ev_h.direction = True, -1
    s = solve_ivp(rhs, [0, 4000], [r0, p0], rtol=1e-11, atol=1e-13,
                  method='DOP853', events=[ev_t, ev_h])
    if len(s.t_events[0]):
        return s.y_events[0][0][0]
    return 0.0             # attraversa l'ergosfera (cattura)

J_g = np.linspace(0.3, 3.2, 30)     # J basso: la separatrice J_c entra in frame
A_g = np.linspace(0.7, 1.5, 25)
num = np.zeros((len(A_g), len(J_g)))
ana = np.zeros_like(num)
for i, Av in enumerate(A_g):
    for j, Jv in enumerate(J_g):
        num[i, j] = r_min_numerico(Jv, Av)
        ana[i, j] = r_inv_analitico(Jv, Av)
mask = np.isfinite(ana) & (num > 0)
res = np.where(mask, np.abs(num - ana) / np.where(ana > 0, ana, 1), np.nan)
print(f"  residuo relativo |r_min^num - r_inv^(*)|: max = "
      f"{np.nanmax(res):.2e}, mediana = {np.nanmedian(res):.2e}")

# ----------------------------------------------------- V3: quasi-statico
print()
print("=" * 72)
print("[V3] dinamico a = -1/(H eta): quasi-statico vs numerico")
print("=" * 72)
H_c = 0.02
eta0 = -1.0 / H_c

def a_dyn(e_):
    return -1.0 / (H_c * e_)

A_sym = sp.Symbol('A', positive=True)
ee = sp.Symbol('ee', negative=True)
H_dyn = H_tau.subs(A, -1 / (H_c * ee))
Hf_d = sp.lambdify((ee, r, pr, pphi), H_dyn, 'numpy')
dHp_d = sp.lambdify((ee, r, pr, pphi), sp.diff(H_dyn, pr), 'numpy')
dHr_d = sp.lambdify((ee, r, pr, pphi), sp.diff(H_dyn, r), 'numpy')

def r_min_dinamico(Jv, eta_l):
    pg = np.linspace(-150, 150, 300001)
    with np.errstate(invalid='ignore'):
        Hg = Hf_d(eta_l, r0, pg, Jv)
    roots = [brentq(lambda p: Hf_d(eta_l, r0, p, Jv), pg[i], pg[i + 1])
             for i in range(len(pg) - 1)
             if np.isfinite(Hg[i]) and np.isfinite(Hg[i + 1])
             and Hg[i] * Hg[i + 1] <= 0]
    cand = [p_ for p_ in roots if dHp_d(eta_l, r0, p_, Jv) < 0]
    if not cand:
        return np.nan, np.nan
    p0 = cand[-1]

    def rhs(e_, y):
        return [dHp_d(e_, y[0], y[1], Jv), -dHr_d(e_, y[0], y[1], Jv)]
    ev_t = lambda e_, y: dHp_d(e_, y[0], y[1], Jv)
    ev_t.terminal, ev_t.direction = True, 1
    ev_h = lambda e_, y: y[0] - r_plus * 1.005
    ev_h.terminal, ev_h.direction = True, -1
    s = solve_ivp(rhs, [eta_l, -1e-2], [r0, p0], rtol=1e-11, atol=1e-13,
                  method='DOP853', events=[ev_t, ev_h])
    if len(s.t_events[0]):
        return s.y_events[0][0][0], s.t_events[0][0]
    if len(s.t_events[1]):
        return 0.0, s.t_events[1][0]
    return np.nan, np.nan

J_g3 = np.linspace(0.3, 3.2, 30)
# A0 fino a 1.3 (oltre ~1.39 il punto di lancio r0 e' gia' congelato)
eta_g = np.linspace(eta0, eta0 / 1.3, 25)
num3 = np.zeros((len(eta_g), len(J_g3)))
qs3 = np.zeros_like(num3)
for i, el in enumerate(eta_g):
    for j, Jv in enumerate(J_g3):
        rm, e_ev = r_min_dinamico(Jv, el)
        num3[i, j] = rm
        qs3[i, j] = r_inv_analitico(Jv, a_dyn(e_ev)) \
            if np.isfinite(e_ev) and rm > 0 else np.nan
mask3 = np.isfinite(qs3) & (num3 > 0)
shift3 = np.where(mask3, (num3 - qs3) / np.where(qs3 > 0, qs3, 1), np.nan)
print(f"  H = {H_c}: SHIFT non adiabatico (r_num - r_qs)/r_qs:")
print(f"    mediana = {np.nanmedian(shift3):+.3f}, max = "
      f"{np.nanmax(shift3):+.3f}, min = {np.nanmin(shift3):+.3f}")
print("    (positivo: l'inversione e' ANTICIPATA — il congelamento in")
print("     avvicinamento decelera la particella prima della svolta")
print("     quasi-statica: effetto genuinamente dinamico)")
n_plunge = int(np.sum(num3 == 0.0))
print(f"  celle catturate (attraversamento ergosfera): {n_plunge}/{num3.size}")

# ------------------------------------------------------------- figura
fig, (axa, axb) = plt.subplots(2, 1, figsize=(COL, 6.6))
ext = [J_g[0], J_g[-1], A_g[0], A_g[-1]]
im = axa.imshow(num, origin='lower', extent=ext, aspect='auto',
                cmap='viridis', vmin=2.0)
plt.colorbar(im, ax=axa, label='$r_{\\min}$ (num. = analytic)')
CS = axa.contour(J_g, A_g, ana, levels=[2.5, 3, 4, 5, 6],
                 colors='w', linewidths=0.8)
axa.clabel(CS, fontsize=6, fmt='%.1f')
Ac = np.linspace(A_g[0], A_g[-1], 200)
axa.plot(s_n * Ac**2 / E_n, Ac, 'r--', lw=2.0,
         label='separatrix $J_c=sA^2/\\hat E$')
axa.text(0.42, 1.32, 'CAPTURE\n$r_{\\min}=2M$', fontsize=6.5,
         color='darkorange', ha='center', va='center', fontweight='bold')
axa.text(2.3, 0.85, 'SCATTER\n(smooth periapsis)', fontsize=6.5,
         color='k', ha='center', va='center')
axa.set_xlabel('$J$ (particle angular momentum)')
axa.set_ylabel('$A$ (conformal factor)')
axa.set_xlim(J_g[0], J_g[-1])
axa.set_ylim(A_g[0], A_g[-1])
axa.set_title('Static: $r_{\\min}$ colormap $=$ analytic contours '
              f'(residual {np.nanmax(res):.0e});\nseparatrix $J_c$ '
              'divides capture ($r_{\\min}{=}2M$) from scatter')
axa.legend(loc='lower right', fontsize=6)

A_ax = np.array([a_dyn(e_) for e_ in eta_g])
im2 = axb.pcolormesh(J_g3, A_ax, 100 * shift3, cmap='magma',
                     shading='auto')
plt.colorbar(im2, ax=axb,
             label='non-adiab. shift $(r_{\\min}-r_{qs})/r_{qs}$ [%]')
axb.plot(s_n * A_ax**2 / E_n, A_ax, 'c--', lw=2.0,
         label='separatrix $J_c(A_0)$')
axb.set_xlabel('$J$ (particle angular momentum)')
axb.set_ylabel('$A_0$ at launch')
axb.set_xlim(J_g3[0], J_g3[-1])
axb.set_title(f'Dynamic ($a=-1/H\\eta$, $H={H_c}$): freezing decel.\n'
              'ANTICIPATES the turn vs quasi-static (median '
              f'{np.nanmedian(shift3):+.0%}, brighter = larger)')
axb.legend(loc='upper right', fontsize=6)
savefig(fig, OUT, 'fig_thakurta_kerr_plunge_map')
print('\nFATTO.')
