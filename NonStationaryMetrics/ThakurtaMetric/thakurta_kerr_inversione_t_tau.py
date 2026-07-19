# -*- coding: utf-8 -*-
"""
############################################################################
## DEPRECATO per la FIGURA. Il piano (J,A) qui usa J in [3.9,5.5], range  ##
## TROPPO STRETTO: nasconde l'inversione conformale, che vive a J~16-24.  ##
## La conclusione "0/484 celle, inversione non accessibile" era ERRATA.   ##
## Figura corretta: inversione_conformale_AJ.py ->                        ##
##   fig_thakurta_kerr_inversione_AJ. Teorema in ThakurtaResults.md R10b   ##
##   (cond(r,A)=0; A_inv in (1, A_freeze), sempre raggiungibile).          ##
## Lo script resta valido come cross-check ρ(r*)²=H_t0/H_τ0 (R10).         ##
############################################################################

Inversione di plunge t vs tau in Thakurta-Kerr (equatoriale).

Riferimento Kerr (evidence_t_plunges_deeper.txt): a parita' di lancio il
ramo t raggiunge r_min MINORE del ramo tau (t affonda di piu', 58/60),
perche' il potenziale centrifugo efficace del ramo t e' soppresso:
V_t/V_tau = [F(J - A_phi)/(r^2 E J)]^2 = (K/J)^2, con K il momento
efficace di doranT; a grande r V_t ~ V_tau/E^2.

THAKURTA-KERR (g = A^2 g_Kerr, rotaia Ehat): il rapporto acquista il
fattore conforme. Nella famiglia brachistocrona (H = 0), le funzioni di
svolta sono

    g_tau(r) = Dl - (J/A)^2  w_A          w_A = (Ehat/A)^2 - f
    g_t(r)   = Dl - K_A^2    w_A          K_A = A (f J + 2Ms/r)/Ehat

e il RAPPORTO DELLE HAMILTONIANE (potenziali centrifughi):

    rho(r) = K_A/(J/A) = A^2 [ f + 2Ms/(rJ) ] / Ehat        (#)

    rho < 1  =>  t affonda di piu' (Kerr: rho ~ f/E < 1)
    rho = 1  =>  INVERSIONE:  r_min^t = r_min^tau
    rho > 1  =>  tau affonda di piu'

Il fattore A^2 in (#) e' il risultato nuovo: l'espansione conforme
ALZA il rapporto — l'inversione, marginale in Kerr, diventa una CURVA
raggiungibile A_inv(J) dentro la finestra accessibile (prima del
congelamento al lancio A < Ehat/sqrt(f(r0))).

Curva analitica di inversione: sistema { g_tau(r)=0 , rho(r)=1 }
(sul luogo di inversione le due radici coincidono).

Verifiche:
  V1  H_eta(p_r=0, radice di g_t) = 0 e H_tau(p_r=0, radice di g_tau)=0
  V2  colormap (J, A) di Delta r = r_min^t - r_min^tau dai flussi di
      Hamilton (famiglia H=0, lancio entrante da r0), con:
      - residui numerico vs radici analitiche (attesi ~1e-10);
      - contorno Delta r = 0 numerico vs curva analitica (#)=1.
Figura: Thakurtafigures/fig_thakurta_kerr_inversione_t_tau
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
r0 = 6.0
r_plus = M_n + np.sqrt(M_n**2 - s_n**2)
A_freeze = E_n / np.sqrt(1 - 2 * M_n / r0)      # lancio congelato oltre

# --------------------------------------------------------- Hamiltoniane
r, pr, A = sp.symbols('r p_r A', positive=True)
pphi = sp.Symbol('p_phi', real=True)
f = 1 - 2 * M_n / r
Dl = r**2 - 2 * M_n * r + s_n**2
P = r**2 + s_n**2 + 2 * M_n * s_n**2 / r
vb2 = 1 - A**2 * f / E_n**2
Pb = P + A**2 * (2 * M_n * s_n / r)**2 / E_n**2
php0 = (2 * M_n * s_n / r) * vb2 / Pb
R_ = sp.sqrt(Dl * vb2 / Pb)
H_eta = pphi * php0 + R_ * sp.sqrt((Dl / r**2) * pr**2 + pphi**2 / Pb) - 1
ptp = pphi - 2 * M_n * s_n * A**2 / (r * E_n)
H_tau = ptp * php0 + R_ * sp.sqrt((Dl / r**2) * pr**2 + ptp**2 / Pb) \
    - A**2 * f / E_n

lam = lambda ex: sp.lambdify((r, pr, A, pphi), ex, 'numpy')
Hs = {'tau': (lam(H_tau), lam(sp.diff(H_tau, pr)), lam(sp.diff(H_tau, r))),
      'eta': (lam(H_eta), lam(sp.diff(H_eta, pr)), lam(sp.diff(H_eta, r)))}

# ------------------------------------------------ svolte analitiche
def f_n(rv):
    return 1 - 2 * M_n / rv

def w_A(rv, Av):
    return (E_n / Av)**2 - f_n(rv)

def g_branch(rv, Jv, Av, ramo):
    Dv = rv**2 - 2 * M_n * rv + s_n**2
    if ramo == 'tau':
        K = Jv / Av
    else:
        K = Av * (f_n(rv) * Jv + 2 * M_n * s_n / rv) / E_n
    return Dv - K**2 * w_A(rv, Av)

def r_turn_analitico(Jv, Av, ramo):
    rg = np.linspace(r0, r_plus * 1.001, 3000)
    gv = g_branch(rg, Jv, Av, ramo)
    for i in range(len(rg) - 1):
        if gv[i] > 0 >= gv[i + 1]:
            return brentq(lambda x: g_branch(x, Jv, Av, ramo),
                          rg[i + 1], rg[i])
    return np.nan

print("=" * 72)
print("[V1] radici analitiche = condizione p_r=0, H=0")
print("=" * 72)
rng = np.random.default_rng(4)
ok = True
for _ in range(10):
    Av = rng.uniform(0.9, 1.35)
    Jv = rng.uniform(4.0, 5.5)
    for ramo, Hx in (('tau', Hs['tau'][0]), ('eta', Hs['eta'][0])):
        ra = r_turn_analitico(Jv, Av, ramo)
        if np.isfinite(ra):
            ok = ok and abs(Hx(ra, 0.0, Av, Jv)) < 1e-9
print("  H(p_r=0, r_turn) = 0 per entrambi i rami:", ok)

def rho(rv, Jv, Av):
    return Av**2 * (f_n(rv) + 2 * M_n * s_n / (rv * Jv)) / E_n

def curva_inversione(Av):
    """J_inv(A): sistema g_tau(r)=0 e rho(r)=1."""
    def eq(Jv):
        rt = r_turn_analitico(Jv, Av, 'tau')
        if not np.isfinite(rt):
            return np.nan
        return rho(rt, Jv, Av) - 1.0
    Jg = np.linspace(3.8, 5.6, 60)
    vals = np.array([eq(Jv) for Jv in Jg])
    for i in range(len(Jg) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) \
                and vals[i] * vals[i + 1] <= 0:
            return brentq(eq, Jg[i], Jg[i + 1])
    return np.nan

# --------------------------------------------------- flussi numerici
def r_min_num(Jv, Av, ramo):
    Hx, dHp, dHr = Hs[ramo]
    pg = np.linspace(-200, 200, 200001)
    with np.errstate(invalid='ignore'):
        Hg = Hx(r0, pg, Av, Jv)
    roots = [brentq(lambda p: Hx(r0, p, Av, Jv), pg[i], pg[i + 1])
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
    return 0.0

print()
print("=" * 72)
print(f"[V2] colormap (J, A): Delta r = r_min^t - r_min^tau  "
      f"(r0={r0}, A_freeze={A_freeze:.3f})")
print("=" * 72)
J_g = np.linspace(3.9, 5.5, 23)
A_g = np.linspace(0.9, min(1.42, A_freeze * 0.985), 23)
d_num = np.full((len(A_g), len(J_g)), np.nan)
res_max = 0.0
for i, Av in enumerate(A_g):
    for j, Jv in enumerate(J_g):
        rt_n = r_min_num(Jv, Av, 'tau')
        re_n = r_min_num(Jv, Av, 'eta')
        rt_a = r_turn_analitico(Jv, Av, 'tau')
        re_a = r_turn_analitico(Jv, Av, 'eta')
        if np.isfinite(rt_n) and np.isfinite(re_n) and rt_n > 0 \
                and re_n > 0 and np.isfinite(rt_a) and np.isfinite(re_a):
            d_num[i, j] = re_n - rt_n
            res_max = max(res_max, abs(rt_n - rt_a) / rt_a,
                          abs(re_n - re_a) / re_a)
print(f"  residuo relativo max (num vs analitico, entrambi i rami): "
      f"{res_max:.2e}")
A_curve = np.linspace(A_g[0], A_g[-1], 40)
J_curve = np.array([curva_inversione(Av) for Av in A_curve])
n_inv = int(np.sum(d_num > 0))
print(f"  celle con INVERSIONE (tau piu' profondo): {n_inv}"
      f"/{np.sum(np.isfinite(d_num))}")
jj = np.isfinite(J_curve)
if jj.any():
    print(f"  curva analitica rho=1: A_inv da {A_curve[jj][0]:.3f} "
          f"(J={J_curve[jj][0]:.2f}) a {A_curve[jj][-1]:.3f} "
          f"(J={J_curve[jj][-1]:.2f})")

# --------------------------------------------------------------- figura
fig, ax = plt.subplots(figsize=(COL, COL * 0.92))
vmax = np.nanmax(np.abs(d_num))
im = ax.pcolormesh(J_g, A_g, d_num, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                   shading='auto')
plt.colorbar(im, ax=ax, label='$\\Delta r = r_{min}^{t} - r_{min}^{\\tau}$')
try:
    CS = ax.contour(J_g, A_g, d_num, levels=[0.0], colors='k',
                    linewidths=1.4)
except Exception:
    pass
ax.plot(J_curve, A_curve, 'w--', lw=1.8,
        label='analytic: $\\rho=A^2[f+2Ms/(rJ)]/\\hat E = 1$')
ax.axhline(1.0, color='gray', lw=0.7, ls=':')
ax.text(J_g[0] + 0.03, 1.005, 'Kerr ($A=1$): $t$ sinks deeper',
        fontsize=6.5, color='gray')
ax.set_xlabel('$J$')
ax.set_ylabel('$A$')
ax.set_title('plunge inversion $t$ vs $\\tau$ (Thakurta-Kerr)\n'
             f'($s={s_n}$, $\\hat E={E_n}$, $r_0={r0}$; blue: $t$ '
             'deeper, red: $\\tau$ deeper)\n'
             f'max num-analytic residual {res_max:.0e}')
ax.legend(loc='lower right')
savefig(fig, OUT, 'fig_thakurta_kerr_inversione_t_tau')
print('\nFATTO.')
