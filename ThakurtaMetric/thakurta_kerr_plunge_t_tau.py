# -*- coding: utf-8 -*-
"""
Plunge inversion t vs tau nel piano (a, J) — Kerr e Thakurta-Kerr.

PROTOCOLLO (evidence_t_plunges_deeper.txt): stesse condizioni iniziali
(r0, p_r0, J) nei flussi delle Hamiltoniane ottiche QUADRATICHE dei due
rami equatoriali; integrazione fino a p_r = 0; confronto del raggio
minimo raggiunto rispetto alla singolarita'. In Kerr: il ramo t affonda
di piu' (soppressione centrifuga 1/E^2), TRANNE a J piccolo dove il
dragging inverte l'ordinamento: l'inversione vive nel piano
(a = spin del BH, J = momento angolare z della particella).

Hamiltoniane quadratiche (equatoriali; a A=1 identiche al progetto Kerr):

  H_tau = (1/2k^2)[(f Dl/r^2) p_r^2 + (f^2/Dl) J^2]
  H_t   = (1/2n^2)[(f Dl/r^2) p_r^2 + (f^2/Dl) X^2],  X = J + 2Ma/(rf)
  n^2 = Ehat^2/(Ehat^2 - A^2 f),   k = A^2 f n/Ehat,   Dl = r^2-2Mr+a^2

RAPPORTO DELLE HAMILTONIANE (potenziali centrifughi a p_r=0):

  V_t/V_tau = rho^2 ,   rho = A^2 [ f + 2Ma/(rJ) ] / Ehat

  rho < 1: t piu' profondo;  rho > 1: tau piu' profondo.
  La curva di inversione J_inv(a) cresce con lo spin (dragging);
  il fattore conforme A^2 la SPOSTA verso J piu' alti: a A > sqrt(Ehat)
  l'inversione diventa generica (E_eff = Ehat/A < 1: vince la
  dilatazione temporale, crossover tipo TOV).

Con A costante H e' conservata: r_min = radice algebrica di
V(r) = H(r0,p_r0): luogo esatto di inversione algebrico; le ODE lo
validano. Figura: (a) colormap (J, a) a A=1 con luogo esatto + formula
rho=1; (b) drift del luogo esatto con A = 1, 1.1, 1.2.
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

M_n, E_n = 1.0, 1.2
r0, pr0 = 10.0, -0.5

r, pr, A, J, a_ = sp.symbols('r p_r A J a', positive=True)
f = 1 - 2 * M_n / r
Dl = r**2 - 2 * M_n * r + a_**2
n2 = E_n**2 / (E_n**2 - A**2 * f)
k2 = A**4 * f**2 * n2 / E_n**2
X = J + 2 * M_n * a_ / (r * f)

H_tau = ((f * Dl / r**2) * pr**2 + (f**2 / Dl) * J**2) / (2 * k2)
H_t = ((f * Dl / r**2) * pr**2 + (f**2 / Dl) * X**2) / (2 * n2)

lamb = lambda ex: sp.lambdify((r, pr, A, J, a_), ex, 'numpy')
Hs = {'tau': (lamb(H_tau), lamb(sp.diff(H_tau, pr)), lamb(sp.diff(H_tau, r))),
      't': (lamb(H_t), lamb(sp.diff(H_t, pr)), lamb(sp.diff(H_t, r)))}

def r_plus(av):
    return M_n + np.sqrt(max(M_n**2 - av**2, 1e-12))

def r_min_alg(Jv, av, Av, ramo):
    Hx = Hs[ramo][0]
    H0 = Hx(r0, pr0, Av, Jv, av)
    def g(rv):
        return Hx(rv, 0.0, Av, Jv, av) - H0
    rg = np.linspace(r0, r_plus(av) * 1.001, 4000)
    with np.errstate(invalid='ignore', divide='ignore'):
        gv = np.array([g(x) for x in rg])
    for i in range(len(rg) - 1):
        if np.isfinite(gv[i]) and np.isfinite(gv[i + 1]) \
                and gv[i] < 0 <= gv[i + 1]:
            return brentq(g, rg[i], rg[i + 1])
    return 0.0

def r_min_ode(Jv, av, Av, ramo):
    Hx, dHp, dHr = Hs[ramo]
    def rhs(t_, y):
        return [dHp(y[0], y[1], Av, Jv, av), -dHr(y[0], y[1], Av, Jv, av)]
    ev_t = lambda t_, y: y[1]
    ev_t.terminal, ev_t.direction = True, 1
    ev_h = lambda t_, y: y[0] - r_plus(av) * 1.002
    ev_h.terminal, ev_h.direction = True, -1
    s = solve_ivp(rhs, [0, 20000], [r0, pr0], rtol=1e-10, atol=1e-12,
                  method='DOP853', events=[ev_t, ev_h])
    if len(s.t_events[0]):
        return s.y_events[0][0][0]
    return 0.0

print("=" * 72)
print(f"[V1] Kerr A=1, a=0.4: confronto con l'evidence (r0={r0}, pr0={pr0})")
print("=" * 72)
for Jv in (0.5, 1.5, 2.5, 3.5, 4.5):
    rt = r_min_ode(Jv, 0.4, 1.0, 'tau')
    rtt = r_min_ode(Jv, 0.4, 1.0, 't')
    print(f"  J={Jv:4.1f}:  r_min(t)={rtt:7.4f}  r_min(tau)={rt:7.4f}"
          f"  Delta={rtt - rt:+8.4f}")
print("  evidence (a=0.4, r0=10): +0.06, -0.45(J=1) ... -0.42(J=4.5)")

print()
print("=" * 72)
print("[V2] colormap (J, a) a A=1 + luogo esatto + formula rho=1")
print("=" * 72)
J_g = np.linspace(0.05, 3.0, 28)
a_g = np.linspace(0.05, 0.98, 28)

def mappa_ode(Av):
    """colormap INTERAMENTE numerica (flussi ODE, nessuna algebra)."""
    d = np.full((len(a_g), len(J_g)), np.nan)
    for i, av in enumerate(a_g):
        for j, Jv in enumerate(J_g):
            rt_ = r_min_ode(Jv, av, Av, 't')
            ru_ = r_min_ode(Jv, av, Av, 'tau')
            if rt_ > 0 and ru_ > 0:
                d[i, j] = rt_ - ru_
    return d

print("  ... colormap ODE 28x28x2 in corso ...")
d1 = mappa_ode(1.0)

# validazione ODE a campione
res_max = 0.0
for av in (0.2, 0.5, 0.9):
    for Jv in (0.5, 1.5, 2.5):
        for ramo in ('t', 'tau'):
            ra = r_min_alg(Jv, av, 1.0, ramo)
            ro = r_min_ode(Jv, av, 1.0, ramo)
            if ra > 0 and ro > 0:
                res_max = max(res_max, abs(ra - ro) / ra)
print(f"  residuo ODE vs algebrico (campione): {res_max:.2e}")

def J_inv_esatta(av, Av):
    """luogo esatto: Delta r(J) = 0 (bisezione in J)."""
    def dd(Jv):
        rt_ = r_min_alg(Jv, av, Av, 't')
        ru_ = r_min_alg(Jv, av, Av, 'tau')
        if rt_ <= 0 or ru_ <= 0:
            return np.nan
        return rt_ - ru_
    Jg = np.linspace(0.02, 4.5, 120)
    vals = np.array([dd(Jv) for Jv in Jg])
    for i in range(len(Jg) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) \
                and vals[i] > 0 >= vals[i + 1]:
            return brentq(dd, Jg[i], Jg[i + 1])
    return np.nan

def rho_r(rv, Jv, av, Av):
    fv = 1 - 2 * M_n / rv
    return Av**2 * (fv + 2 * M_n * av / (rv * Jv)) / E_n

def J_inv_rho(av, Av):
    """FORMULA CORRETTA:  rho(r*)^2 = H_t0/H_tau0  alla svolta tau.

    I due rami partono con valori hamiltoniani DIVERSI (stesse
    (r0, p_r0, J) ma coefficienti diversi): il criterio del rapporto dei
    potenziali va normalizzato al rapporto di lancio, non a 1.
    In campo lontano H_t0/H_tau0 ~ rho(r0)^2: inversione dove
    rho(r*) = rho(r0)."""
    Ht0f = Hs['t'][0]
    Hu0f = Hs['tau'][0]
    def eq(Jv):
        ru_ = r_min_alg(Jv, av, Av, 'tau')
        if ru_ <= 0:
            return np.nan
        rapporto = Ht0f(r0, pr0, Av, Jv, av) / Hu0f(r0, pr0, Av, Jv, av)
        return rho_r(ru_, Jv, av, Av)**2 - rapporto
    Jg = np.linspace(0.02, 4.5, 120)
    vals = np.array([eq(Jv) for Jv in Jg])
    for i in range(len(Jg) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) \
                and vals[i] * vals[i + 1] <= 0:
            return brentq(eq, Jg[i], Jg[i + 1])
    return np.nan

a_c = np.linspace(0.05, 0.98, 30)
Ji_1 = np.array([J_inv_esatta(av, 1.0) for av in a_c])
Jr_1 = np.array([J_inv_rho(av, 1.0) for av in a_c])
print(f"  A=1: J_inv esatta da {np.nanmin(Ji_1):.2f} (a piccolo) a "
      f"{np.nanmax(Ji_1):.2f} (a={a_c[np.nanargmax(Ji_1)]:.2f})")
dev = np.nanmax(np.abs(Ji_1 - Jr_1) / Ji_1)
print(f"  formula corretta rho(r*)^2 = H_t0/H_tau0 vs luogo esatto: "
      f"deviazione max = {dev:.2e}")

# luogo di inversione NUMERICO (bisezione sui flussi ODE)
def J_inv_ode(av, Av):
    def dd(Jv):
        rt_ = r_min_ode(Jv, av, Av, 't')
        ru_ = r_min_ode(Jv, av, Av, 'tau')
        if rt_ <= 0 or ru_ <= 0:
            return np.nan
        return rt_ - ru_
    Jg = np.linspace(0.02, 4.5, 100)
    vals = np.array([dd(Jv) for Jv in Jg])
    for i in range(len(Jg) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) \
                and vals[i] > 0 >= vals[i + 1]:
            return brentq(dd, Jg[i], Jg[i + 1], xtol=1e-8)
    return np.nan

a_pts = np.linspace(0.05, 0.96, 14)
J_ode_pts = np.array([J_inv_ode(av, 1.0) for av in a_pts])
J_ana_pts = np.array([J_inv_rho(av, 1.0) for av in a_pts])
devs = np.abs(J_ode_pts - J_ana_pts)
print("  luogo NUMERICO (ODE) vs FORMULA analitica, punto per punto:")
for av, jo, ja in zip(a_pts, J_ode_pts, J_ana_pts):
    print(f"    a={av:.2f}:  J_inv ODE = {jo:.6f}   analitica = {ja:.6f}"
          f"   |diff| = {abs(jo - ja):.2e}")
print(f"  deviazione massima ODE vs analitica: {np.nanmax(devs):.2e}")

print()
print("=" * 72)
print("[V3] drift conforme del luogo di inversione nel piano (a, J)")
print("=" * 72)
curves = {1.0: Ji_1}
for Av in (1.1, 1.2):
    curves[Av] = np.array([J_inv_esatta(av, Av) for av in a_c])
    print(f"  A={Av}: J_inv(a=0.9) = {curves[Av][np.argmin(np.abs(a_c-0.9))]:.2f}"
          f"   (A=1: {Ji_1[np.argmin(np.abs(a_c-0.9))]:.2f})")

# --------------------------------------------------------------- figura
fig, (axa, axb) = plt.subplots(2, 1, figsize=(COL, 6.6))
vmax = np.nanmax(np.abs(d1))
im = axa.pcolormesh(J_g, a_g, d1, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                    shading='auto')
plt.colorbar(im, ax=axa, label='$\\Delta r = r_{min}^{t} - r_{min}^{\\tau}$')
CS = axa.contour(J_g, a_g, d1, levels=[0.0], colors='k', linewidths=1.8)
axa.plot(Jr_1, a_c, 'w--', lw=1.4,
         label='analytic $\\rho(r_*)^2 = H_{t0}/H_{\\tau 0}$')
axa.plot(J_ode_pts, a_pts, 'o', mfc='none', mec='lime', ms=6, mew=1.3,
         label='$\\Delta r=0$ numerical (ODE bisection)')
axa.set_xlabel('$J$ (particle angular momentum)')
axa.set_ylabel('$a$ (black-hole spin)')
axa.set_title(f'Kerr ($A=1$, $E={E_n}$, $r_0={r0}$, $p_r={pr0}$):\n'
              'colormap + black contour = PURE ODE; '
              f'ODE vs analytic max {np.nanmax(devs):.1e}')
axa.legend(loc='lower right')

axb.plot(J_ana_pts, a_pts, 'k-', lw=1.5, label='analytic')
axb.plot(J_ode_pts, a_pts, 'o', mfc='none', mec='C3', ms=7, mew=1.4,
         label='numerical (ODE)')
axb.set_xlabel('$J_{inv}(a)$')
axb.set_ylabel('$a$')
axb.set_title('numerical/analytic overlap of the inversion\n'
              f'locus: max $|\\Delta J| = {np.nanmax(devs):.1e}$')
axb.grid(alpha=0.3)
axb.legend()
savefig(fig, OUT, 'fig_thakurta_kerr_plunge_t_tau')
print('\nFATTO.')
