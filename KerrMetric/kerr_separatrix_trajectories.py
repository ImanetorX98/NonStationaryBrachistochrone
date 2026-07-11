# -*- coding: utf-8 -*-
"""
TEST FINALE della forma chiusa R12a: tre traiettorie della separatrice
J = J_c = a/E (Kerr equatoriale, ramo tau), stesse condizioni iniziali,
evolute fin dentro l'ergosfera:

  T1  ODE: flusso di Hamilton PMP (integratore, nessuna forma chiusa)
  T2  quadratura diretta di dphi/dr = (a/E) sqrt(Q4)/(Dl P2)
  T3  FORMULA DI WEIERSTRASS valutata davvero: sigma/zeta/P costruite
      con le theta di Jacobi sul reticolo VERO (rombico, tau = 1/2+iy)

Se si sovrappongono, la forma chiusa e' verificata end-to-end,
inclusa la valutazione delle funzioni speciali.
Figura: fig_separatrix_3traiettorie (in KerrMetric/).
"""

import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig

mp.mp.dps = 25
M, a, E = 1.0, 0.9, 1.2
c = a / E
Jc = a / E
rp = M + np.sqrt(M**2 - a**2)
rm = M - np.sqrt(M**2 - a**2)
r0 = 9.0
r_end = rp + 0.05

Q4 = lambda r: r * (2 * M + (E**2 - 1) * r) * (r**2 + c**2)
Dl = lambda r: r**2 - 2 * M * r + a**2
P2 = lambda r: r**2 + c**2
dphidr = lambda r: (a / E) * np.sqrt(Q4(r)) / (Dl(r) * P2(r))

# ------------------------------------------------- T2: quadratura
rg = np.linspace(r_end, r0, 400)
phi2 = np.array([quad(dphidr, rv, r0, limit=300)[0] for rv in rg])

# ------------------------------------------------- T1: ODE Hamilton
# H_tau PMP equatoriale Kerr (R7 con A=1)
import sympy as sp
r_s, pr_s = sp.symbols('r p_r', real=True)
f_s = 1 - 2 * M / r_s
Dl_s = r_s**2 - 2 * M * r_s + a**2
P_s = r_s**2 + a**2 + 2 * M * a**2 / r_s
vb2 = 1 - f_s / E**2
Pb = P_s + (2 * M * a / r_s)**2 / E**2
php0 = (2 * M * a / r_s) * vb2 / Pb
R2s = Dl_s * vb2 / Pb
ptp = Jc - 2 * M * a / (r_s * E)
H = ptp * php0 + sp.sqrt(R2s) * sp.sqrt((Dl_s / r_s**2) * pr_s**2
                                        + ptp**2 / Pb) - f_s / E
Hf = sp.lambdify((r_s, pr_s), H, 'numpy')
dHp = sp.lambdify((r_s, pr_s), sp.diff(H, pr_s), 'numpy')
dHr = sp.lambdify((r_s, pr_s), sp.diff(H, r_s), 'numpy')
dphi_val = sp.lambdify((r_s, pr_s),
                       ptp * (php0 * 0 + 1) * 0 + sp.diff(H, sp.Symbol('x'))
                       if False else
                       (php0 * ptp / ptp), 'numpy')  # placeholder
# dphi/deta = dH/dp_phi: ricostruita simbolicamente
pphi_s = sp.Symbol('p_phi', real=True)
ptp_g = pphi_s - 2 * M * a / (r_s * E)
H_g = ptp_g * php0 + sp.sqrt(R2s) * sp.sqrt((Dl_s / r_s**2) * pr_s**2
                                            + ptp_g**2 / Pb) - f_s / E
dHdJ = sp.lambdify((r_s, pr_s),
                   sp.diff(H_g, pphi_s).subs(pphi_s, Jc), 'numpy')

def lancio_ode(r_start, span, ev_list):
    with np.errstate(invalid='ignore'):
        pg = np.linspace(-200, 200, 400001)
        Hg = Hf(r_start, pg)
    roots = [brentq(lambda p: Hf(r_start, p), pg[i], pg[i + 1])
             for i in range(len(pg) - 1)
             if np.isfinite(Hg[i]) and np.isfinite(Hg[i + 1])
             and Hg[i] * Hg[i + 1] <= 0]
    p0 = [p_ for p_ in roots if dHp(r_start, p_) < 0][0]   # ingoing

    def rhs(t_, y):
        return [dHp(y[0], y[1]), -dHr(y[0], y[1]), dHdJ(y[0], y[1])]
    return solve_ivp(rhs, span, [r_start, p0, 0.0], rtol=1e-12,
                     atol=1e-14, method='DOP853', events=ev_list,
                     dense_output=True, max_step=2.0)

# ramo ESTERNO: da r0=9 in caduta, fino alla svolta marginale
ev_turn = lambda t_, y: y[1]
ev_turn.terminal, ev_turn.direction = True, 1
sA = lancio_ode(r0, [0, 4000], [ev_turn])
tA = np.linspace(0, sA.t[-1], 600)
rA = sA.sol(tA)[0]
phiA = sA.sol(tA)[2]
print(f"ODE ramo esterno: svolta marginale a r = {rA.min():.6f} (r_e = 2)")

# NOTA: dentro l'ergoregione H(p) > 0 per ogni p: la descrizione PMP-eta
# del ramo tau TERMINA a r_e (il costo F = f + b phi' cambia segno sul
# ramo spurio dell'ellisse): versione hamiltoniana della singolarita'
# intrinseca della riduzione (doranTau). L'ODE puo' solo sfiorare r_e;
# la forma chiusa segue l'orbita del worldline vincolato attraverso
# (validata dentro dal confronto Doran di doranTau: max|dphi| = 4.4e-6).

# ------------------------------------------------- T3: Weierstrass vero
A_ = M * a**2 / (2 * E**2)
B_ = a**2 * (E**2 - 1) / (12 * E**2)
g2 = a**4 * (E**2 - 1)**2 / (12 * E**4) - a**2 * M**2 / E**2
g3 = a**4 * (36 * E**2 * M**2 * (1 - E**2)
             - a**2 * (E**2 - 1)**3) / (216 * E**6)

Q4m = lambda r: r * (2 * M + (E**2 - 1) * r) * (r**2 + c**2)
z_of_r = lambda rv: mp.quad(lambda t: 1 / mp.sqrt(Q4m(t)), [0, rv])

# periodo reale: ramo r>0 fino a oo + ramo r<r_neg
r_neg = -2 * M / (E**2 - 1)
w1 = mp.quad(lambda t: 1 / mp.sqrt(Q4m(t)), [0, mp.inf]) \
    + mp.quad(lambda t: 1 / mp.sqrt(abs(Q4m(t))), [-mp.inf, r_neg])
print(f"omega_1 (semiperiodo reale) = {mp.nstr(w1, 10)}")

# reticolo rombico: tau = 1/2 + i y — trova y col fit dell'identita'
def make_wz(y_):
    tau = mp.mpc(0.5, y_)
    q = mp.exp(1j * mp.pi * tau)
    th1 = lambda v: mp.jtheta(1, v, q)
    th1p = lambda v: mp.jtheta(1, v, q, 1)
    th1ppp = lambda v: mp.jtheta(1, v, q, 3)
    eta1 = -(mp.pi**2 / (12 * w1)) * th1ppp(0) / th1p(0)

    def zeta_(z):
        v = mp.pi * z / (2 * w1)
        return eta1 * z / w1 + (mp.pi / (2 * w1)) * th1p(v) / th1(v)

    def sigma_(z):
        v = mp.pi * z / (2 * w1)
        return (2 * w1 / mp.pi) * mp.exp(eta1 * z**2 / (2 * w1)) \
            * th1(v) / th1p(0)

    wp_ = lambda z: -mp.diff(zeta_, z)
    return zeta_, sigma_, wp_

z_test = mp.mpf(0.6) * w1

def resid(y_):
    _, _, wp_ = make_wz(y_)
    p = wp_(z_test)
    return float(abs(mp.re(p)**2 * 0 + (mp.diff(wp_, z_test))**2
                     - (4 * p**3 - g2 * p - g3)))

ys = np.linspace(0.2, 2.0, 37)
res = [resid(y_) for y_ in ys]
y0 = ys[int(np.argmin(res))]
y_opt = float(mp.findroot(lambda y_: mp.mpf(resid(float(y_))), y0,
                          solver='mnewton', tol=1e-18)) \
    if False else y0
# raffina con sezione aurea sul residuo
lo, hi = y0 - 0.06, y0 + 0.06
for _ in range(60):
    m1 = lo + 0.382 * (hi - lo)
    m2 = lo + 0.618 * (hi - lo)
    if resid(m1) < resid(m2):
        hi = m2
    else:
        lo = m1
y_opt = 0.5 * (lo + hi)
zeta_w, sigma_w, wp_w = make_wz(y_opt)
print(f"tau = 1/2 + {y_opt:.10f} i;  residuo identita' Weierstrass: "
      f"{resid(y_opt):.2e}")
# check indipendente: P(z(r*)) = A/r* + B
r_star = 5.0
z_star = z_of_r(r_star)
err_p = abs(wp_w(z_star) - (A_ / r_star + B_))
print(f"check uniformizzazione: |P(z(r))-(A/r+B)| a r=5: {mp.nstr(err_p,3)}")

# costanti della formula
alpha = {}
lam = {}
vv = {}
for rk, sgn in ((rp, +1), (rm, -1)):
    other = rm if rk == rp else rp
    alpha[rk] = rk * ((E**2 - 1) * rk + 2 * M) / (rk - other)
    lam[rk] = alpha[rk] / np.sqrt(Q4(rk))
    vv[rk] = z_of_r(rk)
Lam0 = (E**2 - 1) - alpha[rp] / rp - alpha[rm] / rm \
    + 2 * float(mp.re(zeta_w(vv[rp]))) * lam[rp] \
    + 2 * float(mp.re(zeta_w(vv[rm]))) * lam[rm]

def Phi(z):
    s_ = Lam0 * z
    for rk in (rp, rm):
        s_ += lam[rk] * mp.log(sigma_w(z - vv[rk]) / sigma_w(z + vv[rk]))
    return (a / E) * s_

z0 = z_of_r(r0)
phi3 = []
for rv in rg:
    val = Phi(z0) - Phi(z_of_r(rv))
    phi3.append(float(mp.re(val)))
phi3 = np.array(phi3)

# ------------------------------------------------- confronto e figura
devA = np.abs(phiA - np.interp(rA, rg, phi2))
print(f"max |phi_Weierstrass - phi_quad| = {np.max(np.abs(phi3 - phi2)):.2e}")
print(f"ODE vs forma chiusa (r>2.2): {np.max(devA[rA > 2.2]):.2e}")
print(f"ODE vs forma chiusa (r>2.05): {np.max(devA[rA > 2.05]):.2e}")

set_style()
fig, (ax, axb) = plt.subplots(2, 1, figsize=(COL, 6.6))
ax.plot(rA * np.cos(phiA), rA * np.sin(phiA), 'k-', lw=2.6,
        label='T1: Hamilton ODE (grazes $r_e$ at 2.000000)')
ax.plot(rg * np.cos(phi2), rg * np.sin(phi2), 'C1--', lw=1.6,
        label='T2: closed-form quadrature')
ax.plot(rg * np.cos(phi3), rg * np.sin(phi3), 'C2:', lw=1.8,
        label='T3: Weierstrass ($\\sigma,\\zeta,\\wp$)')
th = np.linspace(0, 2 * np.pi, 200)
ax.plot(2 * M * np.cos(th), 2 * M * np.sin(th), 'b--', lw=1.1,
        label='ergosphere $r_e=2M$')
ax.plot(rp * np.cos(th), rp * np.sin(th), 'k-', lw=1.0)
ax.set_aspect('equal')
ax.set_xlim(-6, 9.6)
ax.set_ylim(-4, 9.6)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title(f'separatrix $J=J_c=a/E$ ($a={a}$, $E={E}$): three methods\n'
             'closed form crosses $r_e$; PMP-$\\eta$ terminates there')
ax.legend(loc='lower left')

axb.semilogy(rg, np.abs(phi3 - phi2) + 1e-18, 'C2-',
             label='|Weierstrass $-$ quadrature|')
axb.semilogy(rA, devA + 1e-18, 'k-', label='|ODE $-$ closed form|')
axb.axvline(2 * M, color='b', ls='--', lw=1.0)
axb.text(2.05, 1e-4, 'ergosphere', fontsize=6.5, color='b')
axb.set_xlabel('$r$')
axb.set_ylabel('deviation in $\\phi(r)$')
axb.set_xlim(1.4, 9)
axb.set_title('quantitative overlap\n(ODE growth near $r_e$ = '
              'separatrix instability)')
axb.legend()
HERE = os.path.dirname(os.path.abspath(__file__))
savefig(fig, HERE, 'fig_separatrix_3traiettorie')
print('salvata fig_separatrix_3traiettorie')
