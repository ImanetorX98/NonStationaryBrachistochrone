# -*- coding: utf-8 -*-
"""
Figura: penetrazione dell'ergosfera da parte delle orbite PROGRADE del ramo t
nella finestra di cattura [J_c^-, J_+^t], GENERE 2 (Kleinian).

Il valore soglia J_+^t=(4M^2+2a^2+a^2/E^2)/(2a) ha il turning point ESATTAMENTE
su r_e=2M (radice semplice di R6 a r=2): tocca l'ergosfera e torna.  Per J
appena dentro la finestra (qui J=3.0) il turning si sposta DENTRO l'ergoregione
(r_in=1.627<r_e): l'orbita PENETRA l'ergosfera.  Il radicando resta sestico a 6
radici semplici -> genere 2 (curva iperellittica y^2=R6, non Weierstrass).

Pannello (a): traiettorie -- J_+^t marginale (tocca r_e) e J=3.0 penetrante,
  flusso di Hamilton (ODE) vs forma chiusa genere-2 (integrali abeliani di 1a+3a
  specie), valida ATTRAVERSO e DENTRO l'ergosfera.
Pannello (b): deviazione |phi_ODE - phi_riduzione| per J=3.0, con la regione
  r<r_e evidenziata: la forma chiusa resta valida ~2e-13 dentro l'ergosfera.

dphi/dr = K(r)/sqrt(R6),  K = r[(E^2-1)r+2M](J(r-2M)+2Ma)/Delta   (GENERALE)
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
# paper_style.py lives in NonStationaryMetrics/; a fresh clone must find it
# without relying on an unarchived sibling directory (CQG-116884, major 10)
sys.path.insert(0, os.path.join(_ROOT, "NonStationaryMetrics"))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a, E = 1.0, 0.9, 1.2
rp = M + np.sqrt(M**2 - a**2)
re = 2*M                                       # ergosfera equatoriale
Jpt = (4*M**2 + 2*a**2 + a**2/E**2)/(2*a)      # J_+^t marginale (tocca r_e)


def Q2(r, J):
    return (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r
            + E**2*a**2*r**2 + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2
            - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)


def R6(r, J):
    return r*Q2(r, J)*((E**2 - 1)*r + 2*M)


def K_gen(r, J):
    return r*((E**2 - 1)*r + 2*M)*(J*(r - 2*M) + 2*M*a)/(r**2 - 2*M*r + a**2)


def flow(J, r_stop):
    """flusso di Hamilton del ramo t da r0=6 fino a r_stop."""
    r, pr, pphi = sp.symbols('r pr pphi')
    f = 1 - 2*M/r
    Dl = r**2 - 2*M*r + a**2
    b = 2*M*a/r
    v2 = 1 - f/E**2
    P = r**2 + a**2 + 2*M*a**2/r
    Pb = P + b**2/E**2
    H = (pphi*b*v2/Pb + sp.sqrt(Dl*v2/Pb)*sp.sqrt((Dl/r**2)*pr**2 + pphi**2/Pb)
         - 1)
    Hn = sp.lambdify((r, pr), H.subs(pphi, J), 'numpy')
    dHp = sp.lambdify((r, pr), sp.diff(H, pr).subs(pphi, J), 'numpy')
    dHr = sp.lambdify((r, pr), sp.diff(H, r).subs(pphi, J), 'numpy')
    dHj = sp.lambdify((r, pr), sp.diff(H, pphi).subs(pphi, J), 'numpy')

    def prof(rv):
        pg = np.linspace(-100, 100, 10001)
        Hv = Hn(rv, pg)
        rts = [brentq(lambda p: Hn(rv, p), pg[i], pg[i+1])
               for i in range(len(pg)-1)
               if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1])
               and Hv[i]*Hv[i+1] < 0]
        ing = [p for p in rts if dHp(rv, p) < 0]
        return min(ing) if ing else np.nan

    r0 = 6.0
    p0 = prof(r0)
    s = solve_ivp(lambda t, y: [dHp(y[0], y[1]), -dHr(y[0], y[1]),
                                dHj(y[0], y[1])],
                  [0, 600], [r0, p0, 0.0], rtol=1e-12, atol=1e-14,
                  max_step=0.005, dense_output=True,
                  events=lambda t, y: y[0] - r_stop)
    te = s.t_events[0][0]
    tt = np.linspace(0, te, 500)
    rr = s.sol(tt)[0]
    ph = s.sol(tt)[2]
    # forma chiusa: phi(r) = int_r^{r0} K/sqrt(R6)
    integ = lambda rv: K_gen(rv, J)/np.sqrt(abs(R6(rv, J)))
    phc = np.array([quad(integ, rv, r0, limit=400)[0] for rv in rr])
    return r0, rr, ph, phc


# --- J_+^t marginale (tocca r_e) e J=3.0 penetrante ---
_, rM, phM, _ = flow(Jpt, re + 0.02)
Jp = 3.0
r_in = 1.6272
_, rP, phP, phPc = flow(Jp, r_in + 0.02)
devP = np.abs(phP - phPc)
print(f"J_+^t={Jpt:.4f} tocca r_e; J={Jp} penetra a r_in~{rP.min():.3f}<r_e={re}")
print(f"max|phi_ODE - phi_genus2| (J={Jp}) = {devP.max():.2e}, "
      f"dentro ergosfera = {devP[rP<re].max():.2e}")

# ----------------------------------------------------------- figura
fig, (ax, axb) = plt.subplots(2, 1, figsize=(COL, 6.7))
th = np.linspace(0, 2*np.pi, 400)
# (a) traiettorie
ax.fill(re*np.cos(th), re*np.sin(th), color='0.9', zorder=0)          # ergoregione
ax.plot(re*np.cos(th), re*np.sin(th), 'b--', lw=1.0, label=r'ergosphere $r_e=2M$')
ax.plot(rp*np.cos(th), rp*np.sin(th), 'C3--', lw=0.9, label=r'horizon $r_+$')
ax.plot(rM*np.cos(phM), rM*np.sin(phM), color='0.55', lw=1.3,
        label=rf'$J_+^t={Jpt:.2f}$ (touches $r_e$)')
ax.plot(rP*np.cos(phP), rP*np.sin(phP), 'k-', lw=2.4,
        label=rf'$J={Jp}$ ODE (penetrates)')
ax.plot(rP*np.cos(phPc), rP*np.sin(phPc), 'C2:', lw=1.9,
        label=r'genus-2 closed form')
ax.set_aspect('equal')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title('prograde $t$-orbits penetrate the ergosphere\n'
             '(genus 2; closed form valid across $r_e$)')
ax.legend(fontsize=5.6, loc='lower left', framealpha=0.9)

# (b) deviazione con ergoregione evidenziata
axb.axvspan(rP.min(), re, color='0.9', label='inside ergosphere')
axb.semilogy(rP, devP + 1e-18, 'k-')
axb.axvline(re, color='b', ls='--', lw=1.0)
axb.axvline(rp, color='C3', ls='--', lw=0.9)
axb.text(re + 0.05, devP.max()*0.4, '$r_e$', color='b', fontsize=7)
axb.text(rp - 0.13, devP.max()*0.4, '$r_+$', color='C3', fontsize=7)
axb.set_xlabel('$r$')
axb.set_ylabel(r'$|\phi_{\rm ODE}-\phi_{\rm closed}|$')
axb.set_title(rf'genus-2 closed form vs flow ($J={Jp}$): '
              r'valid inside $r_e$')
axb.legend(fontsize=6, loc='upper right', framealpha=0.9)
savefig(fig, HERE, 'fig_jpt_genus2')
print('FATTO.')
