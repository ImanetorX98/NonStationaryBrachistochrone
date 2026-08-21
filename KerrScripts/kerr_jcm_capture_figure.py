# -*- coding: utf-8 -*-
"""
Figura informativa: cattura dall'ergosfera per J=0.95 J_c^- (ramo t retrogrado).

Al saddle-node J_c^- il sestico R6 ha una radice DOPPIA reale a r_*=3.514>r_e:
barriera a phi infinito, l'orbita spirala su r_* e NON penetra.  Spostando J
appena dentro la finestra (J=0.95 J_c^-=-7.65) la radice doppia si stacca in
COPPIA COMPLESSA (3.393+-0.601i): la barriera sparisce, l'orbita fa PLUNGE
attraverso l'ergosfera fino all'orizzonte -> CATTURATA.

Pannello (a): traiettoria J=0.95 J_c^- (plunge, penetra r_e) -- flusso di
  Hamilton vs forma chiusa genere-2; separatrice J_c^- (spirale su r_*) grigia.
Pannello (b): meccanismo saddle-node -- le due radici rilevanti di R6 nel piano
  complesso al variare di J da J_c^- (doppia reale) a 0.95 J_c^- (coppia
  complessa): collisione -> distacco = cattura.
"""

import os
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq, fsolve

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
re = 2*M


def Q2(r, J):
    return (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r
            + E**2*a**2*r**2 + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2
            - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)


def R6(r, J):
    return r*Q2(r, J)*((E**2 - 1)*r + 2*M)


def K_gen(r, J):
    return r*((E**2 - 1)*r + 2*M)*(J*(r - 2*M) + 2*M*a)/(r**2 - 2*M*r + a**2)


def roots6(J):
    rs = sp.Symbol('r')
    c = [float(x) for x in sp.Poly(
        sp.expand(rs*Q2(rs, J)*((E**2-1)*rs+2*M)), rs).all_coeffs()]
    return np.roots(c)


# saddle-node J_c^-
rst, Jst = fsolve(lambda x: [R6(x[0], x[1]),
                             (R6(x[0]+1e-7, x[1]) - R6(x[0]-1e-7, x[1]))/2e-7],
                  [3.51, -8.05])
J = 0.95*Jst
print(f"J_c^-={Jst:.4f} (r_*={rst:.4f}); J=0.95 J_c^-={J:.4f}")

# --- flusso di Hamilton (ramo t) per J=0.95 J_c^- ---
r, pr, pphi = sp.symbols('r pr pphi')
f = 1 - 2*M/r
Dl = r**2 - 2*M*r + a**2
b = 2*M*a/r
v2 = 1 - f/E**2
P = r**2 + a**2 + 2*M*a**2/r
Pb = P + b**2/E**2
H = pphi*b*v2/Pb + sp.sqrt(Dl*v2/Pb)*sp.sqrt((Dl/r**2)*pr**2 + pphi**2/Pb) - 1
Hn = sp.lambdify((r, pr), H.subs(pphi, J), 'numpy')
dHp = sp.lambdify((r, pr), sp.diff(H, pr).subs(pphi, J), 'numpy')
dHr = sp.lambdify((r, pr), sp.diff(H, r).subs(pphi, J), 'numpy')
dHj = sp.lambdify((r, pr), sp.diff(H, pphi).subs(pphi, J), 'numpy')


def prof(rv):
    pg = np.linspace(-200, 200, 20001)
    Hv = Hn(rv, pg)
    rts = [brentq(lambda p: Hn(rv, p), pg[i], pg[i+1]) for i in range(len(pg)-1)
           if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
    ing = [p for p in rts if dHp(rv, p) < 0]
    return min(ing) if ing else np.nan


r0 = 6.0
p0 = prof(r0)
sA = solve_ivp(lambda t, y: [dHp(y[0], y[1]), -dHr(y[0], y[1]), dHj(y[0], y[1])],
               [0, 2000], [r0, p0, 0.0], rtol=1e-11, atol=1e-13, max_step=0.005,
               dense_output=True, events=lambda t, y: y[0] - (rp + 0.03))
te = sA.t_events[0][0]
tA = np.linspace(0, te, 700)
rA = sA.sol(tA)[0]
phiA = sA.sol(tA)[2]
integ = lambda rv: K_gen(rv, J)/np.sqrt(abs(R6(rv, J)))
phiR = np.array([quad(integ, rv, r0, limit=500)[0] for rv in rA])
devA = np.abs(phiA - phiR)
print(f"orbita: r {rA[0]:.2f}->{rA[-1]:.3f}; penetra r_e? {'SI' if rA.min()<re else 'NO'}; "
      f"cattura r_+? {'SI' if rA.min()<rp+0.06 else 'no'}")
print(f"max|phi_ODE - phi_genus2| = {devA.max():.2e}")

# --- separatrice J_c^-: spirale su r_* (per confronto) ---
p0s = prof(r0)  # stessa profondita' iniziale non serve; uso il flusso a J_c^-
Hns = sp.lambdify((r, pr), H.subs(pphi, Jst), 'numpy')
dHps = sp.lambdify((r, pr), sp.diff(H, pr).subs(pphi, Jst), 'numpy')
dHrs = sp.lambdify((r, pr), sp.diff(H, r).subs(pphi, Jst), 'numpy')
dHjs = sp.lambdify((r, pr), sp.diff(H, pphi).subs(pphi, Jst), 'numpy')


def profs(rv):
    pg = np.linspace(-200, 200, 20001)
    Hv = Hns(rv, pg)
    rts = [brentq(lambda p: Hns(rv, p), pg[i], pg[i+1])
           for i in range(len(pg)-1)
           if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
    ing = [p for p in rts if dHps(rv, p) < 0]
    return min(ing) if ing else np.nan


ss = solve_ivp(lambda t, y: [dHps(y[0], y[1]), -dHrs(y[0], y[1]),
                             dHjs(y[0], y[1])],
               [0, 4000], [r0, profs(r0), 0.0], rtol=1e-11, atol=1e-13,
               max_step=0.01, dense_output=True,
               events=lambda t, y: y[0] - (rst + 0.02))
ts = np.linspace(0, ss.t_events[0][0], 700)
rS = ss.sol(ts)[0]
phiS = ss.sol(ts)[2]

# --- moto delle due radici rilevanti da J_c^- a 0.95 J_c^- ---
Jscan = np.linspace(Jst, J, 60)
tracks = []
for Jv in Jscan:
    rr = roots6(Jv)
    near = sorted(rr, key=lambda z: abs(z - rst))[:2]
    tracks.append(sorted(near, key=lambda z: z.imag))
tracks = np.array(tracks)

# ----------------------------------------------------------- figura
fig, (ax, axb) = plt.subplots(2, 1, figsize=(COL, 6.8))
th = np.linspace(0, 2*np.pi, 400)
# (a) traiettoria
ax.fill(re*np.cos(th), re*np.sin(th), color='0.9', zorder=0)
ax.plot(re*np.cos(th), re*np.sin(th), 'b--', lw=1.0, label=r'ergosphere $r_e$')
ax.plot(rp*np.cos(th), rp*np.sin(th), 'C3--', lw=0.9, label=r'horizon $r_+$')
ax.plot(rst*np.cos(th), rst*np.sin(th), color='0.6', ls=':', lw=0.9,
        label=r'marginal circle $r_*$')
ax.plot(rS*np.cos(phiS), rS*np.sin(phiS), color='0.6', lw=1.1,
        label=rf'$J_c^-={Jst:.2f}$ (spirals onto $r_*$)')
ax.plot(rA*np.cos(phiA), rA*np.sin(phiA), 'k-', lw=2.3,
        label=rf'$0.95\,J_c^-={J:.2f}$ ODE (captured)')
ax.plot(rA*np.cos(phiR), rA*np.sin(phiR), 'C2:', lw=1.8,
        label='genus-2 closed form')
ax.set_aspect('equal')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title('retrograde capture: $0.95\\,J_c^-$ plunges through $r_e$\n'
             'to the horizon; $J_c^-$ only spirals onto $r_*$')
ax.legend(fontsize=5.3, loc='upper left', framealpha=0.9)

# (b) saddle-node: le due radici nel piano complesso
axb.plot(tracks[:, 0].real, tracks[:, 0].imag, 'C0.-', ms=2, lw=0.9)
axb.plot(tracks[:, 1].real, tracks[:, 1].imag, 'C0.-', ms=2, lw=0.9)
axb.plot(rst, 0, 'ks', ms=7, label=r'$J_c^-$: double root at $r_*$ (real)')
axb.plot([tracks[-1, 0].real, tracks[-1, 1].real],
         [tracks[-1, 0].imag, tracks[-1, 1].imag], 'C3o', ms=6,
         label=r'$0.95\,J_c^-$: complex pair')
axb.axhline(0, color='0.7', lw=0.6)
axb.set_xlabel(r'$\operatorname{Re}\,r$')
axb.set_ylabel(r'$\operatorname{Im}\,r$')
axb.set_title('saddle-node: double root splits into a complex pair\n'
              '(barrier removed $\\Rightarrow$ capture)')
axb.legend(fontsize=6, loc='upper right', framealpha=0.9)
savefig(fig, HERE, 'fig_jcm_capture')
print('FATTO.')
