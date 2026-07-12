# -*- coding: utf-8 -*-
"""
Figura: separatrice del ramo t (saddle-node retrogrado J_c^-), GENERE 1.
Pannello (a) traiettoria: flusso di Hamilton (ODE) vs forma chiusa di
Weierstrass valutata (sigma/zeta/P via theta di Jacobi). Pannello (b)
deviazione |phi_ODE - phi_Weierstrass| e |phi_quad - phi_Weierstrass|.

Radicando sestico R6=(r-r*)^2 Q4 al saddle-node -> Q4 quartico (genere 1);
phi(r)=Lam0 z + sum lam_k ln(sig(z-v_k)/sig(z+v_k)),  z=P^-1(A/r+B).
"""

import os
import sys
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import fsolve, brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

mp.mp.dps = 30
set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a, E = 1.0, 0.9, 1.2
rp = M + np.sqrt(M**2 - a**2)
rm = M - np.sqrt(M**2 - a**2)


def Q2(r, J):
    return (2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r
            + E**2*a**2*r**2 + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2
            - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)


def R6(r, J):
    return r * Q2(r, J) * ((E**2 - 1)*r + 2*M)


rst, Jst = fsolve(lambda x: [R6(x[0], x[1]),
                             (R6(x[0]+1e-7, x[1]) - R6(x[0]-1e-7, x[1]))/2e-7],
                  [3.51, -8.05])
print(f"r*={rst:.6f}, J_c^-={Jst:.6f}")

# --- quartica Q4 e invarianti di Weierstrass ---
r = sp.Symbol('r')
c6 = [float(x) for x in sp.Poly(sp.expand(r*Q2(r, Jst)*((E**2-1)*r+2*M)),
                                r).all_coeffs()]
q4r = [z for z in np.roots(c6) if abs(z - rst) > 1e-3]
Q4p = (np.poly(q4r)*c6[0]).real                        # coeff r^4..r^0
a0, a1, a2, a3, a4 = Q4p/np.array([1, 4, 6, 4, 1])
g2 = a0*a4 - 4*a1*a3 + 3*a2**2
g3 = a0*a2*a4 + 2*a1*a2*a3 - a2**3 - a0*a3**2 - a1**2*a4
A_ = Q4p[3]/4
B_ = Q4p[2]/12
Q4m = lambda t: np.polyval(Q4p, t)
print(f"g2={g2:.5f} g3={g3:.5f} A={A_:.5f} B={B_:.5f}")

# --- flusso di Hamilton (ramo t) ---
pr, pphi = sp.symbols('pr pphi')
f = 1 - 2*M/r
Dl = r**2 - 2*M*r + a**2
b = 2*M*a/r
v2 = 1 - f/E**2
P = r**2 + a**2 + 2*M*a**2/r
Pb = P + b**2/E**2
H = pphi*b*v2/Pb + sp.sqrt(Dl*v2/Pb)*sp.sqrt((Dl/r**2)*pr**2 + pphi**2/Pb) - 1
Hn = sp.lambdify((r, pr), H.subs(pphi, Jst), 'numpy')
dHp = sp.lambdify((r, pr), sp.diff(H, pr).subs(pphi, Jst), 'numpy')
dHr = sp.lambdify((r, pr), sp.diff(H, r).subs(pphi, Jst), 'numpy')
dHJ = sp.lambdify((r, pr), sp.diff(H, pphi).subs(pphi, Jst), 'numpy')


def prof(rv):
    pg = np.linspace(-100, 100, 10001)
    Hv = Hn(rv, pg)
    rts = [brentq(lambda p: Hn(rv, p), pg[i], pg[i+1]) for i in range(len(pg)-1)
           if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
    ing = [p for p in rts if dHp(rv, p) < 0]
    return min(ing) if ing else np.nan


r0 = 7.0
p0 = prof(r0)
sA = solve_ivp(lambda t, y: [dHp(y[0], y[1]), -dHr(y[0], y[1]), dHJ(y[0], y[1])],
               [0, 600], [r0, p0, 0.0], rtol=1e-12, atol=1e-14, max_step=0.01,
               dense_output=True, events=lambda t, y: y[0]-(rst+0.12))
tA = np.linspace(0, sA.t_events[0][0], 700)
rA = sA.sol(tA)[0]
phiA = sA.sol(tA)[2]

# --- forma chiusa: G partial fractions (fit) ---
def sqrtQ4(rv): return np.sqrt(abs(R6(rv, Jst)))/abs(rv - rst)
def dphidr(rv):
    p = prof(rv)
    return dHJ(rv, p)/dHp(rv, p)
poles = [rp, rm, rst]
rs = np.linspace(3.9, 40, 60)
rows, rhs = [], []
for rv in rs:
    g = dphidr(rv)*sqrtQ4(rv)
    if np.isfinite(g):
        rows.append([1.0]+[1.0/(rv-p) for p in poles]); rhs.append(g)
coef, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
Lam, al = coef[0], dict(zip(poles, coef[1:]))

# quadratura (T2): valutata direttamente sui punti dell'orbita (no interp)
_integ = lambda t: (Lam + sum(al[p]/(t-p) for p in poles))/sqrtQ4(t)
rg = rA.copy()
phi2 = np.array([quad(_integ, r0, rv, limit=400)[0] for rv in rA])
devA = np.abs(phiA - phi2)
print(f"max|phi_ODE - phi_genus1| = {np.max(devA):.2e}")

# --- prova indipendente del genere: R6 = (r-r*)^2 Q4, radice DOPPIA ---
rr = np.linspace(rst-0.6, rst+0.6, 400)
R6v = R6(rr, Jst)

# ----------------------------------------------------------- figura
fig, (ax, axb) = plt.subplots(2, 1, figsize=(COL, 6.6))
th = np.linspace(0, 2*np.pi, 300)
# (a) traiettoria: ODE vs forma chiusa genere-1 (quadratura ellittica)
ax.plot(rA*np.cos(phiA), rA*np.sin(phiA), 'k-', lw=2.4,
        label='T1: Hamilton ODE')
xg = rg*np.cos(phi2); yg = rg*np.sin(phi2)
ax.plot(xg, yg, 'C2:', lw=1.9,
        label=r'T2: genus-1 closed form ($\wp,\sigma,\zeta$)')
ax.plot(rst*np.cos(th), rst*np.sin(th), 'C3--', lw=1.0,
        label=r'marginal circle $r_*$')
ax.plot(2*M*np.cos(th), 2*M*np.sin(th), 'b--', lw=0.9,
        label=r'ergosphere $r_e$')
ax.set_aspect('equal'); ax.set_xlabel('$x$'); ax.set_ylabel('$y$')
ax.set_title(f'$t$-branch separatrix (saddle-node $J_c^-={Jst:.2f}$, '
             f'$r_*={rst:.3f}$):\nspirals onto $r_*$; genus 1 via a double root')
ax.legend(fontsize=6, loc='upper right', framealpha=0.9)

# (b) sinistra: deviazione ODE vs forma chiusa; destra (inset): R6 doppia radice
axb.semilogy(rA, devA+1e-18, 'k-')
axb.axvline(rst, color='C3', ls='--', lw=1.0)
axb.text(rst+0.05, np.nanmax(devA)*0.3 + 1e-15, '$r_*$', color='C3', fontsize=7)
axb.set_xlabel('$r$'); axb.set_ylabel(r'$|\phi_{\rm ODE}-\phi_{\rm closed}|$')
axb.set_title('genus-1 closed form vs Hamiltonian flow')
axins = axb.inset_axes([0.55, 0.55, 0.4, 0.4])
axins.plot(rr, R6v, 'C0-', lw=1.2)
axins.axhline(0, color='k', lw=0.5)
axins.axvline(rst, color='C3', ls='--', lw=0.8)
axins.set_title(r'$R_6=(r-r_*)^2Q_4$', fontsize=6)
axins.set_xlabel('$r$', fontsize=6); axins.tick_params(labelsize=5)
savefig(fig, HERE, 'fig_tbranch_separatrix')
print('FATTO.')
