# -*- coding: utf-8 -*-
"""
La curva adiabatica dai pezzi ANALITICI CHIUSI vs la curva completa.
Mostra quanto pesa la parte irriducibile (integrale iterato) che non sappiamo
integrare in forma chiusa.

phi = phi_0 + (A'/A)[ closed + irreducible ] + O((A'/A)^2)
  closed      = -1/2 Ê ∂_E phi_0 · η        (Kleinian, esplicito)
  irreducible = Ê ψ = Ê·1/2(ρ-ρ̃)           (polilog iperellittico, NON chiuso)

Ramo t; parametri M=1,a=0.9,Ê=1.4,J=6; A'/A=eps.
"""
import os, sys
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq
from scipy.interpolate import interp1d
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a, J = 1.0, 0.9, 6.0
Ehat = 1.4
r0 = 12.0
eps = 0.06                       # A'/A (rate), moderato per visibilita'

# forma chiusa F=K/√R6 e ∂_E F
r, Es = sp.symbols('r E'); Dl = r**2-2*M*r+a**2
Q2 = (2*Es**2*J**2*M*r-Es**2*J**2*r**2-4*Es**2*J*M*a*r+2*Es**2*M*a**2*r+Es**2*a**2*r**2
      +Es**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2-8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
R6 = r*Q2*((Es**2-1)*r+2*M); K = r*((Es**2-1)*r+2*M)*(J*(r-2*M)+2*M*a)/Dl
Fn = sp.lambdify(r, (K/sp.sqrt(R6)).subs(Es, Ehat), 'numpy')
dEF = sp.lambdify(r, sp.diff(K/sp.sqrt(R6), Es).subs(Es, Ehat), 'numpy')

# flusso congelato -> eta(r)=t(r), h
rr, pr, Ess, Jss = sp.symbols('r pr E J_')
f2 = 1-2*M/rr; Dl2 = rr**2-2*M*rr+a**2; b2 = 2*M*a/rr; v2 = 1-f2/Ess**2
P2 = rr**2+a**2+2*M*a**2/rr; Pb2 = P2+b2**2/Ess**2
H2 = Jss*b2*v2/Pb2+sp.sqrt(Dl2*v2/Pb2)*sp.sqrt((Dl2/rr**2)*pr**2+Jss**2/Pb2)-1
H2n = sp.lambdify((rr, pr, Ess, Jss), H2, 'numpy')
dHp = sp.lambdify((rr, pr, Ess, Jss), sp.diff(H2, pr), 'numpy')
dHr = sp.lambdify((rr, pr, Ess, Jss), sp.diff(H2, rr), 'numpy')
def prof(rv):
    pg = np.linspace(-80, 80, 3001); Hv = H2n(rv, pg, Ehat, J)
    rts = [brentq(lambda p: H2n(rv, p, Ehat, J), pg[i], pg[i+1]) for i in range(len(pg)-1)
           if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1] < 0]
    ing = [p for p in rts if dHp(rv, p, Ehat, J) < 0]; return min(ing) if ing else np.nan
ev = lambda lam, y: y[1]; ev.terminal = True; ev.direction = 1
s = solve_ivp(lambda lam, y: [dHp(y[0], y[1], Ehat, J), -dHr(y[0], y[1], Ehat, J)],
              [0, 300], [r0, prof(r0)], rtol=1e-11, atol=1e-13, max_step=0.01,
              dense_output=True, events=ev)
lam = np.linspace(0, s.t[-1], 6000); rlam = s.sol(lam)[0]; lam_of_r = interp1d(rlam, lam)

rg = np.linspace(r0-0.02, rlam.min()+0.3, 700)
eta = lam_of_r(rg); h = np.gradient(eta, rg); dEFv = dEF(rg)
phi0 = cumulative_trapezoid(Fn(rg), rg, initial=0)            # phi_0(r) (leading)
dEphi = cumulative_trapezoid(dEFv, rg, initial=0)             # ∂_E phi_0
rho = cumulative_trapezoid(dEphi*h, rg, initial=0)
rhot = cumulative_trapezoid(eta*dEFv, rg, initial=0)
psi = 0.5*(rho-rhot)
dphi_tot = -Ehat*cumulative_trapezoid(dEFv*eta, rg, initial=0)  # δφ_tot (coeff)
Closed = -0.5*Ehat*dEphi*eta                                    # parte chiusa
Irr = Ehat*psi                                                  # parte irriducibile

# curve (x,y) ramo incoming
phiC = phi0 + eps*Closed              # phi_0 + correzione chiusa
phiF = phi0 + eps*dphi_tot            # completa (chiuso + irriducibile)

fig, (ax, axb) = plt.subplots(1, 2, figsize=(2*COL, COL*0.95))
th = np.linspace(0, 2*np.pi, 200)
ax.plot(rg*np.cos(phi0), rg*np.sin(phi0), 'C7-', lw=1.0, label=r'$\phi_0$ (frozen)')
ax.plot(rg*np.cos(phiC), rg*np.sin(phiC), 'C0--', lw=1.8,
        label=r'$\phi_0+$closed correction')
ax.plot(rg*np.cos(phiF), rg*np.sin(phiF), 'C3:', lw=1.8,
        label=r'full (closed$+$irreducible)')
ax.plot(2*M*np.cos(th), 2*M*np.sin(th), 'b:', lw=0.6)
ax.plot(r0, 0, 'ks', ms=3)
ax.set_aspect('equal'); ax.set_xlabel('$x$'); ax.set_ylabel('$y$')
ax.set_title(f'incoming branch ($A^\\prime/A={eps}$): closed pieces\n'
             'vs full; gap $=$ irreducible', fontsize=6.6)
ax.legend(fontsize=5.4, loc='upper left', framealpha=0.9)

# pannello b: le due correzioni e lo scarto (parte irriducibile)
axb.plot(rg, eps*Closed, 'C0-', lw=1.5, label=r'closed correction $\epsilon\,$Closed')
axb.plot(rg, eps*dphi_tot, 'k-', lw=1.0, label=r'total $\epsilon\,\delta\phi_{\rm tot}$')
axb.fill_between(rg, eps*Closed, eps*dphi_tot, color='C3', alpha=0.4,
                 label=r'irreducible gap $\epsilon\,$Irr')
axb.set_xlabel('$r$'); axb.set_ylabel(r'$\Delta\varphi$ correction (rad)')
axb.set_title('correction breakdown:\nclosed (74%) + irreducible (26%)', fontsize=6.6)
axb.legend(fontsize=5.6, loc='upper right', framealpha=0.9)
axb.invert_xaxis()
savefig(fig, HERE, 'fig_adiabatic_pieces')
fr = np.sqrt(np.sum(Irr**2*abs(np.gradient(rg))))/np.sqrt(np.sum(dphi_tot**2*abs(np.gradient(rg))))
print(f"A'/A={eps}: frazione irriducibile (L2) = {fr:.1%}")
print(f"scarto max |phiF-phiC| = {np.max(np.abs(phiF-phiC)):.3f} rad "
      f"(su phi_0 max {np.max(np.abs(phi0)):.2f} rad)")
print('FATTO.')
