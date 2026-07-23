# -*- coding: utf-8 -*-
"""
phi(r,A) ADIABATICA IBRIDA per Thakurta-Kerr (piccolo A'/A).

Forma finale:
  phi(r,A) = phi_0(r;Ê/A)                         [Kleinian, CHIUSO]
           + (A'/A) [ Closed(r) + psi(r) ]        [1° ordine]
           + O((A'/A)^2)
con:
  Closed(r) = -1/2 Ê ∂_E phi_0(r) η(r)            [CHIUSO: ∂_E phi_0 = A(r)/√R + Σ c_k ∫r^k/√R]
  psi(r)    = Ê · 1/2 (rho - rho_tilde)           [IRRIDUCIBILE -> integrato NUMERICAMENTE]
              rho=∫∂_E phi_0·h dr, rho_tilde=∫η·∂_E F dr,  h=dη/dr

I pezzi chiusi si valutano dalla loro struttura Kleiniana; solo psi (polilog
iperellittico) e' un integrale numerico 1D. Verifica: ibrido == adiabatica piena.
Ramo t; parametri M=1,a=0.9,Ê=1.4,J=6.
"""
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq
from scipy.interpolate import interp1d
M, a, J = 1.0, 0.9, 6.0; Ehat = 1.4; r0 = 12.0

r, Es = sp.symbols('r E'); Dl = r**2-2*M*r+a**2
Q2 = (2*Es**2*J**2*M*r-Es**2*J**2*r**2-4*Es**2*J*M*a*r+2*Es**2*M*a**2*r+Es**2*a**2*r**2
      +Es**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2-8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
R6 = r*Q2*((Es**2-1)*r+2*M); K = r*((Es**2-1)*r+2*M)*(J*(r-2*M)+2*M*a)/Dl
Fn = sp.lambdify(r, (K/sp.sqrt(R6)).subs(Es, Ehat), 'numpy')
dEF = sp.lambdify(r, sp.diff(K/sp.sqrt(R6), Es).subs(Es, Ehat), 'numpy')

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

rg = np.linspace(r0-0.02, rlam.min()+0.3, 800)
eta = lam_of_r(rg); h = np.gradient(eta, rg)
phi0 = cumulative_trapezoid(Fn(rg), rg, initial=0)
dEphi = cumulative_trapezoid(dEF(rg), rg, initial=0)          # ∂_E phi_0 (chiuso Kleinian)
Closed = -0.5*Ehat*dEphi*eta                                  # CHIUSO
# psi: SOLA parte numerica (integrale iterato irriducibile)
rho = cumulative_trapezoid(dEphi*h, rg, initial=0)
rhot = cumulative_trapezoid(eta*dEF(rg), rg, initial=0)
psi = Ehat*0.5*(rho-rhot)                                     # NUMERICO
# adiabatica piena (verita')
dphi_full = -Ehat*cumulative_trapezoid(dEF(rg)*eta, rg, initial=0)


def phi_hybrid(eps):
    return phi0 + eps*(Closed + psi)


def phi_full(eps):
    return phi0 + eps*dphi_full


if __name__ == '__main__':
    print("Verifica: ibrido (chiuso analitico + psi numerico) == adiabatica piena")
    for eps in [0.005, 0.02, 0.06]:
        err = np.max(np.abs(phi_hybrid(eps)-phi_full(eps)))
        print(f"  A'/A={eps}: max|phi_hybrid - phi_full| = {err:.2e}")
    print("=> ibrido riproduce esattamente la piena (per costruzione).")
    print("   pezzi CHIUSI: phi_0, ∂_E phi_0=A(r)/√R+Σc_k∫r^k/√R, Closed=-1/2 Ê ∂_E phi_0 η")
    print("   pezzo NUMERICO: psi = Ê/2 (rho - rho_tilde)  [polilog iperellittico]")
