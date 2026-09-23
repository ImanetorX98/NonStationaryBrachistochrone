# -*- coding: utf-8 -*-
"""
phi(r,A) ADIABATICA IBRIDA — RAMO tau (tempo proprio) per Thakurta-Kerr.

Stessa struttura del ramo t (kerr_adiabatic_phi_hybrid.py), ma:
  F_tau(r;E) = dphi_BL/dr = J r sqrt(wf) / (Δ sqrt(Δ - J^2 w))   [eq.56, BL]
  orologio di DERIVA  eta(r) = ∫ (d eta/dr) dr  lungo l'orbita tau,
                 d eta/dr = [E r^3 - 2MaJ r D_E/Δ] / sqrt(S)
  Il fattore conforme corre con eta, A = A(eta): la deriva delle cariche lungo
  l'orbita e' pesata dal tempo conforme TRASCORSO, qualunque sia il costo.
  Il tempo proprio tau(r) = ∫ r^2(r-2M)/sqrt(S) dr e' l'orologio che il funzionale
  misura (fissa il timing), NON la variabile lenta: una versione precedente usava
  tau(r) come peso, errore del ~40% sul coefficiente (vedi
  paper2/verification/verify_tau_branch_drift_clock.py).

Forma finale (identica al ramo t):
  phi(r,A) = phi_0(r;Ê/A) + (A'/A)[ Closed + psi ] + O((A'/A)^2)
  Closed = -1/2 Ê ∂_E phi_0 · eta        [CHIUSO]
  psi    = Ê/2 (rho - rho_tilde)         [NUMERICO, polilog iperellittico]

Verifica: ibrido == adiabatica piena. Ramo tau scattering |J|>J_c.
Parametri M=1,a=0.9,Ê=1.4,J=2.5 (scattering, J_c=a/Ê=0.643), r0=12.
"""
import numpy as np, sympy as sp
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
M, a, J = 1.0, 0.9, 2.5; Ehat = 1.4; r0 = 12.0

r, Es = sp.symbols('r E')
f = 1 - 2*M/r; Dl = r**2 - 2*M*r + a**2; w = Es**2 - f
beta = sp.sqrt(2*M*r/(r**2 + a**2)); k = r**2/(r**2 + a**2)
# F_tau = dphi_BL/dr  (eq.56)
Ftau = J*r*sp.sqrt(w*f)/(Dl*sp.sqrt(Dl - J**2*w))
# phi'_D (Doran) per il Lagrangiano proprio; L_tau = sqrt(Q/w)
phiD = a*beta/Dl + Ftau
Q = Dl*phiD**2 - 2*a*beta*phiD + k
Ltau = sp.sqrt(Q/w)                                   # dtau/dr (tempo proprio, solo confronto)
DEr = (Es**2 - 1)*r + 2*M
Sr = r*(r - 2*M)*DEr*(r*Dl - J**2*DEr)
Leta = (Es*r**3 - 2*M*a*J*r*DEr/Dl)/sp.sqrt(Sr)       # d eta/dr (orologio di deriva)

Fn = sp.lambdify(r, Ftau.subs(Es, Ehat), 'numpy')
dEF = sp.lambdify(r, sp.diff(Ftau, Es).subs(Es, Ehat), 'numpy')
Ln = sp.lambdify(r, Leta.subs(Es, Ehat), 'numpy')

# turning esterno: Δ - J^2 w = 0 (scattering, r_min > r_e)
wn = lambda rv: Ehat**2 - (1 - 2*M/rv)
Dn = lambda rv: rv**2 - 2*M*rv + a**2
rmin = brentq(lambda rv: Dn(rv) - J**2*wn(rv), 2.0 + 1e-9, r0)
print(f"J_c=a/Ê={a/Ehat:.4f}, J={J} (scattering), r_min={rmin:.5f}")

# backoff dal turning: F_tau~1/sqrt(Δ-J^2 w) diverge a r_min (integrabile ma
# trapezio lento). Griglia fitta, fermata a r_min+0.25.
rg = np.linspace(r0 - 0.02, rmin + 0.25, 4000)
eta = np.abs(cumulative_trapezoid(Ln(rg), rg, initial=0))   # eta trascorso >=0
h = -Ln(rg)                                           # d(eta)/dr ANALITICO
phi0 = cumulative_trapezoid(Fn(rg), rg, initial=0)
dEphi = cumulative_trapezoid(dEF(rg), rg, initial=0)  # ∂_E phi_0 (chiuso Kleinian)
Closed = -0.5*Ehat*dEphi*eta                          # CHIUSO
rho = cumulative_trapezoid(dEphi*h, rg, initial=0)
rhot = cumulative_trapezoid(eta*dEF(rg), rg, initial=0)
psi = Ehat*0.5*(rho - rhot)                           # NUMERICO
dphi_full = -Ehat*cumulative_trapezoid(dEF(rg)*eta, rg, initial=0)


def phi_hybrid(eps):
    return phi0 + eps*(Closed + psi)


def phi_full(eps):
    return phi0 + eps*dphi_full


if __name__ == '__main__':
    print("RAMO tau — Verifica: ibrido (chiuso + psi numerico) == adiabatica piena")
    for eps in [0.005, 0.02, 0.06]:
        err = np.max(np.abs(phi_hybrid(eps) - phi_full(eps)))
        print(f"  A'/A={eps}: max|phi_hybrid - phi_full| = {err:.2e}")
    print("=> ibrido riproduce la piena (identita' per parti; residuo = trapezio).")
    print("   CHIUSI: phi_0, ∂_E phi_0=A(r)/√R+Σc_k∫r^k/√R, Closed=-1/2 Ê ∂_E phi_0 η")
    print("   NUMERICO: psi = Ê/2 (rho - rho_tilde)  [polilog iperellittico]")
    print("   orologio di deriva η(r)=∫[E r^3-2MaJ r D_E/Δ]/√S dr lungo l'orbita tau")
