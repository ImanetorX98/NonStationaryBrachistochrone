# -*- coding: utf-8 -*-
"""
HYBRID ADIABATIC phi(r,A) -- tau BRANCH (proper time) for Thakurta-Kerr.

Same structure as the t branch (kerr_adiabatic_phi_hybrid.py), but:
  F_tau(r;E) = dphi_BL/dr = J r sqrt(wf) / (Δ sqrt(Δ - J^2 w))   [eq.56, BL]
  DRIFT clock  eta(r) = ∫ (d eta/dr) dr  along the tau orbit,
                 d eta/dr = [E r^3 - 2MaJ r D_E/Δ] / sqrt(S)
  The conformal factor runs with eta, A = A(eta): the drift of the charges along
  the orbit is weighted by the ELAPSED conformal time, whatever the cost.
  The proper time tau(r) = ∫ r^2(r-2M)/sqrt(S) dr is the clock that the functional
  measures (it fixes the timing), NOT the slow variable: an earlier version used
  tau(r) as the weight, a ~40% error in the coefficient (see
  paper2/verification/verify_tau_branch_drift_clock.py).

Final form (identical to the t branch):
  phi(r,A) = phi_0(r;Ê/A) + (A'/A)[ Closed + psi ] + O((A'/A)^2)
  Closed = -1/2 Ê ∂_E phi_0 · eta        [CLOSED]
  psi    = Ê/2 (rho - rho_tilde)         [NUMERICAL, hyperelliptic polylog]

Check: hybrid == full adiabatic solution. Scattering tau branch |J|>J_c.
Parameters M=1,a=0.9,Ê=1.4,J=2.5 (scattering, J_c=a/Ê=0.643), r0=12.
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
# phi'_D (Doran) for the proper-time Lagrangian; L_tau = sqrt(Q/w)
phiD = a*beta/Dl + Ftau
Q = Dl*phiD**2 - 2*a*beta*phiD + k
Ltau = sp.sqrt(Q/w)                                   # dtau/dr (proper time, comparison only)
DEr = (Es**2 - 1)*r + 2*M
Sr = r*(r - 2*M)*DEr*(r*Dl - J**2*DEr)
Leta = (Es*r**3 - 2*M*a*J*r*DEr/Dl)/sp.sqrt(Sr)       # d eta/dr (drift clock)

Fn = sp.lambdify(r, Ftau.subs(Es, Ehat), 'numpy')
dEF = sp.lambdify(r, sp.diff(Ftau, Es).subs(Es, Ehat), 'numpy')
Ln = sp.lambdify(r, Leta.subs(Es, Ehat), 'numpy')

# outer turning point: Δ - J^2 w = 0 (scattering, r_min > r_e)
wn = lambda rv: Ehat**2 - (1 - 2*M/rv)
Dn = lambda rv: rv**2 - 2*M*rv + a**2
rmin = brentq(lambda rv: Dn(rv) - J**2*wn(rv), 2.0 + 1e-9, r0)
print(f"J_c=a/Ê={a/Ehat:.4f}, J={J} (scattering), r_min={rmin:.5f}")

# back off from the turning point: F_tau~1/sqrt(Δ-J^2 w) diverges at r_min (integrable but
# the trapezoidal rule converges slowly). Fine grid, stopped at r_min+0.25.
rg = np.linspace(r0 - 0.02, rmin + 0.25, 4000)
eta = np.abs(cumulative_trapezoid(Ln(rg), rg, initial=0))   # elapsed eta >= 0
h = -Ln(rg)                                           # ANALYTIC d(eta)/dr
phi0 = cumulative_trapezoid(Fn(rg), rg, initial=0)
dEphi = cumulative_trapezoid(dEF(rg), rg, initial=0)  # ∂_E phi_0 (closed, Kleinian)
Closed = -0.5*Ehat*dEphi*eta                          # CLOSED
rho = cumulative_trapezoid(dEphi*h, rg, initial=0)
rhot = cumulative_trapezoid(eta*dEF(rg), rg, initial=0)
psi = Ehat*0.5*(rho - rhot)                           # NUMERICAL
dphi_full = -Ehat*cumulative_trapezoid(dEF(rg)*eta, rg, initial=0)


def phi_hybrid(eps):
    return phi0 + eps*(Closed + psi)


def phi_full(eps):
    return phi0 + eps*dphi_full


if __name__ == '__main__':
    print("tau BRANCH -- Check: hybrid (closed + numerical psi) == full adiabatic solution")
    for eps in [0.005, 0.02, 0.06]:
        err = np.max(np.abs(phi_hybrid(eps) - phi_full(eps)))
        print(f"  A'/A={eps}: max|phi_hybrid - phi_full| = {err:.2e}")
    print("=> the hybrid reproduces the full solution (identity by parts; residual = trapezoidal rule).")
    print("   CLOSED: phi_0, ∂_E phi_0=A(r)/√R+Σc_k∫r^k/√R, Closed=-1/2 Ê ∂_E phi_0 η")
    print("   NUMERICAL: psi = Ê/2 (rho - rho_tilde)  [hyperelliptic polylog]")
    print("   drift clock η(r)=∫[E r^3-2MaJ r D_E/Δ]/√S dr along the tau orbit")
