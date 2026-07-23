# Validazione decomposizione psi = psi_a (2a specie, zeta) + psi_b (3a specie, dilog)
# vs integrazione diretta -Ehat int dEF*eta dr, per VARI valori di innesco r0.
# Ramo tau. Conferma che psi_b (il pezzo dilog) e' cio' che resta dalla 3a specie.
import numpy as np, sympy as sp
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
M, a, J, Ehat = 1.0, 0.9, 2.5, 1.4

r, Es = sp.symbols('r E'); f = 1-2*M/r; Dl = r**2-2*M*r+a**2; w = Es**2-f
Emu = (Es**2-1)*r+2*M
Ftau = J*r*(r-2*M)*Emu/(Dl*sp.sqrt(r*(r-2*M)*Emu*(r*Dl-J**2*Emu)))  # =K_t/sqrt(S)
S = r*(r-2*M)*Emu*(r*Dl-J**2*Emu)
Fn = sp.lambdify(r, Ftau.subs(Es, Ehat), 'numpy')
dEF = sp.lambdify(r, sp.diff(Ftau, Es).subs(Es, Ehat), 'numpy')
Sn = sp.lambdify(r, S.subs(Es, Ehat), 'numpy')
wn = lambda rv: Ehat**2-(1-2*M/rv)
rmin = brentq(lambda rv: (rv**2-2*M*rv+a**2)-J**2*wn(rv), 2.0+1e-9, 20)
print(f"rmin={rmin:.4f}  (J={J}, Ehat={Ehat})")
# clock split: deta/dr = (r^3 - 2M r^2)/sqrt(S). h_a=r^3/sqrtS (2a), h_b=-2M r^2/sqrtS (3a)
def ha(rg): return rg**3/np.sqrt(Sn(rg))
def hb(rg): return -2*M*rg**2/np.sqrt(Sn(rg))

print("\n r0    delta_phi_direct   psi_a(2a/zeta)   psi_b(3a/dilog)  |psi_b|/|psi|")
for r0 in [14.0, 12.0, 10.0, 8.0, 6.5]:
    rg = np.linspace(r0-0.02, rmin+0.25, 4000)
    dEphi = cumulative_trapezoid(dEF(rg), rg, initial=0)          # partial_E phi_0 = A
    ha_v, hb_v = ha(rg), hb(rg)                                   # segnati (dη/dr)
    eta_a = cumulative_trapezoid(ha_v, rg, initial=0)             # segnato
    eta_b = cumulative_trapezoid(hb_v, rg, initial=0)             # segnato
    eta = eta_a + eta_b
    # psi = 1/2 Ehat (rho - rho~),  rho=int dEphi*h dr, rho~=int eta*dEF dr
    def psi_piece(hv, etapart):
        rho = cumulative_trapezoid(dEphi*hv, rg, initial=0)
        rhot = cumulative_trapezoid(etapart*dEF(rg), rg, initial=0)
        return 0.5*Ehat*(rho-rhot)
    psi_a = psi_piece(ha_v, eta_a); psi_b = psi_piece(hb_v, eta_b)
    psi = psi_a + psi_b
    # delta_phi_direct (coeff O(A'/A)) = -Ehat int dEF*eta dr  (verita' ODE)
    dphi_dir = -Ehat*cumulative_trapezoid(dEF(rg)*eta, rg, initial=0)
    Closed = -0.5*Ehat*dEphi*eta
    assembled = Closed + psi
    err = np.max(np.abs(assembled - dphi_dir))
    wpsi = np.max(np.abs(psi_b))/max(np.max(np.abs(psi)),1e-30)
    print(f"{r0:5.1f}  match_err={err:.2e}   {np.max(np.abs(psi_a)):.4f}      "
          f"{np.max(np.abs(psi_b)):.4f}       {wpsi:.1%}")
print("\n=> assembled (Closed+psi_a+psi_b) == diretto ODE (linearita' verificata);")
print("   psi_b = pezzo 3a specie (dilog iperellittico), peso riportato per innesco.")
