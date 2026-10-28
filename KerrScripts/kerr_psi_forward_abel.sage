# -*- coding: utf-8 -*-
# BUILDING BLOCK, level 3 (corrected): FORWARD Abel map u(r) on the tau sextic +
# theta(u(r)) with derivatives (ingredients of the genus-2 polylog kernel for psi).
#
# ROBUST toolchain (works around the Singular bug in abelfunctions.integralbasis):
#   - periods tau, A|B  -> Sage RiemannSurface (cohomology basis (1,x))
#   - Abel map u(r)     -> direct integration of (1,x)/(2 sqrt(S)) (convention f_y=2y)
#                          normalised: u = A^{-1} I(r)
#   - theta + derivatives -> RiemannTheta from abelfunctions (numpy, any tau)
# psi(r)=1/2 Ehat (rho-rho~) becomes a genus-2 polylog in the coordinate u(r).
# RATIONAL parameters: M=1, a=9/10, E=7/5, J=5/2 (scattering, real orbit r>=r_min).
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np
from scipy.integrate import quad

M, a, E, J = QQ(1), QQ(9)/10, QQ(7)/5, QQ(5)/2
Mf, af, Ef, Jf = 1.0, 0.9, 1.4, 2.5

Rr = PolynomialRing(QQ, ['x', 'y']); x, y = Rr.gens()
Dl = x**2 - 2*M*x + a**2; Emu = (E**2 - 1)*x + 2*M
Sint = (1250*x*(x - 2*M)*Emu*(x*Dl - J**2*Emu)).change_ring(QQ)   # integral model
print('Sint =', Sint)

X = SageRS(y**2 - Sint, prec=80)
tau = np.array(X.riemann_matrix(), dtype=complex)
PM = np.array(X.period_matrix(), dtype=complex)          # 2x4 = [A | B]
Amat, Bmat = PM[:, :2], PM[:, 2:]
print('genus =', X.genus, ' tau sym err =', np.linalg.norm(tau - tau.T))
print('||A^{-1}B - tau|| =', np.linalg.norm(np.linalg.solve(Amat, Bmat) - tau))

# numerical sextic (same integral model) and differentials (1,x)/(2 sqrt(Sint))
cS = [float(c) for c in Sint.polynomial(x).list()]      # coeff crescenti
def Snum(r): return sum(c*r**k for k, c in enumerate(cS))
def w0(r): return 1.0/(2.0*np.sqrt(Snum(r)))            # dx/(2y)
def w1(r): return r/(2.0*np.sqrt(Snum(r)))              # x dx/(2y)

r0 = 12.0
rmin = 4.046197656444178                                # turning point (branch point)
def abel_unnorm(r):    # I(r)=int_{r0}^{r}(w0,w1)dx  (real axis, S>0 on the orbit)
    I0 = quad(w0, r0, r, limit=200)[0]
    I1 = quad(w1, r0, r, limit=200)[0]
    return np.array([I0, I1], dtype=complex)
def u_of_r(r):         # normalised Abel map
    return np.linalg.solve(Amat, abel_unnorm(r))

# theta and derivatives along the orbit
print("\n r      u1            u2           |theta|      |grad|       |Hess|")
for r in [12.0, 10.0, 8.0, 6.0, 5.0, 4.3]:
    u = u_of_r(r)
    th = complex(RiemannTheta(u, tau))
    g = [complex(RiemannTheta(u, tau, derivs=[e])) for e in ([1,0],[0,1])]
    H = [[complex(RiemannTheta(u, tau, derivs=[p,q])) for q in ([1,0],[0,1])]
         for p in ([1,0],[0,1])]
    gn = np.linalg.norm(g); Hn = np.linalg.norm(np.array(H))
    print(f"{r:5.1f}  {u[0]:+.4f}  {u[1]:+.4f}  {abs(th):.4e}  {gn:.3e}  {Hn:.3e}")

print("\nOK: forward coordinate u(r) + theta/grad/Hess available.")
print("=> ingredients of the genus-2 Kronecker-Eisenstein kernel for psi(r) ready.")
print("   (rho-rho~ = length-2 iterated integral in the 1-forms dEF, L; the")
print("    kernel g^(n)(u,tau) is built from theta and the derivatives evaluated here.)")
