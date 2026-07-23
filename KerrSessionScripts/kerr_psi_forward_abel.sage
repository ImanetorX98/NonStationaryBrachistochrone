# -*- coding: utf-8 -*-
# MATTONE livello 3 (corretto): Abel map IN AVANTI u(r) sulla sestica tau +
# theta(u(r)) con derivate (ingredienti kernel polilog genus-2 per psi).
#
# Toolchain ROBUSTO (aggira il bug Singular di abelfunctions.integralbasis):
#   - periodi tau, A|B  -> RiemannSurface di Sage (base coomologia (1,x))
#   - Abel map u(r)     -> integrazione diretta di (1,x)/(2 sqrt(S)) (conv. f_y=2y)
#                          normalizzata: u = A^{-1} I(r)
#   - theta + derivate  -> RiemannTheta di abelfunctions (numpy, ogni tau)
# psi(r)=1/2 Ehat (rho-rho~) diventa polilog genus-2 nella coordinata u(r).
# Params RAZIONALI: M=1, a=9/10, E=7/5, J=5/2 (scattering, orbita reale r>=r_min).
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np
from scipy.integrate import quad

M, a, E, J = QQ(1), QQ(9)/10, QQ(7)/5, QQ(5)/2
Mf, af, Ef, Jf = 1.0, 0.9, 1.4, 2.5

Rr = PolynomialRing(QQ, ['x', 'y']); x, y = Rr.gens()
Dl = x**2 - 2*M*x + a**2; Emu = (E**2 - 1)*x + 2*M
Sint = (1250*x*(x - 2*M)*Emu*(x*Dl - J**2*Emu)).change_ring(QQ)   # modello intero
print('Sint =', Sint)

X = SageRS(y**2 - Sint, prec=80)
tau = np.array(X.riemann_matrix(), dtype=complex)
PM = np.array(X.period_matrix(), dtype=complex)          # 2x4 = [A | B]
Amat, Bmat = PM[:, :2], PM[:, 2:]
print('genus =', X.genus, ' tau sym err =', np.linalg.norm(tau - tau.T))
print('||A^{-1}B - tau|| =', np.linalg.norm(np.linalg.solve(Amat, Bmat) - tau))

# sestica numerica (stesso modello intero) e differenziali (1,x)/(2 sqrt(Sint))
cS = [float(c) for c in Sint.polynomial(x).list()]      # coeff crescenti
def Snum(r): return sum(c*r**k for k, c in enumerate(cS))
def w0(r): return 1.0/(2.0*np.sqrt(Snum(r)))            # dx/(2y)
def w1(r): return r/(2.0*np.sqrt(Snum(r)))              # x dx/(2y)

r0 = 12.0
rmin = 4.046197656444178                                # turning (branch point)
def abel_unnorm(r):    # I(r)=int_{r0}^{r}(w0,w1)dx  (asse reale, S>0 su orbita)
    I0 = quad(w0, r0, r, limit=200)[0]
    I1 = quad(w1, r0, r, limit=200)[0]
    return np.array([I0, I1], dtype=complex)
def u_of_r(r):         # Abel map normalizzata
    return np.linalg.solve(Amat, abel_unnorm(r))

# theta e derivate lungo l'orbita
print("\n r      u1            u2           |theta|      |grad|       |Hess|")
for r in [12.0, 10.0, 8.0, 6.0, 5.0, 4.3]:
    u = u_of_r(r)
    th = complex(RiemannTheta(u, tau))
    g = [complex(RiemannTheta(u, tau, derivs=[e])) for e in ([1,0],[0,1])]
    H = [[complex(RiemannTheta(u, tau, derivs=[p,q])) for q in ([1,0],[0,1])]
         for p in ([1,0],[0,1])]
    gn = np.linalg.norm(g); Hn = np.linalg.norm(np.array(H))
    print(f"{r:5.1f}  {u[0]:+.4f}  {u[1]:+.4f}  {abs(th):.4e}  {gn:.3e}  {Hn:.3e}")

print("\nOK: coordinata u(r) in avanti + theta/grad/Hess disponibili.")
print("=> ingredienti kernel Kronecker-Eisenstein genus-2 per psi(r) pronti.")
print("   (rho-rho~ = integrale iterato lunghezza-2 nelle 1-forme dEF, L; il")
print("    kernel g^(n)(u,tau) si costruisce da theta e derivate qui valutate.)")
