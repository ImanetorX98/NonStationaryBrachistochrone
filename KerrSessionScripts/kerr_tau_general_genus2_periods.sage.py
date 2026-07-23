import time
t0=time.time()
from sage.all import QQ, ZZ, PolynomialRing, lcm, expand
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface
def log(m): print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)

Rx=PolynomialRing(QQ,'x'); x=Rx.gen()
M=QQ(1); a=QQ(9)/10; E=QQ(6)/5; Jc=a/E; J=-QQ(9)/10*Jc   # -0.9 Jc = -27/40
log(f"J = -0.9 Jc = {J}")
num_wf  = (x-2*M)*((E**2-1)*x+2*M)
num_DJw = x*(x**2-2*M*x+a**2) - J**2*((E**2-1)*x+2*M)
sextic  = expand(x*num_wf*num_DJw)
log(f"sextic = {sextic}")
# radici
from sage.all import CDF
rts=sextic.roots(CDF, multiplicities=False)
log("branch points:")
for z in sorted(rts, key=lambda z: z.real()):
    print(f"     {complex(z).real:+.4f} {complex(z).imag:+.4f}i")
# modello intero
L=lcm([c.denominator() for c in sextic.coefficients()])
R6i=expand(L*sextic)
log(f"L={L}; intero? {all(c in ZZ for c in R6i.coefficients())}")
Rxy=PolynomialRing(QQ,['X','Y']); X,Y=Rxy.gens()
f=Y**2-R6i(X)
log("build RiemannSurface prec=40 ...")
S=RiemannSurface(f, prec=40)
log(f"  genere = {S.genus}")
tau=S.riemann_matrix()
log("matrice di Riemann tau =")
for row in tau:
    print("   " + "  ".join(f"{complex(z).real:+.5f}{complex(z).imag:+.5f}i" for z in row))
sym=(tau-tau.transpose()).norm()
from sage.all import RealField, matrix
ImT=matrix(RealField(40),2,2,[complex(tau[i][j]).imag for i in range(2) for j in range(2)])
log(f"|tau-tau^T| = {float(sym):.2e} (simmetrica); autovalori Im(tau) = {[float(e) for e in ImT.eigenvalues()]} (>0)")
log("DONE.")
