# Chiusura di FAY della terza specie d'ORIZZONTE in Vaidya.
# Lettera terza specie: eta_h = dr/(r-2m), polo SEMPLICE (residuo 2) al punto di
# Weierstrass r=2m + poli all'infinito. (NB: L_2m e' invece 2a specie -> zeta di Klein.)
# Claim (Fay, come L in Kerr ma centrato sull'orizzonte):
#   int eta_h = rho * log[ theta[delta](w-e_h) / theta[delta](w-e_inf) ] + holo,
#   delta ODD. Verifica round-trip: d/dr log-ratio = c0/(r-2m) + c1/sqrt(S) + c2 r/sqrt(S)
#   (holo = normalizzazione a-periodi), residuo lstsq ~ 0.
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np
from scipy.integrate import quad
import itertools

M, E, J = 1.0, 1.4, 2.5   # Vaidya frozen (a=0), m=1
# sestico S(r)=r(r-2m)DE[r^2(r-2m)-J^2 DE], DE=(E^2-1)r+2m  (coeff razionali)
Rr = PolynomialRing(QQ, ['x','Y']); x, Y = Rr.gens()
Em, Es, Js = QQ(1), QQ(7)/5, QQ(5)/2
DE = (Es**2-1)*x + 2*Em
S = Rr(x*(x-2*Em)*DE*(x**2*(x-2*Em) - Js**2*DE))
Spoly = S  # y^2 = S(x), sestico -> genere 2
X = SageRS(Y**2 - Spoly, prec=90)
omega = np.array(X.matrix_of_integral_values([Rr(1), x]), dtype=complex)[:, :2]  # a-per (1,x)/Y
tau = np.array(X.riemann_matrix(), dtype=complex)
ominv = np.linalg.inv(omega)

def Sn(rv):
    DEn = (E**2-1)*rv + 2*M
    return rv*(rv-2*M)*DEn*(rv**2*(rv-2*M) - J**2*DEn)
r_h = 2.0*M                      # orizzonte = branch point (base della mappa di Abel)

# Abel map w(r) = ominv @ [int_{r_h}^r 1/sqrt(S), int_{r_h}^r x/sqrt(S)]  (reg. al branch)
def Iu(r_to):   # int_{r_h}^{r_to} (1,x)/sqrt(S) dx, sub x=r_h+u^2 vicino al branch
    U = np.sqrt(abs(r_to - r_h)); sgn = 1.0 if r_to > r_h else -1.0
    g0 = lambda u: 2*u*(1.0)/np.sqrt(abs(Sn(r_h+sgn*u**2)))
    g1 = lambda u: 2*u*(r_h+sgn*u**2)/np.sqrt(abs(Sn(r_h+sgn*u**2)))
    I0 = quad(g0, 0, U, limit=200)[0]; I1 = quad(g1, 0, U, limit=200)[0]
    return sgn*np.array([I0, I1])
def w_of(r_to): return ominv @ Iu(r_to)
e_h = w_of(r_h + 1e-9)           # ~0 (base al branch orizzonte)
# immagine di infinito: int_{r_h}^{inf} (1,x)/sqrt(S) (converge: (1,x)/r^3)
def Iinf():
    g0 = lambda t: (1.0)/np.sqrt(abs(Sn(1.0/t)))*(1.0/t**2)   # x=1/t, dx=-dt/t^2
    g1 = lambda t: (1.0/t)/np.sqrt(abs(Sn(1.0/t)))*(1.0/t**2)
    # from r_h to inf: split r_h..R0 (direct) + R0..inf (sub)
    R0 = 40.0
    d0 = quad(lambda xx: 1.0/np.sqrt(abs(Sn(xx))), r_h+1e-9, R0, limit=300)[0]
    d1 = quad(lambda xx: xx/np.sqrt(abs(Sn(xx))), r_h+1e-9, R0, limit=300)[0]
    t0 = quad(g0, 1e-6, 1.0/R0, limit=300)[0]; t1 = quad(g1, 1e-6, 1.0/R0, limit=300)[0]
    return np.array([d0+t0, d1+t1])
e_inf = ominv @ Iinf()

# caratteristiche half-integer; odd = theta[delta](0)=0
def theta_d(z, av, bv, grad=False):
    av = np.array(av); bv = np.array(bv); zz = z + tau@av + bv
    pref = np.exp(1j*np.pi*(av@tau@av) + 2j*np.pi*(av@(z+bv)))
    th = complex(RiemannTheta(zz, tau))
    if not grad: return pref*th
    g = np.array([complex(RiemannTheta(zz, tau, derivs=[e])) for e in ([1,0],[0,1])])
    return 2j*np.pi*av + g/th   # grad_z log theta[delta]
halfs = [np.array(v)/2 for v in itertools.product([0,1], repeat=2)]
odd = [(av,bv) for av in halfs for bv in halfs if abs(theta_d(np.zeros(2),av,bv))<1e-6]
print("caratteristiche ODD:", len(odd))

# campiona orbita (tra orizzonte e r0), costruisci basi
rg = np.linspace(r_h+0.6, 11.0, 36)
inv_h = 1.0/(rg - 2*M)                       # terza specie letter
sq = np.sqrt(np.abs(Sn(rg)))
hol0 = 1.0/sq; hol1 = rg/sq                  # 2 olomorfi (2a specie base)
W = np.array([w_of(rv) for rv in rg])
dwdr = (ominv @ np.vstack([hol0, hol1])).T   # dw/dr = ominv @ (du_i/dr)
basis = np.vstack([inv_h, hol0, hol1]).T     # (N,3): {1/(r-2m), 1/sqrtS, r/sqrtS}

print("\n delta#  residuo(d log-ratio ~ c0/(r-2m)+holo)   c0")
best = None
for k,(av,bv) in enumerate(odd):
    D = np.array([ (theta_d(W[i]-e_h,av,bv,True)-theta_d(W[i]-e_inf,av,bv,True)) @ dwdr[i]
                   for i in range(len(rg))])
    coef,_,_,_ = np.linalg.lstsq(basis, D, rcond=None)
    resid = np.max(np.abs(basis@coef - D))/max(np.max(np.abs(D)),1e-30)
    print(f"  {k}: {resid:.3e}   c0={coef[0]:+.4f}")
    if best is None or resid<best[0]: best=(resid,k,coef)
print(f"\nMIGLIORE: delta #{best[1]}, residuo rel {best[0]:.2e}, c0={best[2][0]:.4f}")
print("=> se residuo ~0: eta_h (3a specie orizzonte) = log(theta-ratio) di FAY, CONFERMATO.")
