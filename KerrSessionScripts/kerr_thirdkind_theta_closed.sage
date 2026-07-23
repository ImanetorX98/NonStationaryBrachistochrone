# VALIDAZIONE forma chiusa theta del 3a-specie vs diretto (livello differenziale).
# psi_b usa eta_b = int h_b, h_b = -2M r^2/sqrt(S) (dipolo 3a specie).
# Claim: eta_b(r) = rho0 * log[ theta[delta](w-e+) / theta[delta](w-e-) ] + holo,
#   delta = caratteristica ODD (theta[delta](0)=0). Valido:
#   d/dr log-ratio = c0 * h_b + c1 du1/dr + c2 du2/dr   (holo = normalizz. a-periodi).
# Modello dispari (quintica), w normalizzato, base al branch point s_b=1/r_min -> e-=-e+.
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np
from scipy.integrate import quad
import itertools

M, a, J, Ehat = 1.0, 0.9, 2.5, 1.4
Rs = PolynomialRing(QQ, ['s', 'Y']); s, Y = Rs.gens()
lam = [QQ(1200), QQ(-2300), QQ(-11428), QQ(-5519), QQ(24700), QQ(62500)]
qpoly = sum(lam[i]*s**i for i in range(6))
X = SageRS(Y**2 - qpoly, prec=80)
omega = np.array(X.matrix_of_integral_values([Rs(1), s]), dtype=complex)[:, :2]  # a-per 1a specie
tau = np.array(X.riemann_matrix(), dtype=complex)
ominv = np.linalg.inv(omega)

def qn(sv): return sum(float(lam[i])*sv**i for i in range(6))
rmin = 4.046197656444178; s_b = 1.0/rmin                        # branch point base
def Iu(s_to):   # int_{s_to}^{s_b} (1,s)/sqrt(q) ds regolarizzato s=s_b-u^2
    U = np.sqrt(s_b - s_to)
    g0 = lambda u: 2*u*(1.0)/np.sqrt(qn(s_b-u**2))
    g1 = lambda u: 2*u*(s_b-u**2)/np.sqrt(qn(s_b-u**2))
    I0 = quad(g0, 0, U, limit=200)[0]; I1 = quad(g1, 0, U, limit=200)[0]
    return -np.array([I0, I1])                                  # segno per s_b->s_to
def w_of(s_to): return ominv @ Iu(s_to)
e_plus = w_of(0.0); e_minus = -e_plus                           # base branch pt -> e-=-e+

# caratteristiche: 16 half-integer; odd = theta[delta](0)=0
def theta_d(z, av, bv, derivs=None):
    av = np.array(av); bv = np.array(bv)
    zz = z + tau @ av + bv
    pref = np.exp(1j*np.pi*(av @ tau @ av) + 2j*np.pi*(av @ (z + bv)))
    if derivs is None:
        return pref*complex(RiemannTheta(zz, tau))
    # grad_z log theta[delta] = 2pi i a + grad theta(zz)/theta(zz)
    th = complex(RiemannTheta(zz, tau))
    g = np.array([complex(RiemannTheta(zz, tau, derivs=[e])) for e in ([1,0],[0,1])])
    return 2j*np.pi*av + g/th
halfs = [np.array(v)/2 for v in itertools.product([0,1], repeat=2)]
odd = []
for av in halfs:
    for bv in halfs:
        if abs(theta_d(np.zeros(2), av, bv)) < 1e-6:
            odd.append((av, bv))
print("caratteristiche ODD trovate:", len(odd))

# campiona orbita, costruisci h_b(r), du_i/dr(r)
def Sn(rv):
    f=1-2*M/rv; Dl=rv**2-2*M*rv+a**2; Emu=(Ehat**2-1)*rv+2*M
    return rv*(rv-2*M)*Emu*(rv*Dl-J**2*Emu)
rg = np.linspace(11.5, rmin+0.35, 40)
h_b = -2*M*rg**2/np.sqrt(Sn(rg))
sr = 1.0/rg; dsdr = -1.0/rg**2
du1dr = (1/np.sqrt(qn(sr)))*dsdr; du2dr = (sr/np.sqrt(qn(sr)))*dsdr
W = np.array([w_of(sv) for sv in sr])                          # w(r) (2,) each
dwdr = (ominv @ np.vstack([du1dr, du2dr])).T                   # (N,2)

# base completa 1a+2a+3a specie in x-coord: {r^k/sqrt(S)}, k=0..3
sqrtS = np.sqrt(Sn(rg))
basis = np.vstack([rg**k/sqrtS for k in range(4)]).T            # (N,4)
print("\n delta_id   fit_residuo(D ~ sum_k c_k r^k/sqrtS)   c2(coeff 3a specie)")
best=None
for k,(av,bv) in enumerate(odd):
    D = np.array([ (theta_d(W[i]-e_plus,av,bv,1)-theta_d(W[i]-e_minus,av,bv,1)) @ dwdr[i]
                   for i in range(len(rg))])
    coef,_,_,_ = np.linalg.lstsq(basis, D, rcond=None)
    resid = np.max(np.abs(basis@coef - D))/max(np.max(np.abs(D)),1e-30)
    print(f"  {k}: {resid:.3e}    c2={coef[2]:+.4f}")
    if best is None or resid<best[0]: best=(resid,k,coef)
print(f"\nMIGLIORE: delta #{best[1]}, residuo rel {best[0]:.2e}, c0={best[2][0]:.4f}")
print("=> se residuo ~0: eta_b (3a specie) = forma chiusa theta[delta] CONFERMATA.")
print("   rho0 atteso ~1.06 (residuo 3a specie); c0 lo riflette (a meno di norm.)")
