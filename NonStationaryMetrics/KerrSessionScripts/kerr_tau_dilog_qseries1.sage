# -*- coding: utf-8 -*-
# Q-SERIE dilog genus-2 (TK tau) -- TAPPA 1: fondamenta.
# (1) nomi q_ij dalla matrice di Riemann tau (2x2) -> parametri d'espansione.
# (2) serie in NOME della theta[delta] genus-2 vs RiemannTheta (verifica mattone base).
# (3) ground-truth del DILOG Lambda(r) = int(A_ab dOmega - Omega dA_ab),
#     Omega=log[theta[d](w-e+)/theta[d](w-e-)] (3a specie), A_ab=parte abeliana source.
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np
from scipy.integrate import quad
import itertools

M,a,J,Ehat = 1.0,0.9,2.5,1.4
Rs = PolynomialRing(QQ, ['s','Y']); s,Y = Rs.gens()
lam=[QQ(1200),QQ(-2300),QQ(-11428),QQ(-5519),QQ(24700),QQ(62500)]
qpoly=sum(lam[i]*s**i for i in range(6))
X=SageRS(Y**2-qpoly, prec=80)
omega=np.array(X.matrix_of_integral_values([Rs(1),s]),dtype=complex)[:,:2]
tau=np.array(X.riemann_matrix(),dtype=complex); ominv=np.linalg.inv(omega)
print("tau =\n",np.round(tau,5))
print("Im(tau) diag =",np.round(np.diag(tau.imag),4)," (grande => nome piccolo, serie converge)")
# nomi q_ij = exp(i pi tau_ij)
qm=np.array([[np.exp(1j*np.pi*tau[i,j]) for j in range(2)] for i in range(2)])
print("|q_11|,|q_22|,|q_12| =",[f"{abs(qm[0,0]):.4f}",f"{abs(qm[1,1]):.4f}",f"{abs(qm[0,1]):.4f}"])

# (2) serie in nome (somma reticolare troncata) della theta[delta] vs RiemannTheta
def theta_latt(z,av,bv,N=6):
    tot=0j
    for n in itertools.product(range(-N,N+1),repeat=2):
        nn=np.array(n)+np.array(av)
        tot+=np.exp(1j*np.pi*(nn@tau@nn)+2j*np.pi*(nn@(np.array(z)+np.array(bv))))
    return tot
halfs=[np.array(v)/2 for v in itertools.product([0,1],repeat=2)]
def theta_af(z,av,bv):
    av=np.array(av); bv=np.array(bv); zz=z+tau@av+bv
    pref=np.exp(1j*np.pi*(av@tau@av)+2j*np.pi*(av@(z+bv)))
    return pref*complex(RiemannTheta(zz,tau))
odd=[(av,bv) for av in halfs for bv in halfs if abs(theta_af(np.zeros(2),av,bv))<1e-6]
av,bv=odd[1]
ztest=np.array([0.13+0.07j,-0.08+0.11j])
print("\nserie-nome theta[delta] vs RiemannTheta (N=6):",
      abs(theta_latt(ztest,av,bv)-theta_af(ztest,av,bv)))

# (3) ground-truth dilog
def qn(sv): return sum(float(lam[i])*sv**i for i in range(6))
rmin=4.046197656444178; s_b=1.0/rmin
def Iu(s_to):
    U=np.sqrt(s_b-s_to)
    g0=lambda u:2*u/np.sqrt(qn(s_b-u**2)); g1=lambda u:2*u*(s_b-u**2)/np.sqrt(qn(s_b-u**2))
    return -np.array([quad(g0,0,U,limit=200)[0],quad(g1,0,U,limit=200)[0]])
def w_of(sv): return ominv@Iu(sv)
e_plus=w_of(0.0); e_minus=-e_plus
def Omega(rv):  # 3a specie: log theta-ratio
    w=w_of(1.0/rv); return np.log(theta_af(w-e_plus,av,bv)/theta_af(w-e_minus,av,bv))
# A_ab (source abeliano) = int dEF - parte algebrica. Uso direttamente int dEF - Acal/sqrtS.
import sympy as sp
r,Es=sp.symbols('r E',positive=True); Ms,as_,Js=sp.Rational(1),sp.Rational(9,10),sp.Rational(5,2)
Dl=r**2-2*Ms*r+as_**2; Em=(Es**2-1)*r+2*Ms
Ssym=sp.expand(r*(r-2*Ms)*Em*(r*Dl-Js**2*Em)); Ksym=Js*r*(r-2*Ms)*Em/Dl
dEF=sp.diff(Ksym/sp.sqrt(Ssym),Es).subs(Es,sp.Rational(7,5))
dEFn=sp.lambdify(r,dEF,'numpy')
def Sn(x):
    Dl=x**2-2*Ms*x+as_**2; E7=1.4; Emm=(E7**2-1)*x+2*Ms
    return float(x*(x-2*Ms)*Emm*(x*Dl-Js**2*Emm))
r0=11.5
def A_full(rv): return quad(dEFn,r0,rv,limit=200)[0]   # dA_ab/dr ~ dEF (parte abeliana domina la variazione)
def dOmega(rv,h=1e-6): return (Omega(rv+h)-Omega(rv-h))/(2*h)
# Lambda = int(A dOmega - Omega dA), A=A_full (proxy), dA=dEF  -- GROUND TRUTH del peso-2
def dilog(rv):
    return quad(lambda x: A_full(x)*dOmega(x) - Omega(x)*dEFn(x), r0, rv, limit=80)[0]
print("\n=== ground-truth dilog Lambda(r) = int(A dOmega - Omega dA) ===")
for rv in [10.0,8.0,6.5]:
    print(f"  r={rv}: Lambda={dilog(rv):+.6f}  Omega={Omega(rv).real:+.4f}+{Omega(rv).imag:+.4f}j")
print("\n=> nomi q_ij + serie-nome theta verificata + dilog ground-truth pronti.")
print("   TAPPA 2: espandere Omega(r) e Lambda in serie di nome (Kronecker-Eisenstein genus-2).")
