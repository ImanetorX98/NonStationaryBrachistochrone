# -*- coding: utf-8 -*-
# NAMING KLEINIANO forma-DIFFERENZA (TK tau, genus-2). Riuso la convenzione VALIDATA
# di kerr_thirdkind_theta_closed (delta dispari, e_pm=+-w(r=inf), base al branch point,
# misura 1/sqrt(q)). Nomino i differenziali canonici in r ai punti e_pm:
#   3a specie (U_2): D3 = [zeta_d(w-e+) - zeta_d(w-e-)] . dw/dr   (log theta-ratio)
#   2a specie (U_3,U_4): G_pm,k = [grad zeta_d(w-e_pm)]_k . dw/dr  (zeta shiftata)
# VERIFICA (differenziale, robusto al divisore-theta perche' argomenti SHIFTATI a e_pm):
#   r^k/sqrt(S) = sum coeff * {du1/dr,du2/dr, D3, G+_1,G+_2,G-_1,G-_2}  (residuo ~0).
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
def qn(sv): return sum(float(lam[i])*sv**i for i in range(6))
rmin=4.046197656444178; s_b=1.0/rmin
def Iu(s_to):
    U=np.sqrt(s_b-s_to)
    g0=lambda u:2*u*1.0/np.sqrt(qn(s_b-u**2)); g1=lambda u:2*u*(s_b-u**2)/np.sqrt(qn(s_b-u**2))
    return -np.array([quad(g0,0,U,limit=200)[0],quad(g1,0,U,limit=200)[0]])
def w_of(s_to): return ominv@Iu(s_to)
e_plus=w_of(0.0); e_minus=-e_plus

# caratteristiche odd
def theta_d(z,av,bv,order=0):
    av=np.array(av); bv=np.array(bv); zz=z+tau@av+bv
    pref=np.exp(1j*np.pi*(av@tau@av)+2j*np.pi*(av@(z+bv)))
    th=complex(RiemannTheta(zz,tau))
    if order==0: return pref*th
    g=np.array([complex(RiemannTheta(zz,tau,derivs=[e])) for e in ([1,0],[0,1])])
    if order==1: return 2j*np.pi*av + g/th          # grad log theta[delta]
    H=np.array([[complex(RiemannTheta(zz,tau,derivs=[p,q])) for q in ([1,0],[0,1])] for p in ([1,0],[0,1])])
    return H/th - np.outer(g/th,g/th)                # Hess log theta[delta]
halfs=[np.array(v)/2 for v in itertools.product([0,1],repeat=2)]
odd=[(av,bv) for av in halfs for bv in halfs if abs(theta_d(np.zeros(2),av,bv))<1e-6]
print("odd chars:",len(odd))

def Sn(rv):
    Dl=rv**2-2*M*rv+a**2; Emu=(Ehat**2-1)*rv+2*M
    return rv*(rv-2*M)*Emu*(rv*Dl-J**2*Emu)
rg=np.linspace(11.5,rmin+0.4,22)
sr=1.0/rg; dsdr=-1.0/rg**2
du1dr=(1/np.sqrt(qn(sr)))*dsdr; du2dr=(sr/np.sqrt(qn(sr)))*dsdr
W=np.array([w_of(sv) for sv in sr]); dwdr=(ominv@np.vstack([du1dr,du2dr])).T
sqrtS=np.sqrt(Sn(rg))

print("\n delta  | k=2 resid | k=3 resid | k=4 resid  (r^k/sqrtS su base Kleiniana differenza)")
for di,(av,bv) in enumerate(odd):
    # colonne base ai punti r
    col_du1=du1dr; col_du2=du2dr
    D3=np.array([ (theta_d(W[i]-e_plus,av,bv,1)-theta_d(W[i]-e_minus,av,bv,1))@dwdr[i] for i in range(len(rg))])
    Gp=np.array([ theta_d(W[i]-e_plus,av,bv,2)@dwdr[i] for i in range(len(rg))])   # (N,2)
    Gm=np.array([ theta_d(W[i]-e_minus,av,bv,2)@dwdr[i] for i in range(len(rg))])
    B=np.column_stack([col_du1,col_du2,D3,Gp[:,0],Gp[:,1],Gm[:,0],Gm[:,1]])
    resids=[]
    for k in [2,3,4]:
        target=rg**k/sqrtS
        coef,_,_,_=np.linalg.lstsq(B,target,rcond=None)
        resids.append(np.max(np.abs(B@coef-target))/max(np.max(np.abs(target)),1e-30))
    if min(resids)<1e-2 or di<2 or max(resids)<0.5:
        print(f"  #{di} a={av} b={bv} | {resids[0]:.2e} | {resids[1]:.2e} | {resids[2]:.2e}")
print("\n=> se una delta da' residui ~0 per k=2,3,4: U_2,U_3,U_4 CHIUSI in forma-differenza")
print("   Kleiniana (theta[delta] dispari agli e_pm). Naming peso-1 robusto al divisore-theta.")
