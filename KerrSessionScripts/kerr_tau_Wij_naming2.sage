# -*- coding: utf-8 -*-
# NAMING KLEINIANO v2: costante di Riemann ESATTA + base al branch point.
# Fix rispetto a v1: (1) Abel map basata al branch point (sigma-canonica di Baker),
# (2) ordine cicli coerente con period_matrix, (3) caratteristica sigma = RiemannConstantVector.
# Verifica: d/dr zeta_i(u) = dR_i/dr (BEL 2a specie) -> R_i CHIUSO in zeta_i tabulato.
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np
from scipy.integrate import quad

M,a,J,Ehat = 1.0,0.9,2.5,1.4
Rs = PolynomialRing(QQ, ['s','Y']); s,Y = Rs.gens()
lam=[QQ(1200),QQ(-2300),QQ(-11428),QQ(-5519),QQ(24700),QQ(62500)]
qpoly=sum(lam[i]*s**i for i in range(6))
X=SageRS(Y**2-qpoly, prec=100)
# --- periodi: uso period_matrix (ordine [A|B] verificato in kerr_psi_forward_abel) ---
PM=np.array(X.period_matrix(),dtype=complex)   # 2x4 = [A|B]
Amat,Bmat=PM[:,:2],PM[:,2:]
tau=np.array(X.riemann_matrix(),dtype=complex)
print("A^-1 B - tau =",np.linalg.norm(np.linalg.solve(Amat,Bmat)-tau))
# 2a specie eta: matrix_of_integral_values su [du1,du2,dr1,dr2]; VERIFICO che il blocco
# olomorfo coincida con period_matrix (stesso ordine cicli).
dr1=(lam[3]*s+2*lam[4]*s**2+3*lam[5]*s**3)/2; dr2=lam[5]*s**2/2
MIV=np.array(X.matrix_of_integral_values([Rs(1),s,dr1,dr2]),dtype=complex)  # 4x4 (diff x cicli)
print("blocco olo MIV[:2] vs [A|B]:",np.linalg.norm(MIV[:2]-PM))  # se ~0 stesso ordine
om,omp=MIV[:2,:2],MIV[:2,2:]; eta,etap=MIV[2:,:2],MIV[2:,2:]
ominv=np.linalg.inv(om); kappa=eta@ominv
print("kappa sym err=",np.linalg.norm(kappa-kappa.T))

def qn(sv): return sum(float(lam[i])*sv**i for i in range(6))
# branch point base: radice reale piu' piccola di q (turning in s = 1/r_min)
qroots=np.roots([float(lam[i]) for i in range(5,-1,-1)])
sb=float(min([z.real for z in qroots if abs(z.imag)<1e-8 and z.real>0]))
print("branch point base s_b=",sb," (r_min=",1/sb,")")
# Abel map UNnormalizzata basata a s_b: u(s)=int_{s_b}^{s}(1,t)/(2 sqrt q) dt
# regolarizzo la radice al branch point s=s_b - w^2
def u_unnorm(s_to):
    if abs(s_to-sb)<1e-14: return np.zeros(2,dtype=complex)
    W=np.sqrt(abs(sb-s_to))
    g0=lambda w: 2*w/(2*np.sqrt(abs(qn(sb-w**2))))
    g1=lambda w: 2*w*(sb-w**2)/(2*np.sqrt(abs(qn(sb-w**2))))
    I0=quad(g0,0,W,limit=200)[0]; I1=quad(g1,0,W,limit=200)[0]
    sgn=-1.0 if s_to<sb else 1.0
    return sgn*np.array([I0,I1],dtype=complex)   # int_{s_b}^{s_to}

# costante di Riemann (sigma-canonica). Provo abelfunctions RiemannConstantVector; se crasha,
# uso la caratteristica half-integer per hyperell. odd (Delta = somma Abel di g branch points).
delta=None
try:
    from abelfunctions import RiemannSurface as ARS, RiemannConstantVector
    fA=Y**2-qpoly
    Xa=ARS(fA)
    P0=Xa(sb, 0)     # branch point
    K=np.array(RiemannConstantVector(P0),dtype=complex)
    print("RiemannConstantVector K=",np.round(K,4))
    delta=('vec',K)
except Exception as ex:
    print("RiemannConstantVector fallito:",repr(ex)[:120])
    print("-> uso scan half-integer come fallback")

def logtheta_grad(zz, av, bv):
    z2=zz+tau@av+bv; th=complex(RiemannTheta(z2,tau))
    g=np.array([complex(RiemannTheta(z2,tau,derivs=[e])) for e in ([1,0],[0,1])])
    return 2j*np.pi*av + g/th
def zeta_klein_char(u, av, bv):
    zn=ominv@u; return -(kappa@u)+(ominv.T@logtheta_grad(zn,av,bv))
def zeta_klein_vecK(u, K):
    # theta con shift continuo K (RCV): log theta(ominv u + K); grad rispetto a u
    zn=ominv@u + K
    th=complex(RiemannTheta(zn,tau))
    g=np.array([complex(RiemannTheta(zn,tau,derivs=[e])) for e in ([1,0],[0,1])])
    return -(kappa@u)+(ominv.T@(g/th))

def dRi_dr(i,rv):
    sv=1.0/rv; Y_=np.sqrt(qn(sv)); dsdr=-1.0/rv**2
    Nd=[(lam[3]*sv+2*lam[4]*sv**2+3*lam[5]*sv**3)/4,(lam[5]*sv**2)/4]
    return Nd[i]/Y_*dsdr
def zeta_of_r(rv, spec):
    u=u_unnorm(1.0/rv)
    if spec[0]=='vec': return zeta_klein_vecK(u,spec[1])
    return zeta_klein_char(u,spec[1][0],spec[1][1])

import itertools
halfs=[np.array(v)/2 for v in itertools.product([0,1],repeat=2)]
chars=[('char',(av,bv)) for av in halfs for bv in halfs]
specs=([delta] if delta else [])+chars
rg2=np.linspace(11.0,6.5,7)
DR=np.array([[dRi_dr(0,rv),dRi_dr(1,rv)] for rv in rg2])
print("\n=== verifica d zeta_i/dr vs dR_i/dr (Baker: atteso identita') ===")
best=None
for spec in specs:
    try:
        DZ=np.array([ (zeta_of_r(rv+1e-5,spec)-zeta_of_r(rv-1e-5,spec))/(2e-5) for rv in rg2])
    except Exception: continue
    # confronto diretto DZ vs DR (naming zeta_i = +/- R_i): residuo relativo
    em=min(np.max(np.abs(DZ-DR)),np.max(np.abs(DZ+DR)))/max(np.max(np.abs(DR)),1e-30)
    lbl='RCV' if spec[0]=='vec' else f"char a={spec[1][0]} b={spec[1][1]}"
    if best is None or em<best[0]: best=(em,lbl,DZ)
    if spec[0]=='vec' or em<1e-2: print(f"  {lbl}: |dzeta -/+ dR|_rel = {em:.2e}")
print(f"\nMIGLIORE: {best[1]}  residuo diretto = {best[0]:.2e}")
print("  DZ[0] vs DR[0]:",np.round(best[2][0],4),"  vs",np.round(DR[0],4))
print("=> se ~0: R_i = zeta_i(u) DIRETTO (naming Baker, funzione tabulata). ")
