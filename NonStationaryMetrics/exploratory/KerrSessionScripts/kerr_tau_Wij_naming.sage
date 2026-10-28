# -*- coding: utf-8 -*-
# KLEINIAN NAMING (TK tau, genus-2): the two canonical abelian integrals of the
# W_ij closure ARE Baker's tabulated functions.
#   R_i = int dr_i (BEL second kind)    -> Kleinian zeta_i(u) (differential)
#   L   = int ds/(s Y) (third kind)      -> log[sigma(u-e+)/sigma(u-e-)]
# Check at the DIFFERENTIAL level (robust): d/dr(Kleinian function) = canonical differential.
# Odd model (quintic), params M=1,a=9/10,E=7/5,J=5/2 (same as Bel/thirdkind).
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
# first-kind periods (a,b) + second kind (BEL) -> omega, tau, kappa
dr1=(lam[3]*s+2*lam[4]*s**2+3*lam[5]*s**3)/2; dr2=lam[5]*s**2/2
Pi=np.array(X.matrix_of_integral_values([Rs(1),s,dr1,dr2]),dtype=complex)
omega,omega_p = Pi[:2,:2],Pi[:2,2:]  # WRONG slice? Pi is 4x4: rows=diff, cols=cycles
# matrix_of_integral_values: rows=differentials [du1,du2,dr1,dr2], columns=4 cycles (a1,a2,b1,b2)
om=Pi[:2,:2]; omp=Pi[:2,2:]; eta=Pi[2:,:2]; etap=Pi[2:,2:]
tau=np.linalg.solve(om,omp); ominv=np.linalg.inv(om); kappa=eta@ominv
print("tau sym err=",np.linalg.norm(tau-tau.T)," kappa sym err=",np.linalg.norm(kappa-kappa.T))
Leg=omp@eta.T-om@etap.T; print("Legendre/(2pi i)=\n",np.round(Leg/(2j*np.pi),3)," (canonical -I/2)")

def qn(sv): return sum(float(lam[i])*sv**i for i in range(6))
r0=12.0; s0=1.0/r0
def Sn(rv):
    Dl=rv**2-2*M*rv+a**2; Emu=(Ehat**2-1)*rv+2*M
    return rv*(rv-2*M)*Emu*(rv*Dl-J**2*Emu)
# UNnormalised Abel map u(r)=int_{s0}^{s}(1,t)/(2 sqrt q) dt (measure /2 = Legendre convention -I/2)
def u_unnorm(sv):
    I0=quad(lambda t:1.0/(2*np.sqrt(qn(t))),s0,sv,limit=200)[0]
    I1=quad(lambda t:t/(2*np.sqrt(qn(t))),s0,sv,limit=200)[0]
    return np.array([I0,I1],dtype=complex)
def zval(sv): return u_unnorm(sv)                      # u (Kleinian arg, unnormalised)

# zeta_i(u) = -(kappa u)_i + (omega^{-T} grad_z log theta[delta](omega^{-1} u))_i
halfs=[np.array(v)/2 for v in itertools.product([0,1],repeat=2)]
chars=[(av,bv) for av in halfs for bv in halfs]
def logtheta_grad(zz, av, bv):   # grad_z log theta[av,bv](zz)
    z2=zz+tau@av+bv
    th=complex(RiemannTheta(z2,tau))
    g=np.array([complex(RiemannTheta(z2,tau,derivs=[e])) for e in ([1,0],[0,1])])
    return 2j*np.pi*av + g/th
def zeta_klein(u, av, bv):
    zn=ominv@u
    return -(kappa@u) + (ominv.T @ logtheta_grad(zn,av,bv))

# --- canonical differential dR_i/dr (BEL second kind) and du_i/dr ---
def dRi_dr(i,rv):   # dr_i in r: Ndr_i(s)/Y * ds/dr,  ds/dr=-1/r^2
    sv=1.0/rv; Y_=np.sqrt(qn(sv)); dsdr=-1.0/rv**2
    Nd=[(lam[3]*sv+2*lam[4]*sv**2+3*lam[5]*sv**3)/4,(lam[5]*sv**2)/4]  # BEL numerators /4
    return Nd[i]/Y_*dsdr
def duidr(i,rv):
    sv=1.0/rv; Y_=np.sqrt(qn(sv)); dsdr=-1.0/rv**2
    return (1.0 if i==0 else sv)/(2*np.sqrt(qn(sv)))*dsdr   # /2 convention
# dzeta_i/dr = grad_u zeta_i . du/dr ; numerical through a finite difference in r on zeta_klein
print("\n=== NAMING R_i = zeta_i(u): scan of 16 characteristics, lstsq dzeta_i = A dR + C du ===")
rg2=np.linspace(11.0,6.0,7)
DRr=np.array([[dRi_dr(0,rv),dRi_dr(1,rv)] for rv in rg2])
DUr=np.array([[duidr(0,rv),duidr(1,rv)] for rv in rg2])
Xdes=np.hstack([DRr,DUr])
best=None
for ci,(av,bv) in enumerate(chars):
    DZ=[]
    ok=True
    for rv in rg2:
        h=1e-5
        try:
            dz=(zeta_klein(zval(1.0/(rv+h)),av,bv)-zeta_klein(zval(1.0/(rv-h)),av,bv))/(2*h)
        except Exception: ok=False; break
        DZ.append(dz)
    if not ok: continue
    DZ=np.array(DZ)
    tot=0.0; coefs=[]
    for i in range(2):
        coef,_,_,_=np.linalg.lstsq(Xdes,DZ[:,i],rcond=None)
        rr=float(np.max(np.abs(Xdes@coef-DZ[:,i]))/max(np.max(np.abs(DZ[:,i])),1e-30))
        tot=max(tot,rr); coefs.append(coef)
    cleanness=float(sum(abs(c[2])+abs(c[3]) for c in coefs))  # small |du coeff| = clean
    score=cleanness if tot<1e-4 else 1e6+tot                  # low residual first, then cleanness
    if best is None or score<best[0]:
        best=(score,ci,av,bv,coefs,tot)
_,ci,av,bv,coefs,tot=best
print(f"  BEST char #{ci} a={av} b={bv}  (max lstsq residual={tot:.2e})")
for i in range(2):
    c=coefs[i]
    print(f"   zeta_{i+1}: {c[0]:+.4f} dR1 {c[1]:+.4f} dR2 {c[2]:+.4f} du1 {c[3]:+.4f} du2")
print("=> Baker: expected dzeta_i = dR_i (identity on dR, du small/integer). If residual~0 and the coeffs are")
print("   simple: R_i CLOSED in zeta_i (tabulated function). If the coeffs are messy: the convention is still off.")
