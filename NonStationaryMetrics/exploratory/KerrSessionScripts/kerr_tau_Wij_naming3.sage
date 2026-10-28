# -*- coding: utf-8 -*-
# KLEINIAN NAMING v3: decisive test. R_i = zeta_i(u(r)+C), C = base-point offset
# (Baker bases at the branch point s=inf; my u is based elsewhere -> a constant C on the Jacobian).
# zeta is not translation invariant: if there EXISTS a unique C (2 complex numbers) with
#   d/dr zeta_i(u(r)+C) = dR_i/dr   for i=1,2 and all r  => R_i CLOSED in zeta_i.
# Find C by minimisation, then CHECK the residual on held-out points.
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares

M,a,J,Ehat = 1.0,0.9,2.5,1.4
Rs = PolynomialRing(QQ, ['s','Y']); s,Y = Rs.gens()
lam=[QQ(1200),QQ(-2300),QQ(-11428),QQ(-5519),QQ(24700),QQ(62500)]
qpoly=sum(lam[i]*s**i for i in range(6))
X=SageRS(Y**2-qpoly, prec=100)
PM=np.array(X.period_matrix(),dtype=complex); tau=np.array(X.riemann_matrix(),dtype=complex)
dr1=(lam[3]*s+2*lam[4]*s**2+3*lam[5]*s**3)/2; dr2=lam[5]*s**2/2
MIV=np.array(X.matrix_of_integral_values([Rs(1),s,dr1,dr2]),dtype=complex)
om,eta=MIV[:2,:2],MIV[2:,:2]; ominv=np.linalg.inv(om); kappa=eta@ominv
def qn(sv): return sum(float(lam[i])*sv**i for i in range(6))
r0=12.0; s0=1.0/r0
def u_unnorm(s_to):   # based at s0 (then corrected with C)
    I0=quad(lambda t:1.0/(2*np.sqrt(qn(t))),s0,s_to,limit=200)[0]
    I1=quad(lambda t:t/(2*np.sqrt(qn(t))),s0,s_to,limit=200)[0]
    return np.array([I0,I1],dtype=complex)
def zeta_klein(v):   # v = u + C (full arg); zeta_i=-(kappa (v-?))... use the pure theta form:
    # zeta_i(v) = -(kappa v)_i + (ominv^T grad_z log theta(ominv v))_i   [characteristic 0]
    zn=ominv@v; th=complex(RiemannTheta(zn,tau))
    g=np.array([complex(RiemannTheta(zn,tau,derivs=[e])) for e in ([1,0],[0,1])])
    return -(kappa@v)+(ominv.T@(g/th))
def dRi_dr(i,rv):
    sv=1.0/rv; Y_=np.sqrt(qn(sv)); dsdr=-1.0/rv**2
    Nd=[(lam[3]*sv+2*lam[4]*sv**2+3*lam[5]*sv**3)/4,(lam[5]*sv**2)/4]
    return Nd[i]/Y_*dsdr
def dzeta_dr(rv,C):
    h=1e-6
    vp=u_unnorm(1.0/(rv+h))+C; vm=u_unnorm(1.0/(rv-h))+C
    return (zeta_klein(vp)-zeta_klein(vm))/(2*h)

rg_fit=np.linspace(11.0,7.0,5); rg_test=np.array([10.5,8.5,6.5])
DRfit=np.array([[dRi_dr(0,rv),dRi_dr(1,rv)] for rv in rg_fit])
def resid(Cr):
    C=np.array([Cr[0]+1j*Cr[1],Cr[2]+1j*Cr[3]])
    out=[]
    for k,rv in enumerate(rg_fit):
        dz=dzeta_dr(rv,C)
        for i in range(2):
            d=dz[i]-DRfit[k,i]; out+=[d.real,d.imag]
    return out
print("searching for C (Jacobian base offset) such that dzeta_i(u+C)/dr = dR_i/dr ...")
best=None
for C0 in [ [0,0,0,0],[0.1,0.1,0.1,0.1],[0.5,0,0.5,0],[0,0.5,0,0.5],
            [0.25,0.25,-0.25,0.25],[-0.3,0.2,0.3,-0.2] ]:
    sol=least_squares(resid,C0,method='lm',max_nfev=4000)
    if best is None or sol.cost<best.cost: best=sol
C=np.array([best.x[0]+1j*best.x[1],best.x[2]+1j*best.x[3]])
print("C =",np.round(C,5),"  fit cost =",best.cost)
# CHECK on held-out points
print("\n=== held-out CHECK: dzeta_i(u+C)/dr vs dR_i/dr ===")
mx=0
for rv in rg_test:
    dz=dzeta_dr(rv,C); dR=np.array([dRi_dr(0,rv),dRi_dr(1,rv)])
    e=np.max(np.abs(dz-dR)); mx=max(mx,e)
    print(f"  r={rv}: dzeta={np.round(dz,5)}  dR={np.round(dR,5)}  |diff|={e:.2e}")
print(f"\nMAX held-out diff = {mx:.2e}")
print("=> if ~1e-6: R_i = zeta_i(u+C) CONFIRMED (Baker naming, C=base offset). Naming CLOSED.")
print("   (if not: the second kind is NOT a bare zeta_i -> the difference form at e_pm is needed).")
