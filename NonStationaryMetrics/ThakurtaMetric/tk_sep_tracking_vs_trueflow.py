# -*- coding: utf-8 -*-
# DECISIVE physics test: does the tracked-source off-shell analytic dphi MATCH the true flow
# up to the separatrix turning point, or OVERSHOOT it (spurious power divergence => needs the
# moving-root dr_d/dlambda counterterm)?  Config E=1.2,a=0.9 (physical t-branch, Jc=2.9364, r_d=1.512).
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct
from scipy.optimize import brentq
from scipy.interpolate import interp1d
M,a,Ehat=1.0,0.9,1.2
rr,pr,Ess,Jss=sp.symbols('r pr E J_')
f2=1-2*M/rr; Dl2=rr**2-2*M*rr+a**2; b2=2*M*a/rr; v2=1-f2/Ess**2
P2=rr**2+a**2+2*M*a**2/rr; Pb2=P2+b2**2/Ess**2
H2=Jss*b2*v2/Pb2+sp.sqrt(Dl2*v2/Pb2)*sp.sqrt((Dl2/rr**2)*pr**2+Jss**2/Pb2)-1
H2n=sp.lambdify((rr,pr,Ess,Jss),H2,'numpy')
Hp=sp.diff(H2,pr); HJ=sp.diff(H2,Jss); HE=sp.diff(H2,Ess); G=HJ/Hp
dHp=sp.lambdify((rr,pr,Ess,Jss),Hp,'numpy'); dHr=sp.lambdify((rr,pr,Ess,Jss),sp.diff(H2,rr),'numpy')
dHJ=sp.lambdify((rr,pr,Ess,Jss),HJ,'numpy'); dHE=sp.lambdify((rr,pr,Ess,Jss),HE,'numpy')
dG_pr=sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,pr),'numpy')
dG_E =sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,Ess),'numpy')
dG_J =sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,Jss),'numpy')

# physical separatrix at Ehat=1.2, a=0.9 (from tk_sep_physical_separatrix.py)
Jc=2.93635; r_d=1.51229
print(f"E=1.2 a=0.9: Jc={Jc:.5f}  r_d={r_d:.5f}")

def prof(rv,Ev,Jvv,want_ingoing=True):
    pg=np.linspace(-120,120,4001); Hv=H2n(rv,pg,Ev,Jvv)
    rts=[brentq(lambda p:H2n(rv,p,Ev,Jvv),pg[i],pg[i+1]) for i in range(len(pg)-1)
         if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1]<0]
    ing=[p for p in rts if dHp(rv,p,Ev,Jvv)<0]
    return (min(ing) if ing else np.nan)

ev=lambda lam,y:y[1]; ev.terminal=True; ev.direction=1
def flow(eps,J0,r0):
    def rhs(lam,y):
        rv,pv,ph=y; s=np.exp(-eps*lam); Ee=Ehat*s; Je=J0*s
        return [dHp(rv,pv,Ee,Je),-dHr(rv,pv,Ee,Je)-eps*pv,dHJ(rv,pv,Ee,Je)]
    p0=prof(r0,Ehat,J0)
    so=solve_ivp(rhs,[0,2000],[r0,p0,0.0],rtol=1e-11,atol=1e-13,max_step=0.01,dense_output=True,events=ev)
    lam=np.linspace(0,so.t[-1],20000); Y=so.sol(lam); return lam,Y[0],Y[1],Y[2]

r0=8.0
# evaluate overshoot on the CLEAN arc ABOVE r_d (avoid horizon plunge / pole crossing)
r_ev=np.array([3.0,2.5,2.0])
print(f"r_d={r_d}; eval on clean arc above r_d at r={list(r_ev)}")
print(f"{'J0':>7} {'%Jc':>6} {'r_min':>7} {'ovsh@3.0':>11} {'ovsh@2.5':>11} {'ovsh@2.0':>11}  (overshoot=|analytic-(true-froz)/eps|)")
for J0 in [2.60,2.75,2.85,2.90,2.92]:
    lam0,rF,prF,phiF=flow(0.0,J0,r0)
    EulerH=Ehat*dHE(rF,prF,Ehat,J0)+J0*dHJ(rF,prF,Ehat,J0)
    Hpr_fr=dHp(rF,prF,Ehat,J0)
    S_lam=ct(EulerH+prF*Hpr_fr,lam0,initial=0)
    EulerG=Ehat*dG_E(rF,prF,Ehat,J0)+J0*dG_J(rF,prF,Ehat,J0)
    Gpr=dG_pr(rF,prF,Ehat,J0)
    partA=Gpr*(lam0*EulerH)/Hpr_fr-lam0*EulerG
    partB=-Gpr*S_lam/Hpr_fr
    dphi=ct(partA+partB,rF,initial=0)
    phi0f=interp1d(rF,phiF,fill_value='extrapolate',bounds_error=False)
    dphi_r=interp1d(rF,dphi,fill_value='extrapolate',bounds_error=False)
    eps=0.001
    _,rL,_,pL=flow(eps,J0,r0)
    ptrue=interp1d(rL,pL,fill_value='extrapolate',bounds_error=False)
    ov=[abs((ptrue(rc)-phi0f(rc))/eps-dphi_r(rc)) for rc in r_ev]
    print(f"{J0:7.3f} {100*J0/Jc:6.1f} {rF[-1]:7.4f} {ov[0]:11.4e} {ov[1]:11.4e} {ov[2]:11.4e}")
print("\nOvershoot GROWS as J0->Jc on the CLEAN arc above r_d => analytic (tracked source) has a")
print("spurious separatrix divergence the true flow lacks => OPEN (moving-root dr_d counterterm needed).")
print("Flat/small => tracked source already matches => CLOSED.")
