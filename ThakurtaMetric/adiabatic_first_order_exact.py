# -*- coding: utf-8 -*-
# EXACT first-order adiabatic correction (closes referee 4.6 to O(eps^2), slope ~2).
# Canonical perturbation theory of the non-autonomous optical Hamiltonian flow H2,
# with E_eff=Ehat/A, J_eff=J/A, A=exp(eps*lambda):
#   delta p_r/eps = [lambda*Euler_H - S]/H_pr,   Euler_H=E H_E + J H_J,  S=int_0^lambda Euler_H dlambda
#   delta phi/eps = int[ G_pr*(delta p_r/eps) - lambda*(E G_E + J G_J) ] dr,   G=H_J/H_pr=dphi/dr
# Verified: coeff -0.3298 = true -0.328 ; robust slope ~1.90 (=> O(eps^2), first order CLOSED).
# Decomposition: partA (single integral, ~-full Euler) + partB (OFF-SHELL, S nested = length-2
#   iterated integral, ~+1/2 full Euler). The papers 1/2-Euler is the approximation partB~=-1/2 partA
#   (accurate to ~2%); this file is the exact term. Closed form = partA (clock x algebraic) +
#   partB (length-2 iterated integral, genus-2 dilogarithm structure).
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct
from scipy.optimize import brentq
from scipy.interpolate import interp1d
M,a,Ehat=1.0,0.9,1.4; J0,r0=6.0,12.0
rr,pr,Ess,Jss=sp.symbols('r pr E J_')
f2=1-2*M/rr; Dl2=rr**2-2*M*rr+a**2; b2=2*M*a/rr; v2=1-f2/Ess**2
P2=rr**2+a**2+2*M*a**2/rr; Pb2=P2+b2**2/Ess**2
H2=Jss*b2*v2/Pb2+sp.sqrt(Dl2*v2/Pb2)*sp.sqrt((Dl2/rr**2)*pr**2+Jss**2/Pb2)-1
H2n=sp.lambdify((rr,pr,Ess,Jss),H2,'numpy')
Hp=sp.diff(H2,pr); HJ=sp.diff(H2,Jss); HE=sp.diff(H2,Ess)
dHp=sp.lambdify((rr,pr,Ess,Jss),Hp,'numpy'); dHr=sp.lambdify((rr,pr,Ess,Jss),sp.diff(H2,rr),'numpy')
dHJ=sp.lambdify((rr,pr,Ess,Jss),HJ,'numpy'); dHE=sp.lambdify((rr,pr,Ess,Jss),HE,'numpy')
def prof(rv,E,Jv):
    pg=np.linspace(-80,80,3001); Hv=H2n(rv,pg,E,Jv)
    rts=[brentq(lambda p:H2n(rv,p,E,Jv),pg[i],pg[i+1]) for i in range(len(pg)-1)
         if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1]<0]
    ing=[p for p in rts if dHp(rv,p,E,Jv)<0]; return min(ing) if ing else np.nan
ev=lambda lam,y:y[1]; ev.terminal=True; ev.direction=1
def flow(eps):
    def rhs(lam,y):
        rv,pv,ph=y; s=np.exp(-eps*lam); E=Ehat*s; Jv=J0*s
        return [dHp(rv,pv,E,Jv),-dHr(rv,pv,E,Jv),dHJ(rv,pv,E,Jv)]
    so=solve_ivp(rhs,[0,300],[r0,prof(r0,Ehat,J0),0.0],rtol=1e-12,atol=1e-14,max_step=0.004,dense_output=True,events=ev)
    lam=np.linspace(0,so.t[-1],16000); Y=so.sol(lam); return lam,Y[0],Y[1],Y[2]
lam0,rF,prF,phiF=flow(0.0)
# Euler_H(lam) and S(lam) along frozen orbit
EulerH=Ehat*dHE(rF,prF,Ehat,J0)+J0*dHJ(rF,prF,Ehat,J0)
S_lam=ct(EulerH,lam0,initial=0)
# analytic dp_r/eps(r) = [lam*EulerH - S]/H_pr
Hpr_fr=dHp(rF,prF,Ehat,J0)
dpr_an=(lam0*EulerH - S_lam)/Hpr_fr
dpr_an_r=interp1d(rF,dpr_an,bounds_error=False,fill_value='extrapolate')
pr0_r=interp1d(rF,prF,bounds_error=False,fill_value='extrapolate')
# numeric p_r^1(r)
eps=0.0005; _,rL,prL,_=flow(eps)
pr1_num=(interp1d(rL,prL,bounds_error=False,fill_value='extrapolate')(rF)-prF)/eps
rc=np.linspace(8.0,11.0,1500)
an=dpr_an_r(rc); nu=interp1d(rF,pr1_num,bounds_error=False,fill_value='extrapolate')(rc)
print("STEP1: analytic dp_r/eps vs numeric p_r^1 on [8,11]:")
print(f"  analytic mean={np.mean(an):+.4f}  numeric mean={np.mean(nu):+.4f}  max|diff|={np.nanmax(np.abs(an-nu)):.2e}")
print(f"  match? {np.nanmax(np.abs(an-nu))<0.02*max(abs(np.mean(an)),0.1)}")

# STEP 2: G=H_J/H_pr ; G_pr, G_E, G_J ; delta phi/eps = int[ G_pr*dpr - lam*(E G_E + J G_J) ]dr
G=HJ/Hp
dG_pr=sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,pr),'numpy')
dG_E =sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,Ess),'numpy')
dG_J =sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,Jss),'numpy')
phi0f=interp1d(rF,phiF,bounds_error=False,fill_value='extrapolate')
lam_r=interp1d(rF,lam0,bounds_error=False,fill_value='extrapolate')
# integrate along frozen orbit in r (rF decreasing). integrand as fn of r:
integ=dG_pr(rF,prF,Ehat,J0)*dpr_an - lam0*(Ehat*dG_E(rF,prF,Ehat,J0)+J0*dG_J(rF,prF,Ehat,J0))
# cumulative in r from r0: rF is decreasing, ct handles via x=rF
dphi_eps=ct(integ,rF,initial=0)
dphi_eps_r=interp1d(rF,dphi_eps,bounds_error=False,fill_value='extrapolate')
s_=-1.0; phi0_shape=lambda x:-phi0f(x)
# note: dphi from flow phi is +, shape is - ; the correction integ built in flow-phi convention
print("\nSTEP2: delta phi/eps coeff on [8,11]: %.4f  (true correction ~ -0.328 in flow-phi conv, +0.328 shape)"%np.mean(dphi_eps_r(rc)))
epss=np.array([0.0025,0.005,0.01,0.02])
res=[]
for eps in epss:
    _,rL,_,pL=flow(eps); pt=interp1d(rL,pL,bounds_error=False,fill_value='extrapolate')(rc)
    # flow phi convention directly: phi_true(flow) vs phi0(flow)+eps*dphi_eps
    pred=phi0f(rc)+eps*dphi_eps_r(rc)
    res.append(np.nanmax(np.abs(pt-pred)))
res=np.array(res); print(f"  SLOPE(full PT) = {np.polyfit(np.log(epss),np.log(res),1)[0]:.2f}  res@0.0025={res[0]:.2e}")
print(f"  (compare: 1/2 Euler-only gave slope 1.25, res 8.4e-5)")

# STEP3: robust slope + decompose into (A) single-integral part and (B) off-shell S-part
# dphi/eps = int[ G_pr*(lam*EulerH)/Hpr - lam*EulerG ]dr  +  int[ -G_pr*S/Hpr ]dr
#            \_______ part A (single integral) _______/     \__ part B: off-shell (S = nested) __/
EulerG=Ehat*dG_E(rF,prF,Ehat,J0)+J0*dG_J(rF,prF,Ehat,J0)
Gpr=dG_pr(rF,prF,Ehat,J0)
partA_int = Gpr*(lam0*EulerH)/Hpr_fr - lam0*EulerG
partB_int = -Gpr*S_lam/Hpr_fr
A=ct(partA_int,rF,initial=0); B=ct(partB_int,rF,initial=0)
Ar=interp1d(rF,A,bounds_error=False,fill_value='extrapolate'); Br=interp1d(rF,B,bounds_error=False,fill_value='extrapolate')
print("STEP3 decomposition on [8,11]: partA=%.4f  partB(off-shell S)=%.4f  sum=%.4f"%(np.mean(Ar(rc)),np.mean(Br(rc)),np.mean(Ar(rc)+Br(rc))))
print(f"  => off-shell S-part is {100*abs(np.mean(Br(rc)))/abs(np.mean(Ar(rc)+Br(rc))):.0f}% of total")
# robust slope with 6 eps
epss=np.array([0.002,0.004,0.008,0.016,0.032]); res=[]
for eps in epss:
    _,rL,_,pL=flow(eps); pt=interp1d(rL,pL,bounds_error=False,fill_value='extrapolate')(rc)
    res.append(np.nanmax(np.abs(pt-(phi0f(rc)+eps*(Ar(rc)+Br(rc))))))
res=np.array(res); print(f"  robust SLOPE = {np.polyfit(np.log(epss),np.log(res),1)[0]:.2f}")
