# -*- coding: utf-8 -*-
# PHYSICS ANCHOR for the Vaidya v-branch off-shell closed form.
# closed wrap = -int(Gpr/Hpr) S_D dr  == off-shell sub-piece partB of the TRUE non-autonomous v flow;
# total partA+partB vs true flow slope ~2. Source Theta=m d_m (no dilation), m INCREASES (accretion).
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct, quad
from scipy.optimize import brentq
from scipy.interpolate import interp1d
r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
f=1-2*m/r; w=E**2-f
Hv=pr*(f-E**2)-1+sp.sqrt(w)*sp.sqrt(E**2*pr**2+J**2/r**2)
la=lambda ex: sp.lambdify((r,pr,m,E,J),ex,'numpy')
Hn=la(Hv); Hp=la(sp.diff(Hv,pr)); Hr=la(sp.diff(Hv,r)); Hm=la(sp.diff(Hv,m)); HJ=la(sp.diff(Hv,J))
G=sp.diff(Hv,J)/sp.diff(Hv,pr); Gpr=la(sp.diff(G,pr)); Gm=la(sp.diff(G,m))
E0,J0,r0,m0=1.4,2.5,12.0,1.0
def prof(rv,mv):
    ps=np.linspace(-40,-1e-4,4000); vals=[Hn(rv,p,mv,E0,J0) for p in ps]
    s=np.where(np.diff(np.sign(vals)))[0]; return brentq(lambda p:Hn(rv,p,mv,E0,J0),ps[s[0]],ps[s[0]+1])
ev=lambda lam,y:y[1]; ev.terminal=True; ev.direction=1
def flow(eps):
    def rhs(lam,y):
        rv,pv,ph=y; mv=m0*np.exp(eps*lam)      # accretion
        return [Hp(rv,pv,mv,E0,J0),-Hr(rv,pv,mv,E0,J0),HJ(rv,pv,mv,E0,J0)]
    so=solve_ivp(rhs,[0,300],[r0,prof(r0,m0),0.0],rtol=1e-11,atol=1e-13,max_step=0.01,dense_output=True,events=ev)
    lam=np.linspace(0,so.t[-1],12000); Y=so.sol(lam); return lam,Y[0],Y[1],Y[2]
lam0,rF,prF,phiF=flow(0.0)
ThetaH=m0*Hm(rF,prF,m0,E0,J0); S_D=ct(ThetaH,lam0,initial=0)
Hpr=Hp(rF,prF,m0,E0,J0); Gp=Gpr(rF,prF,m0,E0,J0)
partB=ct(Gp*S_D/Hpr,rF,initial=0)                       # Vaidya sign (m increases)
partA=ct(-Gp*(lam0*ThetaH)/Hpr + lam0*m0*Gm(rF,prF,m0,E0,J0),rF,initial=0)
dphi=partA+partB; dphi_r=interp1d(rF,dphi,fill_value='extrapolate',bounds_error=False)
phi0f=interp1d(rF,phiF,fill_value='extrapolate',bounds_error=False)
partB_r=interp1d(rF,partB,fill_value='extrapolate',bounds_error=False)
rc=np.linspace(8.0,11.0,1000)
print("(1) PHYSICS: total (partA+partB) vs TRUE non-autonomous v flow:")
epss=np.array([0.002,0.004,0.008,0.016]); res=[]
for eps in epss:
    _,rL,_,pL=flow(eps); pt=interp1d(rL,pL,fill_value='extrapolate',bounds_error=False)(rc)
    res.append(np.nanmax(np.abs(pt-(phi0f(rc)+eps*dphi_r(rc)))))
res=np.array(res); print(f"    slope={np.polyfit(np.log(epss),np.log(res),1)[0]:.2f}  (~2)")
# closed wrap (block1+block2) on same orbit == -partB
mv,Ev,Jv=1.0,1.4,2.5; DE=(E**2-1)*r+2*m; K1=sp.simplify(f-E**2); K2=sp.simplify(w*E**2-K1**2)
Q2v=(2*E**2*J**2*m*r-E**2*J**2*r**2+E**2*r**4+4*J**2*m**2-4*J**2*m*r+J**2*r**2); S_v=sp.expand(r*DE*Q2v)
Akv=sp.lambdify(r,(E**2*J*DE*r**4/Q2v).subs({m:mv,E:Ev,J:Jv}),'numpy'); Svn=sp.lambdify(r,S_v.subs({m:mv,E:Ev,J:Jv}),'numpy')
def pr0(R):
    g=lambda p:Hn(R,p,mv,Ev,Jv); ps=np.linspace(-40,-1e-4,4000); vs=[g(p) for p in ps]
    s=np.where(np.diff(np.sign(vs)))[0]; return brentq(g,ps[s[0]],ps[s[0]+1])
def inner(x): p=pr0(x); return mv*Hm(x,p,mv,Ev,Jv)/Hp(x,p,mv,Ev,Jv)
def Sig(x): return quad(inner,r0,x,limit=200)[0]
def wrap_closed(x): return -quad(lambda t:(Akv(t)/np.sqrt(max(Svn(t),0.0)))*Sig(t),r0,x,limit=120)[0]
print("\n(2) closed wrap == -partB (off-shell sub-piece of the true v flow):")
for x in [11.0,10.0,9.0]:
    print(f"    r={x}: closed={wrap_closed(x):+.7f}  (-partB)={-float(partB_r(x)):+.7f}  diff={abs(wrap_closed(x)+float(partB_r(x))):.1e}")
print("\n=> v-branch off-shell physics-anchored (closed==flow sub-piece; total slope ~2).")
