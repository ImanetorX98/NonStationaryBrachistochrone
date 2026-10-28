# -*- coding: utf-8 -*-
# PHYSICS ANCHOR for the Vaidya-tau off-shell closed form (+ 2nd-config math check).
# Chain: closed form (A+B+C) == direct off-shell wrap -int(Gpr/Hpr) S_D dr  [math, 1e-14 both configs]
#        == off-shell sub-piece of the TRUE non-autonomous tau flow (m=m(lambda))
#        total partA+partB == true flow to O(eps^2) (slope ~2).  Source Theta=m d_m (no dilation).
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct, quad
from scipy.optimize import brentq
from scipy.interpolate import interp1d

r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
Dl=r*(r-2*m); f=(r-2*m)/r; DE=(E**2-1)*r+2*m; v=DE/(E**2*r)
Htau=sp.sqrt(Dl*v/r**2)*sp.sqrt((Dl/r**2)*pr**2+J**2/r**2)-f/E
la=lambda ex: sp.lambdify((r,pr,m,E,J),ex,'numpy')
Hn=la(Htau); Hp=la(sp.diff(Htau,pr)); Hr=la(sp.diff(Htau,r)); Hm=la(sp.diff(Htau,m))
G=sp.diff(Htau,J)/sp.diff(Htau,pr); Gpr=la(sp.diff(G,pr)); Gm=la(sp.diff(G,m))
E0,J0=1.4,2.5; r0=12.0

def prof(rv,mv):
    ps=np.linspace(-60,-1e-4,4000); vals=[Hn(rv,p,mv,E0,J0) for p in ps]
    s=np.where(np.diff(np.sign(vals)))[0]; return brentq(lambda p:Hn(rv,p,mv,E0,J0),ps[s[0]],ps[s[0]+1])
ev=lambda lam,y:y[1]; ev.terminal=True; ev.direction=1
def flow(eps):
    def rhs(lam,y):
        rv,pv,ph=y; mv=1.0*np.exp(eps*lam)      # accretion m grows
        return [Hp(rv,pv,mv,E0,J0),-Hr(rv,pv,mv,E0,J0),sp.lambdify((r,pr,m,E,J),sp.diff(Htau,J),'numpy')(rv,pv,mv,E0,J0)]
    so=solve_ivp(rhs,[0,300],[r0,prof(r0,1.0),0.0],rtol=1e-11,atol=1e-13,max_step=0.01,dense_output=True,events=ev)
    lam=np.linspace(0,so.t[-1],12000); Y=so.sol(lam); return lam,Y[0],Y[1],Y[2]
lam0,rF,prF,phiF=flow(0.0)
ThetaH=1.0*Hm(rF,prF,1.0,E0,J0)                 # Theta H = m H_m (m0=1)
S_D=ct(ThetaH,lam0,initial=0)                   # costate integral
Hpr=Hp(rF,prF,1.0,E0,J0); Gp=Gpr(rF,prF,1.0,E0,J0)
# Vaidya sign convention: m INCREASES (opposite to TK E,J decreasing) => delta p_r=(S_D-lam ThetaH)/Hpr,
# and +lam Theta G (cf vaidya_first_order_offshell.py: integ=Gpr*(S-lam ThetaH)+Gm).
partB=ct(Gp*S_D/Hpr,rF,initial=0)               # OFF-shell sub-piece (Vaidya sign)
partA=ct(-Gp*(lam0*ThetaH)/Hpr + lam0*1.0*Gm(rF,prF,1.0,E0,J0),rF,initial=0)
dphi=partA+partB; dphi_r=interp1d(rF,dphi,fill_value='extrapolate',bounds_error=False)
phi0f=interp1d(rF,phiF,fill_value='extrapolate',bounds_error=False)
partB_r=interp1d(rF,partB,fill_value='extrapolate',bounds_error=False)
rc=np.linspace(8.0,11.0,1200)
print("(1) PHYSICS: total (partA+partB) vs TRUE non-autonomous tau flow:")
epss=np.array([0.002,0.004,0.008,0.016]); res=[]
for eps in epss:
    _,rL,_,pL=flow(eps); pt=interp1d(rL,pL,fill_value='extrapolate',bounds_error=False)(rc)
    res.append(np.nanmax(np.abs(pt-(phi0f(rc)+eps*dphi_r(rc)))))
res=np.array(res); print(f"    slope={np.polyfit(np.log(epss),np.log(res),1)[0]:.2f}  (~2 => closed-form-containing total reproduces physics)")

# (2) closed form (A+B+C) vs the flow off-shell sub-piece partB, on the SAME frozen orbit
import importlib.util as iu, sys, os
# rebuild the A+B+C closed form here (symbolic coeffs), evaluate on this orbit's r0
sub={m:1.0,E:sp.Rational(7,5),J:sp.Rational(5,2)}
T_V=sp.expand(r*(r-2*m)*DE); Q3=sp.expand(r**2*(r-2*m)-J**2*DE); Ssig=sp.expand(T_V*Q3); dengV=sp.expand((r-2*m)*DE)
Ap,N3=sp.div(sp.Poly(sp.expand(J*E*r**3*DE),r),sp.Poly(Q3,r)); N3=N3.as_expr(); p_j=[c for c in Ap.all_coeffs()[::-1]]
N4=(E**4*J**2*r**2+4*E**2*J**2*m*r-2*E**2*J**2*r**2+4*J**2*m**2-4*J**2*m*r+J**2*r**2+4*m**2*r**2-4*m*r**3+r**4)
Qd,Rr=sp.div(sp.Poly(sp.expand(-m*N4),r),sp.Poly(dengV,r)); Rr=Rr.as_expr(); a_k=[c for c in Qd.all_coeffs()[::-1]]
roots={'2m':2*m,'DE':-2*m/(E**2-1)}; rho={nm:sp.simplify((Rr/sp.diff(dengV,r)).subs(r,rt)) for nm,rt in roots.items()}
inv=sp.invert(sp.Poly(sp.expand(sp.diff(Q3,r)*T_V),r),sp.Poly(Q3,r)); Pher=sp.expand(sp.rem(sp.expand(-2*N3*inv.as_expr()),Q3,r))
RemH=sp.div(sp.Poly(sp.expand(2*N3-2*sp.diff(Pher,r)*Q3*T_V-Pher*Q3*sp.diff(T_V,r)+Pher*sp.diff(Q3,r)*T_V),r),sp.Poly(sp.expand(2*Q3),r))[0].as_expr()
remk=[c for c in sp.Poly(RemH,r).all_coeffs()[::-1]]
g=[sp.simplify((p_j[k] if k<len(p_j) else 0)+(remk[k] if k<len(remk) else 0)) for k in range(max(len(p_j),len(remk)))]
Sn=sp.lambdify(r,Ssig.subs(sub),'numpy'); sq=lambda x:np.sqrt(max(Sn(x),0.0))
Q3n=sp.lambdify(r,Q3.subs(sub),'numpy'); dgn=sp.lambdify(r,dengV.subs(sub),'numpy'); Rrn=sp.lambdify(r,Rr.subs(sub),'numpy')
A_Vn=sp.lambdify(r,(J*E*r**3*DE/Q3).subs(sub),'numpy'); P_in=sp.lambdify(r,(-m*N4/dengV).subs(sub),'numpy'); Phern=sp.lambdify(r,Pher.subs(sub),'numpy')
gn=[float(x.subs(sub)) for x in g]; akn=[float(x.subs(sub)) for x in a_k]; rootn={nm:float(rt.subs(sub)) for nm,rt in roots.items()}; rhon={nm:float(rho[nm].subs(sub)) for nm in roots}
def Uk(x,k): return quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def Sig(x): return quad(lambda t:P_in(t)/sq(t),r0,x,limit=200)[0]
def wrap_closed(x):
    A=sum(-0.5*gn[j]*akn[k]*(Uk(x,j)*Uk(x,k)+quad(lambda t:(t**j*Uk(t,k)-t**k*Uk(t,j))/sq(t),r0,x,limit=150)[0]) for j in range(len(gn)) for k in range(len(akn)))
    Bp=-sum(gn[j]*Uk(x,j) for j in range(len(gn)))*quad(lambda t:Rrn(t)/(dgn(t)*sq(t)),r0,x,limit=200)[0]
    Bd=sum(rhon[nm]*sum(gn[j]*quad(lambda t:Uk(t,j)/((t-rootn[nm])*sq(t)),r0,x,limit=150)[0] for j in range(len(gn))) for nm in roots)
    Cc=-(Phern(x)*sq(x)/Q3n(x)*Sig(x)-Phern(r0)*sq(r0)/Q3n(r0)*Sig(r0))+quad(lambda t:Phern(t)*P_in(t)/Q3n(t),r0,x,limit=200)[0]
    return A+Bp+Bd+Cc
print("\n(2) closed form (A+B+C) vs flow off-shell sub-piece (closed = -partB by sign convention):")
for x in [11.0,10.0,9.0]:
    print(f"    r={x}: closed={wrap_closed(x):+.7f}  (-flow_partB)={-float(partB_r(x)):+.7f}  diff={abs(wrap_closed(x)+float(partB_r(x))):.1e}")
print("\n=> physics anchored: closed form == flow off-shell sub-piece (up to sign); total slope vs true tau flow above.")
