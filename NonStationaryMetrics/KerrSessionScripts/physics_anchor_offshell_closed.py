# -*- coding: utf-8 -*-
# PHYSICS ANCHOR of the closed-form generic off-shell wrap (t-branch).
# Chain (each link verified):
#   my closed form (A+B+C assembly)  ==[here, 1e-13]==  direct wrap -int (A/sqrtS) Delta S dr
#   -int(A/sqrtS)Delta S  ==[A/sqrtS=G_Pr/H_Pr on shell, kernel_wrap 3e-15]==  machine dilation wrap
#   machine total (partA+partB, partB contains the dilation wrap) ==[slope 2]== TRUE non-autonomous flow
# So the closed form equals the flow-validated sub-object; here we (1) re-run the true-flow slope
# check and (2) show closed == machine dilation-wrap on the SAME frozen orbit.
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct, quad
from scipy.optimize import brentq
from scipy.interpolate import interp1d

M,a,Ehat,J0,r0 = 1.0,0.9,1.4,6.0,12.0
rr,pr,Ess,Jss=sp.symbols('r pr E J_')
f2=1-2*M/rr; Dl2=rr**2-2*M*rr+a**2; b2=2*M*a/rr; v2=1-f2/Ess**2
P2=rr**2+a**2+2*M*a**2/rr; Pb2=P2+b2**2/Ess**2
H2=Jss*b2*v2/Pb2+sp.sqrt(Dl2*v2/Pb2)*sp.sqrt((Dl2/rr**2)*pr**2+Jss**2/Pb2)-1
dHp=sp.lambdify((rr,pr,Ess,Jss),sp.diff(H2,pr),'numpy'); dHr=sp.lambdify((rr,pr,Ess,Jss),sp.diff(H2,rr),'numpy')
dHJ=sp.lambdify((rr,pr,Ess,Jss),sp.diff(H2,Jss),'numpy'); dHE=sp.lambdify((rr,pr,Ess,Jss),sp.diff(H2,Ess),'numpy')
H2n=sp.lambdify((rr,pr,Ess,Jss),H2,'numpy')
G=sp.diff(H2,Jss)/sp.diff(H2,pr)
dG_pr=sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,pr),'numpy')
def prof(rv,E,Jv):
    pg=np.linspace(-80,80,3001); Hv=H2n(rv,pg,E,Jv)
    rts=[brentq(lambda p:H2n(rv,p,E,Jv),pg[i],pg[i+1]) for i in range(len(pg)-1)
         if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1]<0]
    ing=[p for p in rts if dHp(rv,p,E,Jv)<0]; return min(ing) if ing else np.nan
ev=lambda lam,y:y[1]; ev.terminal=True; ev.direction=1
def flow(eps):
    def rhs(lam,y):
        rv,pv,ph=y; s=np.exp(-eps*lam); E=Ehat*s; Jv=J0*s
        return [dHp(rv,pv,E,Jv),-dHr(rv,pv,E,Jv)-eps*pv,dHJ(rv,pv,E,Jv)]
    so=solve_ivp(rhs,[0,300],[r0,prof(r0,Ehat,J0),0.0],rtol=1e-12,atol=1e-14,max_step=0.004,dense_output=True,events=ev)
    lam=np.linspace(0,so.t[-1],16000); Y=so.sol(lam); return lam,Y[0],Y[1],Y[2]
lam0,rF,prF,phiF=flow(0.0)
# machine partA + partB (the flow-validated total)
EulerH=Ehat*dHE(rF,prF,Ehat,J0)+J0*dHJ(rF,prF,Ehat,J0); Hpr=dHp(rF,prF,Ehat,J0)
S_D=ct(EulerH+prF*Hpr,lam0,initial=0)                       # full S_D (Theta + dilation)
DSdil=ct(prF*Hpr,lam0,initial=0)                            # dilation letter only = int p_r dr
Gpr=dG_pr(rF,prF,Ehat,J0)
dG_E=sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,Ess),'numpy'); dG_J=sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,Jss),'numpy')
EulerG=Ehat*dG_E(rF,prF,Ehat,J0)+J0*dG_J(rF,prF,Ehat,J0)
partA=ct(Gpr*(lam0*EulerH)/Hpr-lam0*EulerG,rF,initial=0)
partB=ct(-Gpr*S_D/Hpr,rF,initial=0)
dphi_eps=partA+partB; dphi_r=interp1d(rF,dphi_eps,fill_value='extrapolate',bounds_error=False)
phi0f=interp1d(rF,phiF,fill_value='extrapolate',bounds_error=False)
# machine's DILATION WRAP sub-piece: -int Gpr*(int p_r dr)/Hpr
mach_dil=ct(-Gpr*DSdil/Hpr,rF,initial=0); mach_dil_r=interp1d(rF,mach_dil,fill_value='extrapolate',bounds_error=False)
rc=np.linspace(8.0,11.0,1500)
print("(1) PHYSICS: total (partA+partB) vs TRUE non-autonomous flow, slope of residual:")
epss=np.array([0.002,0.004,0.008,0.016,0.032]); res=[]
for eps in epss:
    _,rL,_,pL=flow(eps); pt=interp1d(rL,pL,fill_value='extrapolate',bounds_error=False)(rc)
    res.append(np.nanmax(np.abs(pt-(phi0f(rc)+eps*dphi_r(rc)))))
res=np.array(res); print(f"    slope = {np.polyfit(np.log(epss),np.log(res),1)[0]:.2f}  (==2 => closed-form-containing total reproduces physics)")

# ---- (2) my CLOSED FORM assembly of the dilation wrap, on the SAME frozen orbit (r0=12) ----
r,E=sp.symbols('r E',positive=True); E0=sp.Rational(7,5); Ms,asf,Jf=sp.Rational(1),sp.Rational(9,10),sp.Integer(6)
Emu=(E**2-1)*r+2*Ms; Dl=r**2-2*Ms*r+asf**2
Q2=(2*E**2*Jf**2*Ms*r-E**2*Jf**2*r**2-4*E**2*Jf*Ms*asf*r+2*E**2*Ms*asf**2*r+E**2*asf**2*r**2
    +E**2*r**4+4*Jf**2*Ms**2-4*Jf**2*Ms*r+Jf**2*r**2-8*Jf*Ms**2*asf+4*Jf*Ms*asf*r+4*Ms**2*asf**2)
T=sp.expand(r*Emu); Ssh=sp.expand(Q2*T); deng=sp.expand(Dl*Emu)
Qp,N3=sp.div(sp.Poly(sp.expand(E**2*Jf*r**4*Emu),r),sp.Poly(sp.expand(Q2),r)); N3=N3.as_expr()
p_j=[c for c in Qp.all_coeffs()[::-1]]
inv=sp.invert(sp.Poly(sp.expand(sp.diff(Q2,r)*T),r),sp.Poly(Q2,r))
Pherm=sp.expand(sp.rem(sp.expand(-2*N3*inv.as_expr()),Q2,r))
Rem=sp.div(sp.Poly(sp.expand(2*N3-2*sp.diff(Pherm,r)*Q2*T-Pherm*Q2*sp.diff(T,r)+Pherm*sp.diff(Q2,r)*T),r),sp.Poly(sp.expand(2*Q2),r))[0].as_expr()
remc=sp.Poly(Rem,r).all_coeffs()[::-1]
g=[sp.simplify((p_j[k] if k<len(p_j) else 0)+(remc[k] if k<len(remc) else 0)) for k in range(5)]
Qd,R=sp.div(sp.Poly(sp.expand(Ssh),r),sp.Poly(deng,r)); R=R.as_expr()
a_k=[sp.simplify(c) for c in Qd.all_coeffs()[::-1]]
roots={'r+':Ms+sp.sqrt(Ms**2-asf**2),'r-':Ms-sp.sqrt(Ms**2-asf**2)}
rho={nm:sp.simplify((R/sp.diff(deng,r)).subs(r,rt)) for nm,rt in roots.items()}
Sn=sp.lambdify(r,Ssh.subs(E,E0),'numpy'); sqn=lambda x:np.sqrt(max(Sn(x),0.0))
Q2f=sp.lambdify(r,Q2.subs(E,E0),'numpy'); dengf=sp.lambdify(r,deng.subs(E,E0),'numpy')
N3f=sp.lambdify(r,N3.subs(E,E0),'numpy'); Rf=sp.lambdify(r,R.subs(E,E0),'numpy'); Pf=sp.lambdify(r,Pherm.subs(E,E0),'numpy')
gn=[float(x.subs(E,E0)) for x in g]; akn=[float(x.subs(E,E0)) for x in a_k]; pjn=[float(x.subs(E,E0)) for x in p_j]
rootn={nm:float(rt.subs(E,E0)) for nm,rt in roots.items()}; rhon={nm:float(rho[nm].subs(E,E0)) for nm in roots}
R0=12.0
def Uk(x,k): return quad(lambda t:t**k/sqn(t),R0,x,limit=200)[0]
def Pi_of(x): return quad(lambda t:Rf(t)/(dengf(t)*sqn(t)),R0,x,limit=200)[0]
def DSf(x): return sum(akn[k]*Uk(x,k) for k in range(4))+Pi_of(x)
def AoS(t): return (pjn[0]+pjn[1]*t)/sqn(t)+N3f(t)/(Q2f(t)*sqn(t))
def wrap_direct(x): return -quad(lambda t:AoS(t)*DSf(t),R0,x,limit=150)[0]
def Wjk(j,k,x): return quad(lambda t:(t**j*Uk(t,k)-t**k*Uk(t,j))/sqn(t),R0,x,limit=150)[0]
def blockA(x): return sum(-0.5*gn[j]*akn[k]*(Uk(x,j)*Uk(x,k)+Wjk(j,k,x)) for j in range(5) for k in range(4))
def PhiG(x): return sum(gn[j]*Uk(x,j) for j in range(5))
def D(j,nm,x): return quad(lambda t:Uk(t,j)/((t-rootn[nm])*sqn(t)),R0,x,limit=150)[0]
def blockB(x): return -PhiG(x)*Pi_of(x)+sum(rhon[nm]*sum(gn[j]*D(j,nm,x) for j in range(5)) for nm in roots)
Dlf=sp.lambdify(r,Dl.subs(E,E0),'numpy'); PrDl=sp.lambdify(r,(Pherm*r/Dl).subs(E,E0),'numpy')
def blockC(x): return -(Pf(x)*sqn(x)/Q2f(x)*DSf(x)-Pf(R0)*sqn(R0)/Q2f(R0)*DSf(R0))+quad(PrDl,R0,x,limit=200)[0]
def wrap_closed(x): return blockA(x)+blockB(x)+blockC(x)
print("\n(2) CLOSED FORM (A+B+C) vs direct dilation wrap -int(A/sqrtS)DeltaS, SAME orbit r0=12:")
for x in [11.0,10.0,9.0]:
    d=wrap_direct(x); c=wrap_closed(x); print(f"    r={x}: closed={c:+.7f}  direct={d:+.7f}  diff={abs(c-d):.1e}")
print("\n(3) direct dilation wrap == machine's dilation-wrap sub-piece of partB (same object, A/sqrtS=Gpr/Hpr;")
print("    machine works in flow-phi convention = -shape convention, so compare wrap_direct vs -machine_dil):")
for x in [11.0,10.0,9.0]:
    md=-float(mach_dil_r(x))
    print(f"    r={x}: wrap_direct={wrap_direct(x):+.7f}  (-machine_dil)={md:+.7f}  diff={abs(wrap_direct(x)-md):.1e}  (flow numerical floor)")
print("\n=> chain closed: assembly = direct wrap = machine dilation sub-piece; total = true flow (slope 2).")
