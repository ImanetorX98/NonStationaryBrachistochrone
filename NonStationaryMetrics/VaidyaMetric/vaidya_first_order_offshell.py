# -*- coding: utf-8 -*-
# Vaidya analogue of the exact first-order adiabatic correction.
# Same canonical perturbation theory as Thakurta Eq.(40), single slow parameter m(v):
#   Theta = m d_m,  S = int_0^lambda Theta H_v dlambda,  H_v the v-branch Hamiltonian (Eq 11),
#   delta p_r/eps = (S - lambda Theta H)/H_pr  [sign: m INCREASES, opposite to Thakurta E,J],
#   delta phi/eps = int[ G_pr (S - lambda Theta H)/H_pr + lambda Theta G ] dr,  G=H_J/H_pr.
# Off-shell costate: H_ext=p_v+H_v=0 => p_v=-H_v=O(mdot). Verified: residual slope ~2.00
# (closes the Vaidya first order to O(eps^2)), so Eqs (17-21) on-shell horizon-dilogarithm
# are the on-shell COMPONENT; the complete first-order term adds this off-shell piece.
import os,sys
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct
from scipy.optimize import brentq
from scipy.interpolate import interp1d
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
set_style(); HERE=os.path.dirname(os.path.abspath(__file__))
E0,J0,r0=1.4,6.0,12.0; m0=1.0
rr,pr,mm=sp.symbols('r pr m',positive=True); Es,Js=sp.symbols('E J_',positive=True)
f=1-2*mm/rr; w=Es**2-f
Hv=pr*(f-Es**2)-1+sp.sqrt(w)*sp.sqrt(Es**2*pr**2+Js**2/rr**2)   # v-branch, Eq (11)
Hvn=sp.lambdify((rr,pr,mm,Es,Js),Hv,'numpy')
Hp=sp.diff(Hv,pr); Hr=sp.diff(Hv,rr); HJ=sp.diff(Hv,Js); Hm=sp.diff(Hv,mm); G=HJ/Hp
dHp=sp.lambdify((rr,pr,mm,Es,Js),Hp,'numpy'); dHr=sp.lambdify((rr,pr,mm,Es,Js),Hr,'numpy')
dHJ=sp.lambdify((rr,pr,mm,Es,Js),HJ,'numpy'); dHm=sp.lambdify((rr,pr,mm,Es,Js),Hm,'numpy')
dGpr=sp.lambdify((rr,pr,mm,Es,Js),sp.diff(G,pr),'numpy'); dGm=sp.lambdify((rr,pr,mm,Es,Js),sp.diff(G,mm),'numpy')
def prof(rv,mv):
    pg=np.linspace(-80,80,4001); Hval=Hvn(rv,pg,mv,E0,J0)
    rts=[brentq(lambda p:Hvn(rv,p,mv,E0,J0),pg[i],pg[i+1]) for i in range(len(pg)-1)
         if np.isfinite(Hval[i]) and np.isfinite(Hval[i+1]) and Hval[i]*Hval[i+1]<0]
    ing=[p for p in rts if dHp(rv,p,mv,E0,J0)<0]; return min(ing) if ing else np.nan
ev=lambda lam,y:y[1]; ev.terminal=True; ev.direction=1
def flow(eps):
    def rhs(lam,y):
        rv,pv,ph=y; mv=m0*np.exp(eps*lam)
        return [dHp(rv,pv,mv,E0,J0),-dHr(rv,pv,mv,E0,J0),dHJ(rv,pv,mv,E0,J0)]
    so=solve_ivp(rhs,[0,300],[r0,prof(r0,m0),0.0],rtol=1e-12,atol=1e-14,max_step=0.005,dense_output=True,events=ev)
    lam=np.linspace(0,so.t[-1],12000); Y=so.sol(lam); return lam,Y[0],Y[1],Y[2]
lam0,rF,prF,phiF=flow(0.0)
phi0=interp1d(rF,phiF,bounds_error=False,fill_value='extrapolate')
lam_r=interp1d(rF,lam0,bounds_error=False,fill_value='extrapolate')
rg=np.linspace(r0,rF.min()+0.3,5000); rc=np.linspace(8.0,11.0,2000)
prg=interp1d(rF,prF,bounds_error=False,fill_value='extrapolate')(rg)
ThetaH=m0*dHm(rg,prg,m0,E0,J0)                    # Theta_m H = m dH/dm
S=ct(ThetaH*np.gradient(lam_r(rg),rg),rg,initial=0)
Gpr=dGpr(rg,prg,m0,E0,J0)/dHp(rg,prg,m0,E0,J0); Gm=lam_r(rg)*m0*dGm(rg,prg,m0,E0,J0)
integ=Gpr*(S-lam_r(rg)*ThetaH)+Gm                 # complete first order (on- + off-shell)
integ_on=Gpr*(-lam_r(rg)*ThetaH)+Gm               # ON-SHELL only (drop costate S)
xc=interp1d(rg,ct(integ,rg,initial=0),bounds_error=False,fill_value='extrapolate')
xc_on=interp1d(rg,ct(integ_on,rg,initial=0),bounds_error=False,fill_value='extrapolate')
# leading (on-shell, m d_m F * clock) for comparison: Theta_m applied to shape
print("Vaidya v-branch off-shell PT:")
epss=np.array([0.0025,0.005,0.01,0.02,0.04]); res=[]; res_on=[]
for eps in epss:
    _,rL,_,pL=flow(eps); pt=interp1d(rL,pL,bounds_error=False,fill_value='extrapolate')(rc)
    res.append(np.nanmax(np.abs(pt-(phi0(rc)+eps*xc(rc)))))
    res_on.append(np.nanmax(np.abs(pt-(phi0(rc)+eps*xc_on(rc)))))
res=np.array(res); res_on=np.array(res_on)
print("  coeff on [8,11]: %.4f"%np.mean(xc(rc)))
for e,r,ro in zip(epss,res,res_on): print(f"  eps={e:.4f}  res_full={r:.2e}  res_onshell={ro:.2e}")
s_full=np.polyfit(np.log(epss),np.log(res),1)[0]
s_on=np.polyfit(np.log(epss),np.log(res_on),1)[0]
print(f"  SLOPE full = {s_full:.2f}  (2 => off-shell PT closes Vaidya too)")
print(f"  SLOPE on-shell only = {s_on:.2f}  (1 => costate piece needed)")

# ---- convergence figure (fig:vaidya-offshell) ----
fig,ax=plt.subplots(1,1,figsize=(COL,COL*0.85))
ax.loglog(epss,res_on,'ko-',ms=4,label=r'on-shell only: slope $\approx%.2f$'%s_on)
ax.loglog(epss,res,'C3s-',ms=4,label=r'complete Eq.(40): slope $\approx%.2f$'%s_full)
ax.loglog(epss,epss**2*res[2]/epss[2]**2,'C3:',lw=0.6,alpha=0.6,label=r'$O(\varepsilon^2)$ guide')
ax.loglog(epss,epss*res_on[2]/epss[2],'k:',lw=0.6,alpha=0.6,label=r'$O(\varepsilon)$ guide')
ax.set_xlabel(r"$\varepsilon=M\dot m/m$"); ax.set_ylabel('residual vs true flow')
ax.set_title('Vaidya first order: off-shell costate closes to $O(\\varepsilon^2)$',fontsize=6.6)
ax.legend(fontsize=5.6,loc='upper left')
savefig(fig,os.path.join(os.path.dirname(HERE),'paper','Immagini'),'fig_vaidya_offshell')

# ---- the proper-time branch H_tau (eq:Htau-vaidya) ----
# The on-shell closed form of Paper I section 4.4 is written on H_tau, not on H_v, so the
# complete first-order term is validated here on that branch too.  Same perturbation
# theory, S = int Theta H_tau dlambda (anchored at the launch); with the mass increasing
# along lambda the off-shell term enters with the + sign used above.
Ht=pr*(f-Es**2)-Es+sp.sqrt(w)*sp.sqrt((Es*pr+1)**2+Js**2/rr**2)
L5=lambda X: sp.lambdify((rr,pr,mm,Es,Js),X,'numpy')
Htn,tHp,tHr,tHJ,tHm=L5(Ht),L5(sp.diff(Ht,pr)),L5(sp.diff(Ht,rr)),L5(sp.diff(Ht,Js)),L5(sp.diff(Ht,mm))
Gt=sp.diff(Ht,Js)/sp.diff(Ht,pr); tGpr,tGm=L5(sp.diff(Gt,pr)),L5(sp.diff(Gt,mm))
Et,Jt=1.4,2.5
def prof_t(rv,mv):
    pg=np.linspace(-80,80,4001); Hval=Htn(rv,pg,mv,Et,Jt)
    rts=[brentq(lambda p:Htn(rv,p,mv,Et,Jt),pg[i],pg[i+1]) for i in range(len(pg)-1)
         if np.isfinite(Hval[i]) and np.isfinite(Hval[i+1]) and Hval[i]*Hval[i+1]<0]
    ing=[p for p in rts if tHp(rv,p,mv,Et,Jt)<0]; return min(ing) if ing else np.nan
def flow_t(eps):
    def rhs(lam,y):
        rv,pv,ph=y; mv=m0*np.exp(eps*lam)
        return [tHp(rv,pv,mv,Et,Jt),-tHr(rv,pv,mv,Et,Jt),tHJ(rv,pv,mv,Et,Jt)]
    so=solve_ivp(rhs,[0,300],[r0,prof_t(r0,m0),0.0],rtol=1e-12,atol=1e-14,max_step=0.005,dense_output=True,events=ev)
    lam=np.linspace(0,so.t[-1],12000); Y=so.sol(lam); return lam,Y[0],Y[1],Y[2]
lt,rT,pT,phT=flow_t(0.0)
phi0t=interp1d(rT,phT,bounds_error=False,fill_value='extrapolate')
lam_t=interp1d(rT,lt,bounds_error=False,fill_value='extrapolate')
rmin_t=rT.min(); rgt=np.linspace(r0,rmin_t+0.3,5000)
rct=np.linspace(rmin_t+0.35*(r0-rmin_t),rmin_t+0.80*(r0-rmin_t),2000)
prt=interp1d(rT,pT,bounds_error=False,fill_value='extrapolate')(rgt)
ThH=m0*tHm(rgt,prt,m0,Et,Jt)
St=ct(ThH*np.gradient(lam_t(rgt),rgt),rgt,initial=0)
Gp_t=tGpr(rgt,prt,m0,Et,Jt)/tHp(rgt,prt,m0,Et,Jt); Gm_t=lam_t(rgt)*m0*tGm(rgt,prt,m0,Et,Jt)
xt=interp1d(rgt,ct(Gp_t*(St-lam_t(rgt)*ThH)+Gm_t,rgt,initial=0),bounds_error=False,fill_value='extrapolate')
xt_on=interp1d(rgt,ct(Gp_t*(-lam_t(rgt)*ThH)+Gm_t,rgt,initial=0),bounds_error=False,fill_value='extrapolate')
rt_=[];rto=[]
for eps in epss:
    _,rL,_,pL=flow_t(eps); pt=interp1d(rL,pL,bounds_error=False,fill_value='extrapolate')(rct)
    rt_.append(np.nanmax(np.abs(pt-(phi0t(rct)+eps*xt(rct)))))
    rto.append(np.nanmax(np.abs(pt-(phi0t(rct)+eps*xt_on(rct)))))
st_full=np.polyfit(np.log(epss),np.log(rt_),1)[0]; st_on=np.polyfit(np.log(epss),np.log(rto),1)[0]
print(f"Vaidya tau-branch (E={Et}, J={Jt}): SLOPE full = {st_full:.2f}, on-shell only = {st_on:.2f}")

# ---- normalisation test at a frozen mass different from 1 ----
# The source uses Theta = m d_m, so the first-order coefficient multiplies the RELATIVE rate
# q = mdot/m, not mdot.  With m0 = 1 the two coincide and the factor is invisible; at m0 = 2
# only q closes the correction to O(q^2).  (M = 1 units; E = 1.4, J = 5, r0 = 24, compact
# window r in [16, 22] on the monotone ingoing arc.)
m2, J2, r02 = 2.0, 5.0, 24.0
def prof2(rv):
    pg=np.linspace(-80,80,4001); Hval=Htn(rv,pg,m2,Et,J2)
    rts=[brentq(lambda p:Htn(rv,p,m2,Et,J2),pg[i],pg[i+1]) for i in range(len(pg)-1)
         if np.isfinite(Hval[i]) and np.isfinite(Hval[i+1]) and Hval[i]*Hval[i+1]<0]
    return min(p for p in rts if tHp(rv,p,m2,Et,J2)<0)
ev15=lambda lam,y:y[0]-15.0; ev15.terminal=True
def flow2(q):
    def rhs(lam,y):
        rv,pv,ph=y; mv=m2*np.exp(q*lam)
        return [tHp(rv,pv,mv,Et,J2),-tHr(rv,pv,mv,Et,J2),tHJ(rv,pv,mv,Et,J2)]
    so=solve_ivp(rhs,[0,500],[r02,prof2(r02),0.0],rtol=1e-12,atol=1e-14,max_step=0.01,dense_output=True,events=ev15)
    lam=np.linspace(0,so.t[-1],12000); Y=so.sol(lam); return lam,Y[0],Y[1],Y[2]
l2,r2_,p2_,ph2=flow2(0.0)
phi02=interp1d(r2_,ph2); lam2=interp1d(r2_,l2)
rg2=np.linspace(r02,15.2,5000); rc2=np.linspace(16.0,22.0,1500)
pr2_=interp1d(r2_,p2_)(rg2)
Th2=m2*tHm(rg2,pr2_,m2,Et,J2)                        # Theta H = m d_m H
S2=ct(Th2*np.gradient(lam2(rg2),rg2),rg2,initial=0)
Gp2=tGpr(rg2,pr2_,m2,Et,J2)/tHp(rg2,pr2_,m2,Et,J2); Gm2=lam2(rg2)*m2*tGm(rg2,pr2_,m2,Et,J2)
on2=interp1d(rg2,ct(Gp2*(-lam2(rg2)*Th2)+Gm2,rg2,initial=0))     # on-shell part, per unit q
off2=interp1d(rg2,ct(Gp2*S2,rg2,initial=0))                      # off-shell part, per unit q
qs=np.array([0.001,0.002,0.004,0.008]); res_q=[]; res_md=[]
for q in qs:
    _,rL,_,pL=flow2(q); pt=interp1d(rL,pL)(rc2)
    res_q.append(np.max(np.abs(pt-(phi02(rc2)+q*(on2(rc2)+off2(rc2))))))
    res_md.append(np.max(np.abs(pt-(phi02(rc2)+q*on2(rc2)+m2*q*off2(rc2)))))   # prefactor mdot = m q
sq_=np.polyfit(np.log(qs),np.log(res_q),1)[0]; smd=np.polyfit(np.log(qs),np.log(res_md),1)[0]
print(f"m0={m2}: off-shell prefactor mdot/m -> slope {sq_:.2f};  prefactor mdot -> slope {smd:.2f}")
assert sq_>1.8 and smd<1.2, "mass-normalisation test failed"
print("FATTO")
