# -*- coding: utf-8 -*-
# Vaidya analogue of the exact first-order adiabatic correction (referee follow-up 4.6).
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

# ---- convergence figure (referee follow-up 4.6, fig:vaidya-offshell) ----
fig,ax=plt.subplots(1,1,figsize=(COL,COL*0.85))
ax.loglog(epss,res_on,'ko-',ms=4,label=r'on-shell only: slope $\approx%.2f$'%s_on)
ax.loglog(epss,res,'C3s-',ms=4,label=r'complete Eq.(40): slope $\approx%.2f$'%s_full)
ax.loglog(epss,epss**2*res[2]/epss[2]**2,'C3:',lw=0.6,alpha=0.6,label=r'$O(\varepsilon^2)$ guide')
ax.loglog(epss,epss*res_on[2]/epss[2],'k:',lw=0.6,alpha=0.6,label=r'$O(\varepsilon)$ guide')
ax.set_xlabel(r"$\varepsilon=M\dot m/m$"); ax.set_ylabel('residual vs true flow')
ax.set_title('Vaidya first order: off-shell costate closes to $O(\\varepsilon^2)$',fontsize=6.6)
ax.legend(fontsize=5.6,loc='upper left')
savefig(fig,os.path.join(os.path.dirname(HERE),'paper','Immagini'),'fig_vaidya_offshell')
print("FATTO")
