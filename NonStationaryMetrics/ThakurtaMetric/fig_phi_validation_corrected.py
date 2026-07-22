# -*- coding: utf-8 -*-
"""
Honest validation figure for the adiabatic correction (referee 4.6/4.7/4.14, follow-up):
 Left  : true non-autonomous optical-metric geodesic (H2 flow, E_eff=Ehat/A, J_eff=J/A)
         vs frozen + leading 1/2-Euler AND vs frozen + the exact Eq.(40) correction, A'/A=0.04.
 Right : residual vs A'/A -- the leading 1/2-Euler is O(eps) (~2% physical error), the exact
         Eq.(40) term is O(eps^2) (slope ~2), both far above the algebraic IBP floor (~1e-6).
"""
import os,sys
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct
from scipy.optimize import brentq
from scipy.interpolate import interp1d
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
set_style(); HERE=os.path.dirname(os.path.abspath(__file__))
M,a,Ehat=1.0,0.9,1.4; J0,r0=6.0,12.0

r,Es,Js=sp.symbols('r E J_'); Dl=r**2-2*M*r+a**2
Q2=(2*Es**2*Js**2*M*r-Es**2*Js**2*r**2-4*Es**2*Js*M*a*r+2*Es**2*M*a**2*r+Es**2*a**2*r**2
    +Es**2*r**4+4*Js**2*M**2-4*Js**2*M*r+Js**2*r**2-8*Js*M**2*a+4*Js*M*a*r+4*M**2*a**2)
R=r*Q2*((Es**2-1)*r+2*M); Kf=r*((Es**2-1)*r+2*M)*(Js*(r-2*M)+2*M*a)/Dl; F=Kf/sp.sqrt(R)
dEF=sp.lambdify((r,Es,Js),sp.diff(F,Es),'numpy'); dJF=sp.lambdify((r,Es,Js),sp.diff(F,Js),'numpy')
rr,pr,Ess,Jss=sp.symbols('r pr E J_')
f2=1-2*M/rr; Dl2=rr**2-2*M*rr+a**2; b2=2*M*a/rr; v2=1-f2/Ess**2
P2=rr**2+a**2+2*M*a**2/rr; Pb2=P2+b2**2/Ess**2
H2=Jss*b2*v2/Pb2+sp.sqrt(Dl2*v2/Pb2)*sp.sqrt((Dl2/rr**2)*pr**2+Jss**2/Pb2)-1
H2n=sp.lambdify((rr,pr,Ess,Jss),H2,'numpy')
Hp=sp.diff(H2,pr); HJ=sp.diff(H2,Jss); HE=sp.diff(H2,Ess); G=HJ/Hp
dHp=sp.lambdify((rr,pr,Ess,Jss),Hp,'numpy'); dHr=sp.lambdify((rr,pr,Ess,Jss),sp.diff(H2,rr),'numpy')
dHJ=sp.lambdify((rr,pr,Ess,Jss),HJ,'numpy'); dHE=sp.lambdify((rr,pr,Ess,Jss),HE,'numpy')
dGpr=sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,pr),'numpy'); dGE=sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,Ess),'numpy')
dGJ=sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,Jss),'numpy')
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
    so=solve_ivp(rhs,[0,300],[r0,prof(r0,Ehat,J0),0.0],rtol=1e-12,atol=1e-14,max_step=0.005,dense_output=True,events=ev)
    lam=np.linspace(0,so.t[-1],12000); Y=so.sol(lam); return lam,Y[0],Y[1],Y[2]
lam0,rF,prF,phiF=flow(0.0)
phi0=interp1d(rF,phiF,bounds_error=False,fill_value='extrapolate')
lam_r=interp1d(rF,lam0,bounds_error=False,fill_value='extrapolate')
rg=np.linspace(r0,7.3,5000); rc=np.linspace(8.0,11.0,2000)
eul=-0.5*ct(lam_r(rg)*(Ehat*dEF(rg,Ehat,J0)+J0*dJF(rg,Ehat,J0)),rg,initial=0)
ec=interp1d(rg,eul,bounds_error=False,fill_value='extrapolate')
prg=interp1d(rF,prF,bounds_error=False,fill_value='extrapolate')(rg)
EulerH=Ehat*dHE(rg,prg,Ehat,J0)+J0*dHJ(rg,prg,Ehat,J0)
S=ct(EulerH*np.gradient(lam_r(rg),rg),rg,initial=0)
integ=dGpr(rg,prg,Ehat,J0)*(lam_r(rg)*EulerH-S)/dHp(rg,prg,Ehat,J0)-lam_r(rg)*(Ehat*dGE(rg,prg,Ehat,J0)+J0*dGJ(rg,prg,Ehat,J0))
xc=interp1d(rg,ct(integ,rg,initial=0),bounds_error=False,fill_value='extrapolate')

fig,ax=plt.subplots(1,2,figsize=(2*COL,COL*0.9))
eps=0.04; _,rL,_,pL=flow(eps); pt=interp1d(rL,pL,bounds_error=False,fill_value='extrapolate')
rr2=np.linspace(r0,rF.min()+0.4,600)
ax[0].plot(rr2*np.cos(pt(rr2)),rr2*np.sin(pt(rr2)),'C0-',lw=2.6,alpha=0.35,label='true non-autonomous geodesic')
ph_h=phi0(rr2)+eps*ec(rr2); ph_x=phi0(rr2)+eps*xc(rr2)
ax[0].plot(rr2*np.cos(ph_h),rr2*np.sin(ph_h),'k--',lw=1.0,label=r'frozen $+\,\varepsilon\cdot$leading (on-shell)')
ax[0].plot(rr2*np.cos(ph_x),rr2*np.sin(ph_x),'C3:',lw=1.3,label=r'frozen $+\,\varepsilon\cdot$Eq.(40) (exact)')
ax[0].set_aspect('equal'); ax[0].set_xlabel('$x$'); ax[0].set_ylabel('$y$')
ax[0].set_title(r'true dynamics vs adiabatic ($t$-branch, $A^\prime/A=0.04$)',fontsize=6.6)
ax[0].legend(fontsize=5.4,loc='upper left',framealpha=0.9)
epss=np.array([0.0025,0.005,0.01,0.02,0.04]); rh=[];rx=[]
for e in epss:
    _,rL,_,pL=flow(e); p=interp1d(rL,pL,bounds_error=False,fill_value='extrapolate')(rc)
    rh.append(np.nanmax(np.abs(p-(phi0(rc)+e*ec(rc))))); rx.append(np.nanmax(np.abs(p-(phi0(rc)+e*xc(rc)))))
rh=np.array(rh);rx=np.array(rx)
ax[1].loglog(epss,rh,'ko-',ms=4,label=r'leading (on-shell): slope $\approx1$ ($\sim2\%$)')
ax[1].loglog(epss,rx,'C3s-',ms=4,label=r'exact Eq.(40): slope $\approx2$')
ax[1].loglog(epss,3e-4*epss,'C7:',lw=0.8,label=r'algebraic IBP floor (Fig.~9c)')
ax[1].loglog(epss,epss**2*rx[2]/epss[2]**2,'C3:',lw=0.5,alpha=0.6)
ax[1].set_xlabel(r"$A^\prime/A$"); ax[1].set_ylabel('residual')
ax[1].set_title('exact term closes the physical error to $O(\\varepsilon^2)$',fontsize=6.6)
ax[1].legend(fontsize=5.4,loc='upper left')
savefig(fig,os.path.join(os.path.dirname(HERE),'paper','Immagini'),'fig_phi_validation_true_dynamic')
print("slopes: leading=%.2f exact=%.2f"%(np.polyfit(np.log(epss),np.log(rh),1)[0],np.polyfit(np.log(epss),np.log(rx),1)[0]))
print("FATTO")
