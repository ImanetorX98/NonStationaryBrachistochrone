# -*- coding: utf-8 -*-
"""
Honest validation figure for the adiabatic correction (referee 4.6/4.7/4.14 and main10 Fig.10):
 Left  : true non-autonomous optical-metric geodesic (frozen E_eff=Ehat/A, J_eff=J/A)
         vs frozen + leading 1/2-Euler AND vs frozen + the exact Eq.(40) correction (t-branch).
 Right : residual vs A'/A for BOTH branches (t and tau) and a SECOND (a,E,J) parameter set.
         The leading 1/2-Euler is O(eps) (slope ~1, ~2% physical error); the exact Eq.(40)
         term is O(eps^2) (slope ~2). Slopes are reported with a least-squares 1sigma
         uncertainty (error bars). The tau-branch Hamiltonian H_tau is verified elsewhere to
         reproduce the Boyer-Lindquist tau shape (kerr_adiabatic_phi_hybrid_tau.py) exactly.
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

rr,pr,Ess,Jss=sp.symbols('r pr E J_')

def H_eta(M,a):
    """conformal/optical (t-family) frozen Hamiltonian H_eta (paper Eq.(29))."""
    f=1-2*M/rr; Dl=rr**2-2*M*rr+a**2; b=2*M*a/rr; v=1-f/Ess**2
    P=rr**2+a**2+2*M*a**2/rr; Pb=P+b**2/Ess**2
    return Jss*b*v/Pb+sp.sqrt(Dl*v/Pb)*sp.sqrt((Dl/rr**2)*pr**2+Jss**2/Pb)-1

def H_tau(M,a):
    """proper-time (tau) frozen Hamiltonian H_tau (paper Eq.(34), A0=1); shape verified
    against the Boyer-Lindquist F_tau to machine precision (see module docstring)."""
    f=1-2*M/rr; Dl=rr**2-2*M*rr+a**2; b=2*M*a/rr; v=1-f/Ess**2
    P=rr**2+a**2+2*M*a**2/rr; pt=Jss-b/Ess; Pb=P+b**2/Ess**2
    return pt*(b*v/Pb)+sp.sqrt(Dl*v/Pb)*sp.sqrt((Dl/rr**2)*pr**2+pt**2/Pb)-f/Ess

def build(Hsym):
    Hn=sp.lambdify((rr,pr,Ess,Jss),Hsym,'numpy')
    Hp=sp.diff(Hsym,pr); HJ=sp.diff(Hsym,Jss); G=HJ/Hp
    d=dict(Hn=Hn,
           Hp=sp.lambdify((rr,pr,Ess,Jss),Hp,'numpy'),
           Hr=sp.lambdify((rr,pr,Ess,Jss),sp.diff(Hsym,rr),'numpy'),
           HJ=sp.lambdify((rr,pr,Ess,Jss),HJ,'numpy'),
           HE=sp.lambdify((rr,pr,Ess,Jss),sp.diff(Hsym,Ess),'numpy'),
           Gpr=sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,pr),'numpy'),
           GE=sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,Ess),'numpy'),
           GJ=sp.lambdify((rr,pr,Ess,Jss),sp.diff(G,Jss),'numpy'))
    return d

def prof(D,rv,E,Jv):
    pg=np.linspace(-120,120,6001); H=D['Hn'](rv,pg,E,Jv)
    rts=[brentq(lambda p:D['Hn'](rv,p,E,Jv),pg[i],pg[i+1]) for i in range(len(pg)-1)
         if np.isfinite(H[i]) and np.isfinite(H[i+1]) and H[i]*H[i+1]<0]
    ing=[p for p in rts if D['Hp'](rv,p,E,Jv)<0]; return min(ing) if ing else np.nan

def analyse(Hbuilder,M,a,Ehat,J0,r0):
    """Return (epss, res_leading, res_exact, rc, phi0, ec, xc, flow) for one config."""
    D=build(Hbuilder(M,a))
    ev=lambda lam,y:y[1]; ev.terminal=True; ev.direction=1
    def flow(eps):
        def rhs(lam,y):
            rv,pv,ph=y; s=np.exp(-eps*lam); E=Ehat*s; Jv=J0*s
            return [D['Hp'](rv,pv,E,Jv),-D['Hr'](rv,pv,E,Jv),D['HJ'](rv,pv,E,Jv)]
        so=solve_ivp(rhs,[0,300],[r0,prof(D,r0,Ehat,J0),0.0],rtol=1e-12,atol=1e-14,
                     max_step=0.005,dense_output=True,events=ev)
        lam=np.linspace(0,so.t[-1],12000); Y=so.sol(lam); return lam,Y[0],Y[1],Y[2]
    lam0,rF,prF,phiF=flow(0.0)
    rturn=rF.min()
    phi0=interp1d(rF,phiF,bounds_error=False,fill_value='extrapolate')
    lam_r=interp1d(rF,lam0,bounds_error=False,fill_value='extrapolate')
    rg=np.linspace(r0,rturn+0.25,5000)
    # compact evaluation window: middle 45% of the arc, away from turning point
    rc=np.linspace(rturn+0.35*(r0-rturn),rturn+0.80*(r0-rturn),2000)
    prg=interp1d(rF,prF,bounds_error=False,fill_value='extrapolate')(rg)
    # leading (on-shell 1/2-Euler), from the Hamiltonian frame: dphi/dr=G, ThetaG on-shell
    ThetaG=Ehat*D['GE'](rg,prg,Ehat,J0)+J0*D['GJ'](rg,prg,Ehat,J0)
    eul=-0.5*ct(lam_r(rg)*ThetaG,rg,initial=0); ec=interp1d(rg,eul,bounds_error=False,fill_value='extrapolate')
    # exact Eq.(40)
    EulerH=Ehat*D['HE'](rg,prg,Ehat,J0)+J0*D['HJ'](rg,prg,Ehat,J0)
    S=ct(EulerH*np.gradient(lam_r(rg),rg),rg,initial=0)
    integ=D['Gpr'](rg,prg,Ehat,J0)*(lam_r(rg)*EulerH-S)/D['Hp'](rg,prg,Ehat,J0)-lam_r(rg)*ThetaG
    xc=interp1d(rg,ct(integ,rg,initial=0),bounds_error=False,fill_value='extrapolate')
    epss=np.array([0.0025,0.005,0.01,0.02,0.04])
    rh=[];rx=[]
    for e in epss:
        _,rL,_,pL=flow(e); p=interp1d(rL,pL,bounds_error=False,fill_value='extrapolate')(rc)
        rh.append(np.nanmax(np.abs(p-(phi0(rc)+e*ec(rc)))))
        rx.append(np.nanmax(np.abs(p-(phi0(rc)+e*xc(rc)))))
    return epss,np.array(rh),np.array(rx),rc,phi0,ec,xc,flow

def slope_pm(x,y):
    """log-log least-squares slope with 1sigma from the covariance."""
    lx,ly=np.log(x),np.log(y); b,cov=np.polyfit(lx,ly,1,cov=True)
    return b[0],np.sqrt(cov[0,0])

# ---- configs: t set1 (main), tau set1 (scattering), t set2 (2nd params) ----
cfgs=[("$t$, $a{=}0.9,\\hat E{=}1.4,J{=}6$",  H_eta,1.0,0.9,1.4,6.0,12.0,'C0','o'),
      ("$\\tau$, $a{=}0.9,\\hat E{=}1.4,J{=}2.5$",H_tau,1.0,0.9,1.4,2.5,12.0,'C1','s'),
      ("$t$, $a{=}0.5,\\hat E{=}1.3,J{=}5$",  H_eta,1.0,0.5,1.3,5.0,10.0,'C2','^')]
out={}
for name,Hb,M,a,E,J,r0,col,mk in cfgs:
    epss,rh,rx,rc,phi0,ec,xc,flow=analyse(Hb,M,a,E,J,r0)
    sh,she=slope_pm(epss,rh); sx,sxe=slope_pm(epss,rx)
    out[name]=(epss,rh,rx,sh,she,sx,sxe,col,mk,flow,phi0,ec,xc,rc)
    print(f"{name}: leading slope {sh:.2f}+/-{she:.2f}  exact slope {sx:.2f}+/-{sxe:.2f}")

# ---- figure ----
fig,ax=plt.subplots(1,2,figsize=(2*COL,COL*0.95))
# left: true dynamics for the main t-branch config
nm0=cfgs[0][0]; epss,rh,rx,sh,she,sx,sxe,col,mk,flow,phi0,ec,xc,rc=out[nm0]
eps=0.04; lam0,rF,prF,phiF=flow(0.0); _,rL,_,pL=flow(eps)
pt=interp1d(rL,pL,bounds_error=False,fill_value='extrapolate')
rr2=np.linspace(12.0,rF.min()+0.4,600)
ax[0].plot(rr2*np.cos(pt(rr2)),rr2*np.sin(pt(rr2)),'C0-',lw=2.6,alpha=0.35,label='true non-autonomous geodesic')
ph_h=phi0(rr2)+eps*ec(rr2); ph_x=phi0(rr2)+eps*xc(rr2)
ax[0].plot(rr2*np.cos(ph_h),rr2*np.sin(ph_h),'k--',lw=1.0,label=r'frozen $+\,\varepsilon\cdot$leading (on-shell)')
ax[0].plot(rr2*np.cos(ph_x),rr2*np.sin(ph_x),'C3:',lw=1.3,label=r'frozen $+\,\varepsilon\cdot$Eq.(40) (exact)')
ax[0].set_aspect('equal'); ax[0].set_xlabel('$x$'); ax[0].set_ylabel('$y$')
ax[0].set_title(r'true dynamics vs adiabatic ($t$-branch, $A^\prime/A=0.04$)',fontsize=6.6)
ax[0].legend(fontsize=5.4,loc='upper left',framealpha=0.9)
# right: residual convergence, both branches + 2nd set, with slope +/- 1sigma
for name,(epss,rh,rx,sh,she,sx,sxe,col,mk,flow,phi0,ec,xc,rc) in out.items():
    ax[1].loglog(epss,rx,col+mk+'-',ms=3.6,lw=1.0,
                 label=f'{name}: exact slope ${sx:.2f}\\pm{sxe:.2f}$')
    ax[1].loglog(epss,rh,col+mk+':',ms=3.0,lw=0.8,alpha=0.55)
ax[1].loglog(out[nm0][0],out[nm0][0]**2*out[nm0][2][2]/out[nm0][0][2]**2,'C7--',lw=0.6,label=r'$O(\varepsilon^2)$ guide')
ax[1].loglog(out[nm0][0],out[nm0][0]*out[nm0][1][2]/out[nm0][0][2],'C7:',lw=0.6,label=r'$O(\varepsilon)$ guide (leading, dotted)')
ax[1].set_xlabel(r"$A^\prime/A$"); ax[1].set_ylabel('residual vs true flow')
ax[1].set_title('exact term closes to $O(\\varepsilon^2)$: both branches, two parameter sets',fontsize=6.4)
ax[1].legend(fontsize=4.9,loc='upper left')
savefig(fig,os.path.join(os.path.dirname(HERE),'paper','Immagini'),'fig_phi_validation_true_dynamic')
print("FATTO")
