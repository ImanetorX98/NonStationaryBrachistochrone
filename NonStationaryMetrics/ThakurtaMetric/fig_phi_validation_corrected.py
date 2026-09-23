# -*- coding: utf-8 -*-
"""
Validation figure for the adiabatic correction against the TRUE canonical flow (main11 fix).
 Left  : true non-autonomous optical-metric geodesic vs frozen + leading (on-shell) AND vs
         frozen + the exact Eq.(40) correction (t-branch).
 Right : residual vs A'/A for BOTH branches (t and tau) and a SECOND (a,E,J) parameter set.

The true flow integrates the PHYSICAL canonical dynamics: in the normalized momentum
P_r=p_r/A the radial equation carries the dilation term -(A'/A) P_r (main11 referee). The
exact source is therefore the FULL Euler operator, S_D = int (Theta H + P_r H_Pr) dlambda;
the extra P_r H_Pr = int p_r dr (radial action) reduces in closed form to the on-shell U_k
basis + third-kind at the seed Kerr null surfaces r_pm (see tests/test_adiabatic_noreg.py and
the pass1-3 checks). Result: on-shell term slope ~1 (O(eps)); exact S_D slope ~2 (O(eps^2)),
both against the true flow. The tau-branch H_tau reproduces the Boyer-Lindquist tau shape
(kerr_adiabatic_phi_hybrid_tau.py) exactly.
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

def analyse(Hbuilder,M,a,Ehat,J0,r0,VBAR2_STOP=1e-9):
    # VBAR2_STOP is the declared tolerance at which the freezing surface counts
    # as reached.  It is a PARAMETER and not a module constant on purpose:
    # provenance/make_provenance.py imports this file by keeping only its
    # imports and function definitions, so a module-level assignment would be
    # dropped and every call from there would die on a NameError.
    """Return (epss, res_leading, res_exact, rc, phi0, ec, xc, flow) for one config."""
    D=build(Hbuilder(M,a))
    # TWO terminal events, because the rail has TWO boundaries and they are not
    # the same thing.  ev_turn is the turning point p_r = 0.  ev_freeze is the
    # FREEZING surface vbar^2 = 1 - f/E^2 = 0, where the indicatrix collapses to
    # a point; under a drift the effective energy E = Ehat exp(-eps lam) decays,
    # so a long enough run reaches it.  Without the second event the solver does
    # not stop, it FAILS there ("Required step size is less than spacing between
    # numbers") and the last sample is an integration endpoint, not a turning
    # point -- at eps = 0.04 it has p_r = -0.577, nowhere near zero.  An earlier
    # version of the left panel marked that endpoint as a turning point.
    ev=lambda lam,y:y[1]; ev.terminal=True; ev.direction=1
    def flow(eps, want_status=False):
        def rhs(lam,y):
            rv,pv,ph=y; s=np.exp(-eps*lam); E=Ehat*s; Jv=J0*s
            # TRUE canonical flow: -eps*pv is the dilation term of the normalized P_r=p_r/A
            return [D['Hp'](rv,pv,E,Jv),-D['Hr'](rv,pv,E,Jv)-eps*pv,D['HJ'](rv,pv,E,Jv)]
        # Trigger slightly ABOVE zero.  The flow stalls as vbar^2 -> 0+ (the square
        # root in the Hamiltonian degenerates) and never crosses, so an event on
        # vbar^2 itself never fires and the solver dies instead: at eps = 0.04 it
        # stops with vbar^2 = 1.1e-14 and "Required step size is less than
        # spacing between numbers".  VBAR2_STOP is the declared tolerance at
        # which we call the freezing surface reached.
        def ev_freeze(lam,y):
            E=Ehat*np.exp(-eps*lam)
            return (1.0-(1.0-2.0*M/y[0])/E**2) - VBAR2_STOP
        ev_freeze.terminal=True; ev_freeze.direction=-1
        so=solve_ivp(rhs,[0,300],[r0,prof(D,r0,Ehat,J0),0.0],rtol=1e-12,atol=1e-14,
                     max_step=0.005,dense_output=True,events=[ev,ev_freeze])
        lam=np.linspace(0,so.t[-1],12000); Y=so.sol(lam)
        if want_status:
            which=('turning' if len(so.t_events[0])
                   else 'freezing' if len(so.t_events[1]) else 'none')
            return lam,Y[0],Y[1],Y[2],which,so.success
        return lam,Y[0],Y[1],Y[2]
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
    # exact Eq.(40) with the CORRECTED terminally-anchored source
    # S_D = int (Theta H + P_r H_Pr) dlam  (the extra P_r*H_Pr = int p_r dr is the radial-momentum
    # dilation of the time-dependent transform P_r=p_r/A; omitting it gives a spurious O(eps) error).
    EulerH=Ehat*D['HE'](rg,prg,Ehat,J0)+J0*D['HJ'](rg,prg,Ehat,J0)
    PrHpr=prg*D['Hp'](rg,prg,Ehat,J0)
    S=ct((EulerH+PrHpr)*np.gradient(lam_r(rg),rg),rg,initial=0)        # S_D
    integ=D['Gpr'](rg,prg,Ehat,J0)*(lam_r(rg)*EulerH-S)/D['Hp'](rg,prg,Ehat,J0)-lam_r(rg)*ThetaG
    xc=interp1d(rg,ct(integ,rg,initial=0),bounds_error=False,fill_value='extrapolate')
    epss=np.array([0.001,0.002,0.004,0.008,0.016])
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

# Dump the raw (eps, residual) pairs so the figure, the table and any refit come
# from ONE dataset (CQG-116884, referee major 7).  paper2/provenance reads this
# file; nothing downstream re-runs the physics or re-types a number.
_raw = {}
for i, (name, *_rest) in enumerate(cfgs):
    e, rh_, rx_, sh_, she_, sx_, sxe_ = out[name][:7]
    _raw[f"eps_{i}"] = e
    _raw[f"res_leading_{i}"] = rh_
    _raw[f"res_exact_{i}"] = rx_
    _raw[f"label_{i}"] = np.array(name)
# make_provenance.py --check re-runs this generator only to compare numbers; it
# must not touch the working tree.  savefig is already guarded in paper_style;
# this raw dataset was still written unconditionally, so the "writes nothing"
# claim was false for the full check.
if os.environ.get('PROVENANCE_NO_WRITE'):
    print(f"  [PROVENANCE_NO_WRITE] skipped adiabatic_convergence_raw.npz "
          f"({len(cfgs)} configs x {len(out[cfgs[0][0]][0])} epsilon values)")
else:
    np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'adiabatic_convergence_raw.npz'), **_raw)
    print(f"  raw convergence data -> adiabatic_convergence_raw.npz "
          f"({len(cfgs)} configs x {len(out[cfgs[0][0]][0])} epsilon values)")

# ---- figure ----
fig,ax=plt.subplots(1,2,figsize=(2*COL,COL*0.95))
# left: true dynamics for the main t-branch config
nm0=cfgs[0][0]; epss,rh,rx,sh,she,sx,sxe,col,mk,flow,phi0,ec,xc,rc=out[nm0]
eps=0.04
lam0,rF,prF,phiF,whichF,okF=flow(0.0, want_status=True)
_,rL,_,pL,whichL,okL=flow(eps, want_status=True)
print(f"  frozen  run: terminated by {whichF:8s} (solver success={okF})")
print(f"  dynamic run: terminated by {whichL:8s} (solver success={okL})")
# The dynamic trajectory turns at a LARGER radius than the frozen one, so a grid
# taken from the frozen r_min pushed 221 of 600 points below rL.min() and the
# curve drawn as "true non-autonomous" was the interpolant extrapolated past the
# integrated domain.  Restrict to the window common to every curve shown, and
# forbid extrapolation so the mistake cannot recur silently.
pt=interp1d(rL,pL,bounds_error=True)
r_lo=max(float(rF.min()),float(rL.min()))+0.05
rr2=np.linspace(12.0,r_lo,600)
assert rr2.min()>=float(rL.min()) and rr2.max()<=float(rL.max()), \
    "true-dynamics panel would extrapolate"
print(f"  true-dynamics panel: r in [{rr2.min():.4f}, {rr2.max():.4f}]; "
      f"frozen r_min={rF.min():.4f}, dynamic r_min={rL.min():.4f} "
      f"(no extrapolation)")
ax[0].plot(rr2*np.cos(pt(rr2)),rr2*np.sin(pt(rr2)),'C0-',lw=2.6,alpha=0.35,label='true non-autonomous geodesic')
# Mark where each run ENDS, and say which boundary ended it -- they are not the
# same boundary.  The frozen run turns (p_r = 0).  The drifting run does not: its
# effective energy decays as E = Ehat exp(-eps lam), and at eps = 0.04 it reaches
# the FREEZING surface vbar^2 = 0 first, with p_r still around -0.58.  Labelling
# that endpoint a turning point, as an earlier version did, misnames the physics
# and misreads a failed integration as a turn.
_lblL = (r'dynamic: turning point ($p_r=0$)' if whichL == 'turning'
         else r'dynamic: freezing surface ($\bar v^2=0$)' if whichL == 'freezing'
         else 'dynamic: integration endpoint (no boundary reached)')
_rt = float(rL.min()); _pt_t = float(pt(_rt))
ax[0].plot([_rt*np.cos(_pt_t)],[_rt*np.sin(_pt_t)],'o',ms=5.0,mfc='none',
           mec='C1',mew=1.4,label=_lblL)
_rtF = float(rF.min())
ax[0].plot([_rtF*np.cos(float(phi0(_rtF)))],[_rtF*np.sin(float(phi0(_rtF)))],
           's',ms=4.0,mfc='none',mec='0.45',mew=1.0,
           label=r'frozen: turning point ($p_r=0$)')
ph_h=phi0(rr2)+eps*ec(rr2); ph_x=phi0(rr2)+eps*xc(rr2)
ax[0].plot(rr2*np.cos(ph_h),rr2*np.sin(ph_h),'k--',lw=1.0,label=r'frozen $+\,\varepsilon\cdot$leading (on-shell)')
ax[0].plot(rr2*np.cos(ph_x),rr2*np.sin(ph_x),'C3:',lw=1.3,label=r'frozen $+\,\varepsilon\cdot$ complete first-order (exact)')
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
