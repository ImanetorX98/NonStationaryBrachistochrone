import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct
from scipy.optimize import brentq
from scipy.interpolate import interp1d
M,a,Ehat=1.0,0.9,1.4; J0,r0=6.0,12.0
r,Es,Js=sp.symbols('r E J_'); Dl=r**2-2*M*r+a**2
Q2=(2*Es**2*Js**2*M*r-Es**2*Js**2*r**2-4*Es**2*Js*M*a*r+2*Es**2*M*a**2*r+Es**2*a**2*r**2
    +Es**2*r**4+4*Js**2*M**2-4*Js**2*M*r+Js**2*r**2-8*Js*M**2*a+4*Js*M*a*r+4*M**2*a**2)
R=r*Q2*((Es**2-1)*r+2*M); Kf=r*((Es**2-1)*r+2*M)*(Js*(r-2*M)+2*M*a)/Dl; F=Kf/sp.sqrt(R)
Fn=sp.lambdify((r,Es,Js),F,'numpy'); dEF=sp.lambdify((r,Es,Js),sp.diff(F,Es),'numpy')
dJF=sp.lambdify((r,Es,Js),sp.diff(F,Js),'numpy')
rr,pr,Ess,Jss=sp.symbols('r pr E J_')
f2=1-2*M/rr; Dl2=rr**2-2*M*rr+a**2; b2=2*M*a/rr; v2=1-f2/Ess**2
P2=rr**2+a**2+2*M*a**2/rr; Pb2=P2+b2**2/Ess**2
H2=Jss*b2*v2/Pb2+sp.sqrt(Dl2*v2/Pb2)*sp.sqrt((Dl2/rr**2)*pr**2+Jss**2/Pb2)-1
H2n=sp.lambdify((rr,pr,Ess,Jss),H2,'numpy')
dHp=sp.lambdify((rr,pr,Ess,Jss),sp.diff(H2,pr),'numpy')
dHr=sp.lambdify((rr,pr,Ess,Jss),sp.diff(H2,rr),'numpy')
dHJ=sp.lambdify((rr,pr,Ess,Jss),sp.diff(H2,Jss),'numpy')
dHE=sp.lambdify((rr,pr,Ess,Jss),sp.diff(H2,Ess),'numpy')
def prof(rv,E,Jv):
    pg=np.linspace(-80,80,3001); Hv=H2n(rv,pg,E,Jv)
    rts=[brentq(lambda p:H2n(rv,p,E,Jv),pg[i],pg[i+1]) for i in range(len(pg)-1)
         if np.isfinite(Hv[i]) and np.isfinite(Hv[i+1]) and Hv[i]*Hv[i+1]<0]
    ing=[p for p in rts if dHp(rv,p,E,Jv)<0]; return min(ing) if ing else np.nan
ev=lambda lam,y:y[1]; ev.terminal=True; ev.direction=1
def flow(eps):
    def rhs(lam,y):
        rv,pv,ph,t=y; s=np.exp(-eps*lam); E=Ehat*s; Jv=J0*s
        return [dHp(rv,pv,E,Jv),-dHr(rv,pv,E,Jv),dHJ(rv,pv,E,Jv),abs(dHE(rv,pv,E,Jv))]
    so=solve_ivp(rhs,[0,300],[r0,prof(r0,Ehat,J0),0.0,0.0],rtol=1e-12,atol=1e-14,
                 max_step=0.005,dense_output=True,events=ev)
    lam=np.linspace(0,so.t[-1],9000); Y=so.sol(lam); return lam,Y[0],Y[2],Y[3]
lam0,rF,phiF,tF=flow(0.0)
phi0f=interp1d(rF,phiF,bounds_error=False,fill_value='extrapolate')
lam_of_r=interp1d(rF,lam0,bounds_error=False,fill_value='extrapolate')
t_of_r=interp1d(rF,tF,bounds_error=False,fill_value='extrapolate')
rg=np.linspace(r0,7.0,4000); rc=np.linspace(8.0,11.0,2000); s_=-1.0; phi0_shape=lambda x:-phi0f(x)
EJ=lambda x: Ehat*dEF(x,Ehat,J0)+J0*dJF(x,Ehat,J0)
for lbl,clock in [("affine λ",lam_of_r),("coord t",t_of_r)]:
    eta=clock(rg); euler=-0.5*ct(EJ(rg)*eta,rg,initial=0); ec=interp1d(rg,euler,bounds_error=False,fill_value='extrapolate')
    epss=np.array([0.0025,0.005,0.01,0.02]); res=[]
    for eps in epss:
        _,rL,pL,_=flow(eps); pt=interp1d(rL,pL,bounds_error=False,fill_value='extrapolate')(rc)
        res.append(np.nanmax(np.abs(pt-s_*(phi0_shape(rc)+eps*ec(rc)))))
    res=np.array(res); sl=np.polyfit(np.log(epss),np.log(res),1)[0]
    print(f"  1/2 Euler + {lbl}: coeff={np.mean(ec(rc)):+.4f}  slope={sl:.2f}  res@0.0025={res[0]:.2e}")

print("\n=== OFF-SHELL costate p_eta (referee 4.6 term): dp_eta/dlam = +(E dH/dE + J dH/dJ) ===")
# along frozen orbit compute Euler_H = E*dHE + J*dHJ, integrate -> p_eta/eps
# need pr along frozen orbit
_,rFa,phiFa,tFa=flow(0.0)
# re-run frozen to get pr
def flow_pr(eps):
    def rhs(lam,y):
        rv,pv,ph,t=y; s=np.exp(-eps*lam); E=Ehat*s; Jv=J0*s
        return [dHp(rv,pv,E,Jv),-dHr(rv,pv,E,Jv),dHJ(rv,pv,E,Jv),abs(dHE(rv,pv,E,Jv))]
    so=solve_ivp(rhs,[0,300],[r0,prof(r0,Ehat,J0),0.0,0.0],rtol=1e-12,atol=1e-14,max_step=0.005,dense_output=True,events=ev)
    lam=np.linspace(0,so.t[-1],9000); Y=so.sol(lam); return lam,Y[0],Y[1]
lam0,rF2,prF=flow_pr(0.0)
EulerH=Ehat*dHE(rF2,prF,Ehat,J0)+J0*dHJ(rF2,prF,Ehat,J0)
p_eta_over_eps=ct(EulerH,lam0,initial=0)   # p_eta/eps (H_ext=0 => H2=-p_eta ~ O(eps))
pe=interp1d(rF2,p_eta_over_eps,bounds_error=False,fill_value='extrapolate')
print(f"  p_eta/eps on [8,11] (=off-shell H2 excursion/eps): range[{pe(rc).min():+.4f},{pe(rc).max():+.4f}]")
print(f"  => at eps=0.0025 the orbit is off-shell H2 by ~{0.0025*np.abs(pe(rc)).max():.2e}")
print(f"  residual of 1/2Euler (affine) at eps=0.0025 was 8.4e-5 -> same O(eps) scale as off-shell term")
print(f"  leading correction coeff ~0.336; off-shell |p_eta/eps|~{np.abs(pe(rc)).mean():.3f} -> ~{100*np.abs(pe(rc)).mean()/0.336:.0f}% of leading (matches ~10% residual slope 1.25)")
