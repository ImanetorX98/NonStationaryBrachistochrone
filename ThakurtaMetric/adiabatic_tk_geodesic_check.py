import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, cumulative_trapezoid as ct
from scipy.interpolate import interp1d
M,a=1.0,0.9
# equatorial Kerr inverse metric (theta=pi/2)
r=sp.symbols('r',positive=True); Dl=r**2-2*M*r+a**2
gtt=-(r**2+a**2+2*M*a**2/r)/Dl; gtp=-2*M*a/(r*Dl); gpp=(1-2*M/r)/Dl; grr=Dl/r**2
gtt_n,gtp_n,gpp_n,grr_n=[sp.lambdify(r,g,'numpy') for g in (gtt,gtp,gpp,grr)]
dgtt,dgtp,dgpp,dgrr=[sp.lambdify(r,sp.diff(g,r),'numpy') for g in (gtt,gtp,gpp,grr)]
def Af(t,eps): return np.exp(eps*t)
def bracket(rv,pt,J,prr): return gtt_n(rv)*pt*pt+2*gtp_n(rv)*pt*J+gpp_n(rv)*J*J+grr_n(rv)*prr*prr
def dbracket_dr(rv,pt,J,prr): return dgtt(rv)*pt*pt+2*dgtp(rv)*pt*J+dgpp(rv)*J*J+dgrr(rv)*prr*prr
E0,J,r0=1.4,6.0,12.0
def geod(eps):
    A0=1.0; pt0=-E0
    # initial p_r from shell H=-1/2 : A^-2*bracket=-1 => bracket=-A^2 => grr*pr^2 = -A^2 - (tphi part)
    A=A0
    rest=gtt_n(r0)*pt0*pt0+2*gtp_n(r0)*pt0*J+gpp_n(r0)*J*J
    pr2=(-A*A-rest)/grr_n(r0)
    pr0=-np.sqrt(pr2)  # ingoing
    ev=lambda s,y:y[3]; ev.terminal=True; ev.direction=1  # turning pr=0
    def rhs(s,y):
        t,rv,ph,prr,pt=y; A=Af(t,eps); Ap=eps*A
        dt=A**-2*(gtt_n(rv)*pt+gtp_n(rv)*J)
        dr=A**-2*grr_n(rv)*prr
        dph=A**-2*(gtp_n(rv)*pt+gpp_n(rv)*J)
        dpr=-0.5*A**-2*dbracket_dr(rv,pt,J,prr)
        dpt=A**-3*Ap*bracket(rv,pt,J,prr)
        return [dt,dr,dph,dpr,dpt]
    s=solve_ivp(rhs,[0,400],[0.0,r0,0.0,pr0,pt0],rtol=1e-12,atol=1e-14,max_step=0.01,
                dense_output=True,events=ev)
    lam=np.linspace(0,s.t[-1],9000); Y=s.sol(lam)
    return dict(s=lam,t=Y[0],r=Y[1],phi=Y[2],pr=Y[3],pt=Y[4])
# sanity: frozen
g0=geod(0.0)
Hs=0.5*Af(g0['t'],0)**-2*bracket(g0['r'],g0['pt'],J,g0['pr'])
print(f"SANITY frozen: H in [{Hs.min():.2e},{Hs.max():.2e}] (should be -0.5); E=-pt const? pt range[{g0['pt'].min():.4f},{g0['pt'].max():.4f}]")
# sanity: dE/ds = A'/A for slow eps
eps=0.01; ge=geod(eps)
dEds=np.gradient(-ge['pt'],ge['s']); ApA=eps*np.ones_like(ge['s'])
print(f"SANITY dE/ds vs A'/A (eps={eps}): mean dE/ds={np.mean(dEds):.4f} vs A'/A={eps} ; E drift over orbit: {-ge['pt'][-1]+ge['pt'][0]:.4f}")
print(f"  J conserved (p_phi fixed by construction). turning r: frozen={g0['r'].min():.3f}")

print("\n=== STRUCTURE: what drifts along the geodesic? (eps=0.01) ===")
ge=geod(0.01)
A_=Af(ge['t'],0.01)
E_eff=-ge['pt']/A_       # effective Kerr energy (unit mass)
J_eff=J/A_               # effective Kerr angular momentum
# restrict to r in [8,11]
m=(ge['r']<11)&(ge['r']>8)
print(f"  along geodesic on r in [8,11]:")
print(f"   E_eff=E_BL/A : range[{E_eff[m].min():.4f},{E_eff[m].max():.4f}]  drift={E_eff[m].max()-E_eff[m].min():.4f}")
print(f"   J_eff=J/A    : range[{J_eff[m].min():.4f},{J_eff[m].max():.4f}]  drift={J_eff[m].max()-J_eff[m].min():.4f}")
print(f"   E_BL=-pt     : range[{(-ge['pt'])[m].min():.4f},{(-ge['pt'])[m].max():.4f}]")
print(f"   A            : range[{A_[m].min():.4f},{A_[m].max():.4f}]")
print(f"  => tells us which parameter(s) the frozen orbit must be differentiated in")

print("\n=== STRUCTURAL CHECK: is dphi/dr = F(r; E_eff, J_eff) EXACT along geodesic? ===")
import sympy as sp2
rs,Es2,Js2=sp2.symbols('r E J_',positive=True); Dl2=rs**2-2*M*rs+a**2
gtp2=-2*M*a/(rs*Dl2); gpp2=(1-2*M/rs)/Dl2; grr2=Dl2/rs**2
gtt2=-(rs**2+a**2+2*M*a**2/rs)/Dl2
# unit-mass shell: gtt E^2? careful: p_t=-E_eff*? For unit mass Kerr: g^ab p_a p_b=-1
# p_t=-E_eff, p_phi=J_eff, solve g^rr pr^2 = -1 - (gtt E^2 -2gtp E J + gpp J^2)
pr_unit=sp2.sqrt((-1-(gtt2*Es2**2-2*gtp2*Es2*Js2+gpp2*Js2**2))/grr2)
F_check=(-gtp2*Es2+gpp2*Js2)/(grr2*(-pr_unit))  # dphi/dr = (dphi/ds)/(dr/ds), ingoing pr<0
Fc=sp2.lambdify((rs,Es2,Js2),F_check,'numpy')
ge=geod(0.01); A_=Af(ge['t'],0.01); Ee=-ge['pt']/A_; Je=J/A_
dphidr_num=np.gradient(ge['phi'],ge['r'])
m=(ge['r']<11)&(ge['r']>8.5)
Fpred=Fc(ge['r'][m],Ee[m],Je[m])
print(f"  max|dphi/dr_numeric - F(r;E_eff,J_eff)| on [8.5,11] = {np.nanmax(np.abs(dphidr_num[m]-Fpred)):.2e}")
print(f"  (if ~0 => geodesic shape IS exactly frozen-family with instantaneous charges)")
print(f"  sample: dphi/dr_num={dphidr_num[m][100]:.5f}  F(E_eff,J_eff)={Fpred[100]:.5f}")
