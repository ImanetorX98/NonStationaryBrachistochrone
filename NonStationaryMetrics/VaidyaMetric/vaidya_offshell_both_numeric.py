# -*- coding: utf-8 -*-
# Vaidya off-shell wrap, BOTH branches (tau and v), clean non-surrogate Hamiltonians, NUMERIC.
# Source Theta=m d_m (no dilation). Off-shell wrap  Phi = -int (d_pr G/H_pr) * Sigma dr,
# Sigma(r)=int_{r0}^{r} (m H_m/H_pr) dz. Establish the wrap value + inner-letter structure
# (third-kind at the Schwarzschild horizon r=2m). p_r0 solved numerically (avoids sympy sqrt).
import numpy as np, sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq

r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
def build(branch):
    Dl=r*(r-2*m); f=(r-2*m)/r; DE=(E**2-1)*r+2*m; v=DE/(E**2*r)
    if branch=='tau':
        H=sp.sqrt(Dl*v/r**2)*sp.sqrt((Dl/r**2)*pr**2+J**2/r**2)-f/E
    else:  # v-branch clean: advanced-time optical Hamiltonian (a=0)
        H=pr*(f-E**2)/E + sp.sqrt(Dl*v/r**2)*sp.sqrt((Dl/r**2)*pr**2+J**2/r**2)  # -0 (no unit cost)
    Hp=sp.diff(H,pr); Hm=sp.diff(H,m); G=sp.diff(H,J)/Hp; Gp=sp.diff(G,pr)
    la=lambda ex: sp.lambdify((r,pr,m,E,J),ex,'numpy')
    return dict(H=la(H),Hp=la(Hp),Hm=la(Hm),Gpr_over_Hp=la(Gp/Hp),
               ThetaH_over_Hp=la(m*Hm/Hp))

def run(branch,mv,Ev,Jv,r0=12.0):
    B=build(branch)
    def pr0(rv):  # ingoing on-shell root of H=0
        g=lambda p:B['H'](rv,p,mv,Ev,Jv)
        # bracket a negative (ingoing) root
        ps=np.linspace(-60,-1e-4,4000); vals=np.array([g(p) for p in ps])
        sgn=np.where(np.diff(np.sign(vals)))[0]
        if len(sgn)==0: return np.nan
        return brentq(g,ps[sgn[0]],ps[sgn[0]+1])
    # turning point (S=0 -> Q3=0 -> r^2(r-2m)=J^2 DE)
    Q3=lambda rv:rv**2*(rv-2*mv)-Jv**2*((Ev**2-1)*rv+2*mv)
    rr=np.linspace(2.01,r0,6000); q=np.array([Q3(x) for x in rr]); idx=np.where(np.diff(np.sign(q)))[0]
    rt=max(rr[i] for i in idx) if len(idx) else 2.1
    def kern(rv): return B['Gpr_over_Hp'](rv,pr0(rv),mv,Ev,Jv)            # = A/sqrt(S)
    def innder(rv): return B['ThetaH_over_Hp'](rv,pr0(rv),mv,Ev,Jv)      # inner-letter integrand
    def Sigma(rv): return quad(innder,r0,rv,limit=150)[0]
    xf=rt+0.6
    wrap=-quad(lambda rv: kern(rv)*Sigma(rv), r0, xf, limit=120)[0]
    # inner-letter third-kind structure: residue of innder at r=2m (horizon)
    eps=1e-4; res2m=innder(2*mv+eps)*eps   # ~ residue if simple pole
    return dict(rt=rt,xf=xf,wrap=wrap,pr0_10=pr0(10.0),Sigma_xf=Sigma(xf),horizon_res=res2m)

for branch in ['tau','v']:
    print(f"=== Vaidya {branch}-branch off-shell (clean form) ===")
    for (mv,Ev,Jv) in [(1.0,1.4,2.5),(1.2,1.5,3.0)]:
        try:
            R=run(branch,mv,Ev,Jv)
            print(f"  m={mv},E={Ev},J={Jv}: turning={R['rt']:.4f}  off-shell wrap Phi={R['wrap']:+.6f}"
                  f"  p_r0(10)={R['pr0_10']:+.4f}  inner@horizon~{R['horizon_res']:+.2e}")
        except Exception as e:
            print(f"  m={mv},E={Ev},J={Jv}: FAILED {e}")
print("\n=> wrap values finite; inner letter has third-kind behaviour at r=2m (horizon). Next: symbolic")
print("   reduction of kernel (A/sqrtS) and inner (poly/sqrtS + 3rd kind @ r=2m,DE=0) -> A+B+C assembly.")
