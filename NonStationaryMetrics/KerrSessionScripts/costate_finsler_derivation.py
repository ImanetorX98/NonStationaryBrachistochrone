# -*- coding: utf-8 -*-
# Referee Issue 2 (BLOCKER) -- RIGOROUS Legendre map for the rail brachistochrone.
# Min-time Finsler function F(r,r',phi') = dt/dsigma for the controlled worldline x'^mu=(t',r',phi')
# with rail u_t = g_{t mu}x'^mu / N = -E,  N=sqrt(-g_{ab}x'^a x'^b).  Solve the rail quadratic for
# t'=F; the Pontryagin costate is p_phi=dF/dphi'.  Compare to the MECHANICAL angular momentum
# u_phi = g_{phi mu}u^mu = (g_{t phi}t'+g_{phi phi}phi')/N.  Then test conformal weights.
import numpy as np

M,a=1.0,0.9
def metric(r,A=1.0):
    Dl=r**2-2*M*r+a**2
    g_tt=-(1-2*M/r); g_tp=-2*M*a/r; g_rr=r**2/Dl; g_pp=r**2+a**2+2*M*a**2/r
    return A**2*np.array([[g_tt,0,g_tp],[0,g_rr,0],[g_tp,0,g_pp]])   # (t,r,phi), conformal A^2

def solve_tprime(r,rp,phip,E,A=1.0):
    # rail: g_{t mu}x'^mu = -E * N,  N^2=-(g_ab x'^a x'^b). x'=(t',rp,phip). Solve quadratic for t'.
    g=metric(r,A); gtt,gtp,gpp,grr=g[0,0],g[0,2],g[2,2],g[1,1]
    # let a2 t'^2 + a1 t' + a0 relation from (g_tt t' + g_tp phip)^2 = E^2 * (-(gtt t'^2+2gtp t' phip+gpp phip^2+grr rp^2))
    L=gtt; Kc=gtp*phip                     # u_t*N = L t' + Kc
    # N^2 = -(gtt t'^2 + 2 gtp t' phip + gpp phip^2 + grr rp^2)
    # (L t'+Kc)^2 = E^2 * N^2
    A2=L**2 + E**2*gtt
    A1=2*L*Kc + E**2*2*gtp*phip
    A0=Kc**2 + E**2*(gpp*phip**2+grr*rp**2)
    disc=A1**2-4*A2*A0
    if disc<0: return None
    roots=[(-A1+np.sqrt(disc))/(2*A2),(-A1-np.sqrt(disc))/(2*A2)]
    # future-directed: t'>0
    tp=[t for t in roots if t>0]
    return max(tp) if tp else None

def F(r,rp,phip,E,A=1.0): return solve_tprime(r,rp,phip,E,A)

def costate_pphi(r,rp,phip,E,A=1.0,h=1e-6):
    fp,fm=F(r,rp,phip+h,E,A),F(r,rp,phip-h,E,A)
    if fp is None or fm is None: return None
    return (fp-fm)/(2*h)   # dF/dphi'
def costate_pr(r,rp,phip,E,A=1.0,h=1e-6):
    fp,fm=F(r,rp+h,phip,E,A),F(r,rp-h,phip,E,A)
    if fp is None or fm is None: return None
    return (fp-fm)/(2*h)

def mech_uphi(r,rp,phip,E,A=1.0):
    g=metric(r,A); tp=F(r,rp,phip,E,A)
    if tp is None: return None
    xup=np.array([tp,rp,phip]); N=np.sqrt(-xup@g@xup)
    u=xup/N                                    # u^mu
    return (g@u)[2]                            # u_phi = g_{phi mu}u^mu

def check_ut(r,rp,phip,E,A=1.0):
    g=metric(r,A); tp=F(r,rp,phip,E,A); xup=np.array([tp,rp,phip]); N=np.sqrt(-xup@g@xup)
    u=xup/N; return (g@u)[0]                    # should be -E

# ---- test at a representative state/velocity ----
r,rp,phip,E=6.0,-1.0,0.2,1.4
print(f"state r={r}, velocity (r',phi')=({rp},{phip}), rail E={E}")
Fv=F(r,rp,phip,E); pphi=costate_pphi(r,rp,phip,E); uphi=mech_uphi(r,rp,phip,E)
print(f"  F=dt/dsigma      = {Fv:+.6f}")
print(f"  u_t (rail check) = {check_ut(r,rp,phip,E):+.6f}  (=-E)")
print(f"  costate p_phi=dF/dphi' = {pphi:+.6f}")
print(f"  mechanical u_phi       = {uphi:+.6f}")
print(f"  ratio p_phi/u_phi      = {pphi/uphi:+.6f}")
print(f"  => costate {'==' if abs(pphi-uphi)<1e-4 else '!='} mechanical")

# ---- Euler homogeneity check: F = r' p_r + phi' p_phi (F degree-1 homogeneous) ----
prr=costate_pr(r,rp,phip,E); euler=rp*prr+phip*pphi
print(f"\nEuler check: r' p_r + phi' p_phi = {euler:+.6f}  vs F={Fv:+.6f}  (equal => F degree-1 homog OK)")

# ---- DECISIVE conformal-weight test: min-COORDINATE-TIME is conformally invariant ----
# Physical problem (metric A^2 g, rail Ehat) vs frozen-Kerr description (metric g, rail Ehat/A).
# If F_phys(Ehat,A) == F_frozen(Ehat/A,1) then the costate p_phi is conformally INVARIANT (weight 0);
# the mechanical u_phi carries the A-weight. This settles which object the frozen J_eff=J/A is.
print("\nDECISIVE test  F(Ehat,A) vs F(Ehat/A, 1)  [min-coordinate-time conformal invariance]:")
Ehat=1.4
for A in [1.5,2.0]:
    fa=F(r,rp,phip,Ehat,A); fb=F(r,rp,phip,Ehat/A,1.0)
    pa=costate_pphi(r,rp,phip,Ehat,A); pb=costate_pphi(r,rp,phip,Ehat/A,1.0)
    ua=mech_uphi(r,rp,phip,Ehat,A); ub=mech_uphi(r,rp,phip,Ehat/A,1.0)
    if None in (fa,fb,pa,pb,ua,ub):
        print(f"  A={A}: velocity not admissible under A^2 g (skip)"); continue
    print(f"  A={A}: F_phys={fa:+.5f} F_frozen(E/A)={fb:+.5f} diff={abs(fa-fb):.1e}")
    print(f"        p_phi: phys={pa:+.5f} frozen(E/A)={pb:+.5f}  (equal=>costate weight 0)")
    print(f"        u_phi: phys={ua:+.5f} frozen(E/A)={ub:+.5f}  ratio phys/frozen={ua/ub:.4f} (=A =>u_phi weight -1 frozen)")
print("\nInterpretation: if costate p_phi is conformally invariant but the paper's J_eff scales as 1/A,")
print("then the paper's 'J' is the MECHANICAL angular momentum, NOT the conserved Pontryagin costate.")

# ---- FINAL: is the paper's (p_r0,J) the Finsler costate or the mechanical geodesic momentum? ----
# paper on-shell p_r0 = sqrt(S)/(Delta*Emu),  S=r*Emu*Q2.
import numpy as _np
def paper_pr0(rv,Ev,Jv):
    M_,a_=1.0,0.9; Emu=(Ev**2-1)*rv+2*M_; Dl=rv**2-2*M_*rv+a_**2
    Q2=(2*Ev**2*Jv**2*M_*rv-Ev**2*Jv**2*rv**2-4*Ev**2*Jv*M_*a_*rv+2*Ev**2*M_*a_**2*rv+Ev**2*a_**2*rv**2
        +Ev**2*rv**4+4*Jv**2*M_**2-4*Jv**2*M_*rv+Jv**2*rv**2-8*Jv*M_**2*a_+4*Jv*M_*a_*rv+4*M_**2*a_**2)
    S=rv*Emu*Q2; return _np.sqrt(S)/(Dl*Emu)
def geodesic_ur(rv,Ev,Jv):   # covariant u_r for a Kerr geodesic with charges (E,J): mass-shell
    M_,a_=1.0,0.9; Dl=rv**2-2*M_*rv+a_**2
    g_tt=-(1-2*M_/rv); g_tp=-2*M_*a_/rv; g_pp=rv**2+a_**2+2*M_*a_**2/rv
    Det=g_tt*g_pp-g_tp**2; gtt=g_pp/Det; gpp=g_tt/Det; gtp=-g_tp/Det; grr=Dl/rv**2
    pr2=(-1-(gtt*Ev**2-2*gtp*Ev*Jv+gpp*Jv**2))/grr; return _np.sqrt(pr2) if pr2>0 else float('nan')
# Finsler figuratrix: vary velocity direction, find where p_phi=J, read p_r
def finsler_pr_at(rv,Ev,Jtarget):
    best=None
    for phip in _np.linspace(-2.0,2.0,20001):
        pp=costate_pphi(rv,-1.0,phip,Ev)
        if pp is None: continue
        if best is None or abs(pp-Jtarget)<best[0]:
            best=(abs(pp-Jtarget),costate_pr(rv,-1.0,phip,Ev),phip)
    return best
rv,Ev,Jv=6.0,1.4,6.0
p0=paper_pr0(rv,Ev,Jv); gur=geodesic_ur(rv,Ev,Jv); fb=finsler_pr_at(rv,Ev,Jv)
print("\n=== FINAL: identify the paper's J ===")
print(f"  paper p_r0(E,J)          = {p0:+.6f}")
print(f"  geodesic (mechanical) u_r = {gur:+.6f}  (if paper==this => J is mechanical u_phi, weight -1)")
print(f"  Finsler costate p_r @ p_phi=J = {fb[1]:+.6f} (match |p_phi-J|={fb[0]:.1e}) (if paper==this => J is costate, weight 0)")
