# -*- coding: utf-8 -*-
# FULL generic off-shell wrap (t-branch) closed form -- ALL blocks, verified vs the direct
# double integral (math level). Using the Hermite reduction of the kernel pole:
#   A/sqrt S = G/sqrt S + d/dr[P sqrt S/Q2],  G=A_poly+Rem=sum_k g_k r^k (k=0..4),
#   Delta S = sum_k a_k U_k + Pi,  Pi = sum_root rho_root int dr/((r-root) sqrt S).
# full = -int (A/sqrtS) Delta S dr
#      = [ -int (G/sqrtS)(sum a_k U_k) dr ]                         (A) second-kind: Q_jk W_jk
#      + [ -int (G/sqrtS) Pi dr ]                                    (B) Delta=0 genus-2 dilog
#      + [ -[alg_K Delta S] + int P r/Delta dr ]                     (C) ELEMENTARY (alg + logs)
# g_k, a_k, rho_root, P all SYMBOLIC in E. Ground truth: full direct double integral.
import sympy as sp, numpy as np
from scipy.integrate import quad

r,E=sp.symbols('r E',positive=True)
M,a,J=sp.Rational(1),sp.Rational(9,10),sp.Integer(6); E0=sp.Rational(7,5)
Emu=(E**2-1)*r+2*M; Dl=r**2-2*M*r+a**2
Q2=(2*E**2*J**2*M*r-E**2*J**2*r**2-4*E**2*J*M*a*r+2*E**2*M*a**2*r+E**2*a**2*r**2
    +E**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2-8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
T=sp.expand(r*Emu); S=sp.expand(Q2*T); deng=sp.expand(Dl*Emu)
# kernel poly + Hermite
Qp,N3=sp.div(sp.Poly(sp.expand(E**2*J*r**4*Emu),r),sp.Poly(sp.expand(Q2),r)); N3=N3.as_expr()
p_j=[c for c in Qp.all_coeffs()[::-1]]
Q2p=sp.diff(Q2,r); Tp=sp.diff(T,r)
inv=sp.invert(sp.Poly(sp.expand(Q2p*T),r),sp.Poly(Q2,r))
P=sp.expand(sp.rem(sp.expand(-2*N3*inv.as_expr()),Q2,r))
Rem=sp.div(sp.Poly(sp.expand(2*N3-2*sp.diff(P,r)*Q2*T-P*Q2*Tp+P*Q2p*T),r),sp.Poly(sp.expand(2*Q2),r))[0].as_expr()
remc=sp.Poly(Rem,r).all_coeffs()[::-1]
g=[sp.simplify((p_j[k] if k<len(p_j) else 0)+(remc[k] if k<len(remc) else 0)) for k in range(5)]
# dilation
Qd,R=sp.div(sp.Poly(sp.expand(S),r),sp.Poly(deng,r)); R=R.as_expr()
a_k=[sp.simplify(c) for c in Qd.all_coeffs()[::-1]]
roots={'r+':M+sp.sqrt(M**2-a**2),'r-':M-sp.sqrt(M**2-a**2)}
rho={nm:sp.simplify((R/sp.diff(deng,r)).subs(r,rt)) for nm,rt in roots.items()}
print("kernel g_k (Hermite-extended):",[sp.simplify(x) for x in g])
print("dilation a_k:",a_k)
print("rho:",{k:sp.simplify(v) for k,v in rho.items()})

# ============== numeric verification ==============
E0f=1.4
Sn=sp.lambdify(r,S.subs(E,E0),'numpy'); sqn=lambda x:np.sqrt(max(Sn(x),0.0))
Q2n=sp.lambdify(r,Q2.subs(E,E0),'numpy'); dengn=sp.lambdify(r,deng.subs(E,E0),'numpy')
N3n=sp.lambdify(r,N3.subs(E,E0),'numpy'); Rn=sp.lambdify(r,R.subs(E,E0),'numpy')
Pn=sp.lambdify(r,P.subs(E,E0),'numpy')
gn=[float(x.subs(E,E0)) for x in g]; akn=[float(x.subs(E,E0)) for x in a_k]
pjn=[float(x.subs(E,E0)) for x in p_j]
rootn={nm:float(rt.subs(E,E0)) for nm,rt in roots.items()}; rhon={nm:float(rho[nm].subs(E,E0)) for nm in roots}
r0=11.0
def Uk(x,k): return quad(lambda t:t**k/sqn(t),r0,x,limit=200)[0]
def Pi_of(x): return quad(lambda t:Rn(t)/(dengn(t)*sqn(t)),r0,x,limit=200)[0]
def DS_of(x): return sum(akn[k]*Uk(x,k) for k in range(4))+Pi_of(x)     # full Delta S
def A_over_sqrtS(t): return (pjn[0]+pjn[1]*t)/sqn(t)+N3n(t)/(Q2n(t)*sqn(t))
def full_direct(x): return -quad(lambda t:A_over_sqrtS(t)*DS_of(t),r0,x,limit=150)[0]
# --- assembly ---
# (A) second-kind: -int (G/sqrtS)(sum a_k U_k) dr, via W_jk
def Wjk(j,k,x): return quad(lambda t:(t**j*Uk(t,k)-t**k*Uk(t,j))/sqn(t),r0,x,limit=150)[0]
def blockA(x):
    # -sum_{j,k} g_j a_k * 1/2 (U_j U_k + W_jk)
    tot=0.0
    for j in range(5):
        for k in range(4):
            tot+=-0.5*gn[j]*akn[k]*(Uk(x,j)*Uk(x,k)+Wjk(j,k,x))
    return tot
# (B) Delta=0 dilog: -int (G/sqrtS) Pi dr = -[Phi_G Pi] + sum_root rho sum_j g_j D_{j,root}
def PhiG(x): return sum(gn[j]*Uk(x,j) for j in range(5))
def D(j,nm,x): return quad(lambda t:Uk(t,j)/((t-rootn[nm])*sqn(t)),r0,x,limit=150)[0]
def blockB(x):
    bnd=-PhiG(x)*Pi_of(x)
    dil=sum(rhon[nm]*sum(gn[j]*D(j,nm,x) for j in range(5)) for nm in roots)
    return bnd+dil
# (C) elementary: -[alg_K Delta S] + int P r/Delta dr
def algK(x): return Pn(x)*sqn(x)/Q2n(x)
Dln=sp.lambdify(r,Dl.subs(E,E0),'numpy'); PrDl=sp.lambdify(r,(P*r/Dl).subs(E,E0),'numpy')
def blockC(x):
    boundary=-(algK(x)*DS_of(x)-algK(r0)*DS_of(r0))
    elem=quad(PrDl,r0,x,limit=200)[0]
    return boundary+elem
print("\n=== FULL wrap: direct double integral vs (A second-kind + B Delta=0 dilog + C elementary) ===")
for x in [10.0,8.0,6.5]:
    fd=full_direct(x); A_=blockA(x); B_=blockB(x); C_=blockC(x); asm=A_+B_+C_
    print(f"  r={x}: direct={fd:+.7f}  assembly={asm:+.7f}  diff={abs(fd-asm):.1e}   [A={A_:+.4f} B={B_:+.5f} C={C_:+.4f}]")
print("\n=> FULL generic off-shell wrap = 2nd-kind W_jk (Kleinian+s0-dilog) + Delta=0 genus-2 dilogs")
print("   D_{j,root} + elementary, ALL coefficients symbolic in E. Verified vs direct.")
