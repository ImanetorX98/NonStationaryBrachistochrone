# -*- coding: utf-8 -*-
# Vaidya tau-branch off-shell FULL closed form (shell-polynomial machine), verified vs direct.
# kernel A_V/sqrt(S), A_V=J E r^3 DE/Q3 = A_poly + N3/Q3 (Q3 branch poles -> Hermite);
# inner Sigma=int P_inner/sqrt(S), P_inner=-m N4/deng_V, deng_V=(r-2m)DE (third-kind @ r=2m, DE=0).
#   Wrap Phi=-int(A_V/sqrtS)*Sigma dr = A(2nd-kind W_jk)+B(horizon/DE dilogs)+C(Hermite elementary).
# S=r(r-2m)DE*Q3, T_V=r(r-2m)DE, Q3=r^2(r-2m)-J^2 DE, DE=(E^2-1)r+2m.
import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

r,m,E,J=sp.symbols('r m E J',positive=True)
mv,Ev,Jv,r0=1.0,sp.Rational(7,5),sp.Rational(5,2),12.0
DE=(E**2-1)*r+2*m; T_V=sp.expand(r*(r-2*m)*DE); Q3=sp.expand(r**2*(r-2*m)-J**2*DE)
S=sp.expand(T_V*Q3); deng_V=sp.expand((r-2*m)*DE)
# kernel A_V = J E r^3 DE / Q3 = A_poly + N3/Q3
Ap,N3=sp.div(sp.Poly(sp.expand(J*E*r**3*DE),r),sp.Poly(Q3,r)); N3=N3.as_expr()
p_j=[c for c in Ap.all_coeffs()[::-1]]
# inner numerator: P_inner = -m*N4/deng_V ; N4 the quartic
N4=(E**4*J**2*r**2+4*E**2*J**2*m*r-2*E**2*J**2*r**2+4*J**2*m**2-4*J**2*m*r+J**2*r**2+4*m**2*r**2-4*m*r**3+r**4)
Qd,Rrem=sp.div(sp.Poly(sp.expand(-m*N4),r),sp.Poly(deng_V,r)); Rrem=Rrem.as_expr()
a_k=[sp.simplify(c) for c in Qd.all_coeffs()[::-1]]          # inner poly part: Sigma_poly=sum a_k U_k
# third-kind residues of (-m N4)/deng_V at roots of deng_V (r=2m, DE=0)
roots={'2m':2*m,'DE':sp.simplify(-2*m/(E**2-1))}
dengp=sp.diff(deng_V,r); rho={nm:sp.simplify((Rrem/dengp).subs(r,rt)) for nm,rt in roots.items()}
# Hermite for N3/(Q3 sqrtS): S=Q3 T_V ; P=-2 N3 (Q3' T_V)^{-1} mod Q3
inv=sp.invert(sp.Poly(sp.expand(sp.diff(Q3,r)*T_V),r),sp.Poly(Q3,r))
Pher=sp.expand(sp.rem(sp.expand(-2*N3*inv.as_expr()),Q3,r))
RemH=sp.div(sp.Poly(sp.expand(2*N3-2*sp.diff(Pher,r)*Q3*T_V-Pher*Q3*sp.diff(T_V,r)+Pher*sp.diff(Q3,r)*T_V),r),sp.Poly(sp.expand(2*Q3),r))
print("Hermite exact-div remainder:",sp.simplify(RemH[1].as_expr()))
remc=RemH[0].as_expr(); remk=[sp.simplify(c) for c in sp.Poly(remc,r).all_coeffs()[::-1]]
g=[sp.simplify((p_j[k] if k<len(p_j) else 0)+(remk[k] if k<len(remk) else 0)) for k in range(max(len(p_j),len(remk)))]
print("kernel g_k (poly+Hermite) =",g)
print("inner a_k =",a_k,"  rho:",{k:sp.simplify(v) for k,v in rho.items()})

# ================= numeric verification =================
sub={m:mv,E:Ev,J:Jv}
Sn=sp.lambdify(r,S.subs(sub),'numpy'); sq=lambda x:np.sqrt(max(Sn(x),0.0))
Q3n=sp.lambdify(r,Q3.subs(sub),'numpy'); dengn=sp.lambdify(r,deng_V.subs(sub),'numpy')
N3n=sp.lambdify(r,N3.subs(sub),'numpy'); Rremn=sp.lambdify(r,Rrem.subs(sub),'numpy'); Phern=sp.lambdify(r,Pher.subs(sub),'numpy')
A_Vn=sp.lambdify(r,(J*E*r**3*DE/Q3).subs(sub),'numpy'); P_inn=sp.lambdify(r,(-m*N4/deng_V).subs(sub),'numpy')
gn=[float(x.subs(sub)) for x in g]; akn=[float(x.subs(sub)) for x in a_k]; pjn=[float(x.subs(sub)) for x in p_j]
rootn={nm:float(rt.subs(sub)) for nm,rt in roots.items()}; rhon={nm:float(rho[nm].subs(sub)) for nm in roots}
rt_turn=max([x for x in np.roots([float(c) for c in sp.Poly(Q3.subs(sub),r).all_coeffs()]) if abs(np.imag(x))<1e-9 and 2.0<np.real(x)<r0],key=lambda z:np.real(z)).real
xf=rt_turn+0.6
def Uk(x,k): return quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def Sigma(x): return quad(lambda t:P_inn(t)/sq(t),r0,x,limit=200)[0]
def wrap_direct(x): return -quad(lambda t:(A_Vn(t)/sq(t))*Sigma(t),r0,x,limit=150)[0]
def Wjk(j,k,x): return quad(lambda t:(t**j*Uk(t,k)-t**k*Uk(t,j))/sq(t),r0,x,limit=150)[0]
def blockA(x): return sum(-0.5*gn[j]*akn[k]*(Uk(x,j)*Uk(x,k)+Wjk(j,k,x)) for j in range(len(gn)) for k in range(len(akn)))
def PhiG(x): return sum(gn[j]*Uk(x,j) for j in range(len(gn)))
def Pi_of(x): return quad(lambda t:Rremn(t)/(dengn(t)*sq(t)),r0,x,limit=200)[0]
def D(j,nm,x): return quad(lambda t:Uk(t,j)/((t-rootn[nm])*sq(t)),r0,x,limit=150)[0]
def blockB(x): return -PhiG(x)*Pi_of(x)+sum(rhon[nm]*sum(gn[j]*D(j,nm,x) for j in range(len(gn))) for nm in roots)
# elementary: Phi_C = -[algK*Sigma] + int algK*Sigma' dr, algK=Pher*sqrt(S)/Q3, Sigma'=P_inner/sqrt(S)
#   => int (Pher*sqrt S/Q3)(P_inner/sqrt S) dr = int Pher*P_inner/Q3 dr  (rational -> elementary)
def algK(x): return Phern(x)*sq(x)/Q3n(x)
def blockC(x):
    boundary=-(algK(x)*Sigma(x)-algK(r0)*Sigma(r0))
    elem=quad(lambda t: Phern(t)*P_inn(t)/Q3n(t), r0, x, limit=200)[0]
    return boundary+elem
print(f"\nturning={rt_turn:.4f}, verify wrap at x<xf={xf:.3f}:")
for x in [11.0,10.0,9.0]:
    d=wrap_direct(x); asm=blockA(x)+blockB(x)+blockC(x)
    print(f"  r={x}: direct={d:+.7f}  A+B+C={asm:+.7f}  diff={abs(d-asm):.1e}")
print("\n(if diff not ~0, block C elementary cross-term needs the Hermite-boundary integral; see TK template)")
