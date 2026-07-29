# -*- coding: utf-8 -*-
# Vaidya v-branch off-shell FULL closure attempt: wrap = block1[genus-2 A+B+C] + block2[2nd-kind + log-Abelian].
# kernel=A_kernel_v/sqrt(S_v), A_kernel_v=E^2 J DE r^4/Q2v (pure 2nd-kind, verified).
# inner = elem_inn(rational) + P_inner_v/sqrt(S_v).  Sigma_elem=int elem_inn = E_rat + c1 log(r-2m)+c2 log(DE).
# block1 = -int(A_kernel_v/sqrtS)Sigma_2nd  (Sigma_2nd=int P_inner_v/sqrtS) -> A+B+C.
# block2 = -int(A_kernel_v/sqrtS)Sigma_elem = -int(A_kernel_v E_rat/sqrtS) - c1 M2m - c2 MDE,
#          M^{pole}=int (A_kernel_v/sqrtS) log(r-pole) dr  (weight-2 log-Abelian letters).
import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
f=1-2*m/r; w=E**2-f; DE=(E**2-1)*r+2*m; K1=sp.simplify(f-E**2); K2=sp.simplify(w*E**2-K1**2)
Q2v=(2*E**2*J**2*m*r-E**2*J**2*r**2+E**2*r**4+4*J**2*m**2-4*J**2*m*r+J**2*r**2)
S_v=sp.expand(r*DE*Q2v); D_v=sp.simplify(K1**2-K2*(w*J**2/r**2-1)); P=sp.simplify(K2+K1**2)
ELEMp=sp.simplify(-2*w*(P+K1**2)-2*P*K1); SECp=sp.simplify(2*w*K1*D_v+2*w*K1*P+P**2+K1**2*D_v)
elem_inn=sp.simplify(m*ELEMp/(r*K2**2*w))
P_inner_v=sp.simplify(-m*SECp*r**2/(K2**2*w))       # inner 2nd-kind numerator over sqrt(S_v) (sign(y)=-1)
A_kernel_v=sp.simplify(E**2*J*DE*r**4/Q2v)
print("P_inner_v factor:",sp.factor(P_inner_v))
# Sigma_elem = int elem_inn dr  (elementary)
Sig_elem=sp.integrate(elem_inn,r)
has_log=Sig_elem.has(sp.log)
dchk=sp.simplify(sp.diff(Sig_elem,r)-elem_inn)
print("d(Sigma_elem)/dr - elem_inn =",dchk," ; Sigma_elem has log?",has_log)
print("Sigma_elem =",sp.simplify(Sig_elem))

sub={m:1.0,E:sp.Rational(7,5),J:sp.Rational(5,2)}; r0=12.0
Svn=sp.lambdify(r,S_v.subs(sub),'numpy'); sq=lambda x:np.sqrt(max(Svn(x),0.0))
Akv=sp.lambdify(r,A_kernel_v.subs(sub),'numpy'); Pinv=sp.lambdify(r,P_inner_v.subs(sub),'numpy')
eleminn_n=sp.lambdify(r,elem_inn.subs(sub),'numpy'); Sigelem_n=sp.lambdify(r,Sig_elem.subs(sub),'numpy')
# direct wrap (ground truth) via full kernel/inner
Hv=pr*(f-E**2)-1+sp.sqrt(w)*sp.sqrt(E**2*pr**2+J**2/r**2); la=lambda ex: sp.lambdify((r,pr,m,E,J),ex,'numpy')
Hn=la(Hv); Hp=la(sp.diff(Hv,pr)); Hm=la(sp.diff(Hv,m)); Gpr=la(sp.diff(sp.diff(Hv,J)/sp.diff(Hv,pr),pr))
mv,Ev,Jv=1.0,1.4,2.5
def pr0(R):
    g=lambda p:Hn(R,p,mv,Ev,Jv); ps=np.linspace(-40,-1e-4,4000); vs=[g(p) for p in ps]
    s=np.where(np.diff(np.sign(vs)))[0]; return brentq(g,ps[s[0]],ps[s[0]+1])
def kern(x): p=pr0(x); return Gpr(x,p,mv,Ev,Jv)/Hp(x,p,mv,Ev,Jv)
def inner(x): p=pr0(x); return mv*Hm(x,p,mv,Ev,Jv)/Hp(x,p,mv,Ev,Jv)
def Sig(x): return quad(inner,r0,x,limit=200)[0]
def wrap_direct(x): return -quad(lambda t:kern(t)*Sig(t),r0,x,limit=120)[0]
# block1 = -int(A_kernel_v/sqrtS) Sigma_2nd, Sigma_2nd=int P_inner_v/sqrtS
def Sig2nd(x): return quad(lambda t:Pinv(t)/sq(t),r0,x,limit=200)[0]
def block1(x): return -quad(lambda t:(Akv(t)/sq(t))*Sig2nd(t),r0,x,limit=120)[0]
# block2 = -int(A_kernel_v/sqrtS) Sigma_elem, Sigma_elem anchored at r0
S0=Sigelem_n(r0)
def block2(x): return -quad(lambda t:(Akv(t)/sq(t))*(Sigelem_n(t)-S0),r0,x,limit=120)[0]
print("\nFULL closure check: direct == block1 + block2 (block1 genus-2, block2 log-Abelian):")
for x in [11.0,10.0,9.0]:
    d=wrap_direct(x); b1=block1(x); b2=block2(x)
    print(f"  r={x}: direct={d:+.7f}  block1={b1:+.7f}  block2={b2:+.7f}  diff={abs(d-b1-b2):.1e}")
print("\n=> block1 reduces via A+B+C (tau machine); block2 via 2nd-kind(E_rat)+c1 M^2m+c2 M^DE (log-Abelian).")
