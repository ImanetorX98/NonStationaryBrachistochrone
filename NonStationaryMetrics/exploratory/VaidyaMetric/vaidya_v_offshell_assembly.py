# -*- coding: utf-8 -*-
# Vaidya v-branch off-shell: FULL wrap = second-kind block (A_kernel_v/sqrt S_v) + elementary block.
# S_v=r DE Q2v (genus-2), A_kernel_v=-E^2 J DE r^4/Q2v, elementary_v=-J DE r/Q2v.
# inner letter also splits: m H_m/H_pr = P2/sqrt(S_v) + P_elem (rational). Verify direct wrap == assembly.
import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
f=1-2*m/r; w=E**2-f; DE=(E**2-1)*r+2*m
Q2v=(2*E**2*J**2*m*r-E**2*J**2*r**2+E**2*r**4+4*J**2*m**2-4*J**2*m*r+J**2*r**2)
S_v=sp.expand(r*DE*Q2v)
Hv=pr*(f-E**2)-1+sp.sqrt(w)*sp.sqrt(E**2*pr**2+J**2/r**2)
la=lambda ex: sp.lambdify((r,pr,m,E,J),ex,'numpy')
Hn=la(Hv); Hp=la(sp.diff(Hv,pr)); Hm=la(sp.diff(Hv,m))
Gv=sp.diff(Hv,J)/sp.diff(Hv,pr); Gpr=la(sp.diff(Gv,pr))
mv,Ev,Jv,r0=1.0,1.4,2.5,12.0
def pr0(R):
    g=lambda p:Hn(R,p,mv,Ev,Jv); ps=np.linspace(-40,-1e-4,4000); vs=[g(p) for p in ps]
    s=np.where(np.diff(np.sign(vs)))[0]; return brentq(g,ps[s[0]],ps[s[0]+1])
Svn=la(S_v)(0,0,0,0) if False else sp.lambdify((r,m,E,J),S_v,'numpy')
sq=lambda x: np.sqrt(max(Svn(x,mv,Ev,Jv),0.0))
def kern(x): return Gpr(x,pr0(x),mv,Ev,Jv)/Hp(x,pr0(x),mv,Ev,Jv)
def inner(x): return mv*Hm(x,pr0(x),mv,Ev,Jv)/Hp(x,pr0(x),mv,Ev,Jv)
# analytic kernel split (physical ingoing branch): kernel = A_kernel_v/sqrt(S_v) + elementary_v
A_kern=sp.lambdify((r,m,E,J),(-E**2*J*DE*r**4/Q2v),'numpy')
elem=sp.lambdify((r,m,E,J),(-J*DE*r/Q2v),'numpy')
print("verify kernel = A_kernel_v/sqrt(S_v) + elementary_v (ingoing):")
for x in [10.0,8.0,6.0]:
    kn=kern(x); ka=A_kern(x,mv,Ev,Jv)/sq(x)+elem(x,mv,Ev,Jv)
    print(f"  r={x}: kernel num={kn:+.6f} split={ka:+.6f} d={abs(kn-ka):.1e}")
# turning point (Q2v=0 physical root)
rr=np.linspace(2.01,r0,6000); q=[Q2v.subs({m:mv,E:Ev,J:Jv,r:x}) for x in rr]
qn=np.array([float(x) for x in q]); idx=np.where(np.diff(np.sign(qn)))[0]
rt=max(rr[i] for i in idx) if len(idx) else 2.1; xf=rt+0.6
def Sigma(x): return quad(inner,r0,x,limit=200)[0]
def wrap_direct(x): return -quad(lambda t:kern(t)*Sigma(t),r0,x,limit=120)[0]
# second-kind-only assembly proxy: -int (A_kernel_v/sqrtS) Sigma  (the genus-2 part)
def wrap_2nd(x): return -quad(lambda t:(A_kern(t,mv,Ev,Jv)/sq(t))*Sigma(t),r0,x,limit=120)[0]
def wrap_elem(x): return -quad(lambda t:elem(t,mv,Ev,Jv)*Sigma(t),r0,x,limit=120)[0]
print(f"\nturning={rt:.4f}. wrap decomposition (direct = 2nd-kind + elementary blocks):")
for x in [11.0,10.0,9.0]:
    d=wrap_direct(x); w2=wrap_2nd(x); we=wrap_elem(x)
    print(f"  r={x}: direct={d:+.6f}  2nd-kind={w2:+.6f}  elem-block={we:+.6f}  sum={w2+we:+.6f}  d={abs(d-(w2+we)):.1e}")
print("\n=> elem-block = -int elementary_v * Sigma dr (rational x transcendental). Next: reduce elem-block")
print("   by IBP -> boundary + int (elementary-log) x (inner) = log x Abelian weight-2 letters.")
