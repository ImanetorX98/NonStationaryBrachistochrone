import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
f=1-2*m/r; w=E**2-f; DE=(E**2-1)*r+2*m; K1=sp.simplify(f-E**2); K2=sp.simplify(w*E**2-K1**2)
Q2v=(2*E**2*J**2*m*r-E**2*J**2*r**2+E**2*r**4+4*J**2*m**2-4*J**2*m*r+J**2*r**2)
S_v=sp.expand(r*DE*Q2v); D_v=sp.simplify(K1**2-K2*(w*J**2/r**2-1)); P=sp.simplify(K2+K1**2)
Hv=pr*(f-E**2)-1+sp.sqrt(w)*sp.sqrt(E**2*pr**2+J**2/r**2)
la=lambda ex: sp.lambdify((r,pr,m,E,J),ex,'numpy')
Hn=la(Hv); Hp=la(sp.diff(Hv,pr)); Hm=la(sp.diff(Hv,m)); Gpr=la(sp.diff(sp.diff(Hv,J)/sp.diff(Hv,pr),pr))
# verified pieces
Akern=sp.lambdify((r,m,E,J),E**2*J*DE*r**4/Q2v,'numpy')     # kernel=Akern/sqrt(S_v) (pure 2nd-kind)
ELEMp=sp.simplify(-2*w*(P+K1**2)-2*P*K1); SECp=sp.simplify(2*w*K1*D_v+2*w*K1*P+P**2+K1**2*D_v)
elem_inn=sp.lambdify((r,m,E,J),m*ELEMp/(r*K2**2*w),'numpy')
A_inn=sp.lambdify((r,m,E,J),m*SECp/(r*K2**2*w),'numpy')
Svn=sp.lambdify((r,m,E,J),S_v,'numpy'); Dvn=sp.lambdify((r,m,E,J),D_v,'numpy'); K1n=sp.lambdify((r,m,E,J),K1,'numpy'); K2n=sp.lambdify((r,m,E,J),K2,'numpy')
mv,Ev,Jv,r0=1.0,1.4,2.5,12.0
def pr0(R):
    g=lambda p:Hn(R,p,mv,Ev,Jv); ps=np.linspace(-40,-1e-4,4000); vs=[g(p) for p in ps]
    s=np.where(np.diff(np.sign(vs)))[0]; return brentq(g,ps[s[0]],ps[s[0]+1])
sq=lambda x:np.sqrt(max(Svn(x,mv,Ev,Jv),0.0))
def kern(x): p=pr0(x); return Gpr(x,p,mv,Ev,Jv)/Hp(x,p,mv,Ev,Jv)
def inner(x): p=pr0(x); return mv*Hm(x,p,mv,Ev,Jv)/Hp(x,p,mv,Ev,Jv)
# verify kernel pure 2nd-kind
print("kernel = Akern/sqrt(S_v)? ", [f"{abs(kern(x)-Akern(x,mv,Ev,Jv)/sq(x)):.1e}" for x in [10,8,6]])
def sgn(x): p=pr0(x); return np.sign(K1n(x,mv,Ev,Jv)+p*K2n(x,mv,Ev,Jv))
def inner_elem(x): return elem_inn(x,mv,Ev,Jv)
def inner_2nd(x): return A_inn(x,mv,Ev,Jv)*sgn(x)/np.sqrt(Dvn(x,mv,Ev,Jv))
print("inner = elem+2nd? ", [f"{abs(inner(x)-(inner_elem(x)+inner_2nd(x))):.1e}" for x in [10,8,6]])
def Sig_elem(x): return quad(inner_elem,r0,x,limit=200)[0]
def Sig_2nd(x): return quad(inner_2nd,r0,x,limit=200)[0]
def Sig(x): return quad(inner,r0,x,limit=200)[0]
def wrap_direct(x): return -quad(lambda t:kern(t)*Sig(t),r0,x,limit=120)[0]
def block1(x): return -quad(lambda t:kern(t)*Sig_2nd(t),r0,x,limit=120)[0]   # 2nd x 2nd (genus-2)
def block2(x): return -quad(lambda t:kern(t)*Sig_elem(t),r0,x,limit=120)[0]  # 2nd x elem (log x Abelian)
print("\nwrap decomposition (direct = block1[genus-2] + block2[log x Abelian]):")
for x in [11.0,10.0,9.0]:
    d=wrap_direct(x); b1=block1(x); b2=block2(x)
    print(f"  r={x}: direct={d:+.6f}  block1={b1:+.6f}  block2={b2:+.6f}  sum-diff={abs(d-b1-b2):.1e}  block2/direct={b2/d*100:+.1f}%")
