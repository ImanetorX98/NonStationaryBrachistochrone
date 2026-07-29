import sympy as sp, numpy as np
from scipy.optimize import brentq
r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
f=1-2*m/r; w=E**2-f; DE=(E**2-1)*r+2*m
Q2v=(2*E**2*J**2*m*r-E**2*J**2*r**2+E**2*r**4+4*J**2*m**2-4*J**2*m*r+J**2*r**2)
S_v=sp.expand(r*DE*Q2v)
Hv=pr*(f-E**2)-1+sp.sqrt(w)*sp.sqrt(E**2*pr**2+J**2/r**2)
la=lambda ex: sp.lambdify((r,pr,m,E,J),ex,'numpy')
Hn=la(Hv); Hp=la(sp.diff(Hv,pr)); Gpr=la(sp.diff(sp.diff(Hv,J)/sp.diff(Hv,pr),pr))
Svn=sp.lambdify((r,m,E,J),S_v,'numpy'); Q2vn=sp.lambdify((r,m,E,J),Q2v,'numpy'); DEn=sp.lambdify((r,m,E,J),DE,'numpy')
mv,Ev,Jv=1.0,1.4,2.5
def pr0(R):
    g=lambda p:Hn(R,p,mv,Ev,Jv); ps=np.linspace(-40,-1e-4,4000); vs=[g(p) for p in ps]
    s=np.where(np.diff(np.sign(vs)))[0]; return brentq(g,ps[s[0]],ps[s[0]+1])
print("A(r)=kernel*sqrt(S_v), and candidate ratios (pure 2nd-kind => A(r) rational):")
for R in [10.0,8.0,6.0,5.0]:
    p=pr0(R); k=Gpr(R,p,mv,Ev,Jv)/Hp(R,p,mv,Ev,Jv); A=k*np.sqrt(Svn(R,mv,Ev,Jv))
    q=Q2vn(R,mv,Ev,Jv); de=DEn(R,mv,Ev,Jv)
    print(f"  r={R}: A={A:+.5f}  A*Q2v/(DE)={A*q/de:+.4f}  A*Q2v/(E^2 DE)={A*q/(Ev**2*de):+.4f}  /r^2={A*q/(Ev**2*de)/R**2:+.5f} /r^3={A*q/(Ev**2*de)/R**3:+.6f}")
