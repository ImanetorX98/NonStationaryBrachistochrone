# -*- coding: utf-8 -*-
# Vaidya v-branch STEP 1 (clean re-derivation): kernel AND inner-letter split, verified vs numeric.
# On shell y=K1+pr K2=+/-sqrt(D_v), rho=1-pr K1, pr=(y-K1)/K2. K1=-DE/r, K2=(r-2m)DE/r^2, w=DE/r.
#   H_pr=y/rho ;  G=wJ/(r^2 y) ;  d_pr G=-wJK2/(r^2 y^2)
#   kernel=d_pr G/H_pr=-wJK2 rho/(r^2 y^3) = -wJ(K2+K1^2)/(r^2 D_v)(1/y) + wJK1/(r^2 D_v)
#   inner =m H_m/H_pr = m/(r K2^2 w)[ ELEMpart + sign(y) SECONDpart/sqrt(D_v) ]
# with P=K2+K1^2, ELEMpart=-2w(P+K1^2)-2PK1, SECONDpart=2wK1 D_v+2wK1 P+P^2+K1^2 D_v.
# sign(y)=-1 for the ingoing branch (pr<0). Verify both splits to ~1e-14.
import sympy as sp, numpy as np
from scipy.optimize import brentq
r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
f=1-2*m/r; w=E**2-f; DE=(E**2-1)*r+2*m
K1=sp.simplify(f-E**2); K2=sp.simplify(w*E**2-K1**2)
Q2v=(2*E**2*J**2*m*r-E**2*J**2*r**2+E**2*r**4+4*J**2*m**2-4*J**2*m*r+J**2*r**2)
S_v=sp.expand(r*DE*Q2v); D_v=sp.simplify(K1**2-K2*(w*J**2/r**2-1))
P=sp.simplify(K2+K1**2)
# CORRECTED kernel (off-shell d_pr G on shell): kernel = -w^2 J E^2/(r^4 y^3) = -sign(y) w^2 J E^2/(r^4 D_v^{3/2})
# => PURE second-kind (no elementary part). A_kern_coeff = -w^2 J E^2/r^4 (over D_v^{3/2}, * sign(y)).
A_kern=sp.simplify(-w**2*J*E**2/r**4)         # kernel 2nd-kind coeff (over D_v^{3/2}, * sign(y))
elem_kern=sp.Integer(0)                        # NO elementary part in the kernel
ELEMp=sp.simplify(-2*w*(P+K1**2)-2*P*K1)
SECp=sp.simplify(2*w*K1*D_v+2*w*K1*P+P**2+K1**2*D_v)
elem_inn=sp.simplify(m*ELEMp/(r*K2**2*w))    # inner elementary
A_inn=sp.simplify(m*SECp/(r*K2**2*w))        # inner 2nd-kind: coeff of sign(y)/sqrt(D_v)
# numeric verification
Hv=pr*(f-E**2)-1+sp.sqrt(w)*sp.sqrt(E**2*pr**2+J**2/r**2)
la=lambda ex: sp.lambdify((r,pr,m,E,J),ex,'numpy')
Hn=la(Hv); Hp=la(sp.diff(Hv,pr)); Hm=la(sp.diff(Hv,m)); Gpr=la(sp.diff(sp.diff(Hv,J)/sp.diff(Hv,pr),pr))
Dvn=sp.lambdify((r,m,E,J),D_v,'numpy'); K1n=sp.lambdify((r,m,E,J),K1,'numpy'); K2n=sp.lambdify((r,m,E,J),K2,'numpy')
Akn=sp.lambdify((r,m,E,J),A_kern,'numpy'); ekn=sp.lambdify((r,m,E,J),elem_kern,'numpy')
Ain=sp.lambdify((r,m,E,J),A_inn,'numpy'); ein=sp.lambdify((r,m,E,J),elem_inn,'numpy')
mv,Ev,Jv=1.0,1.4,2.5
def pr0(R):
    g=lambda p:Hn(R,p,mv,Ev,Jv); ps=np.linspace(-40,-1e-4,4000); vs=[g(p) for p in ps]
    s=np.where(np.diff(np.sign(vs)))[0]; return brentq(g,ps[s[0]],ps[s[0]+1])
print("verify kernel and inner splits (ingoing, sign(y)=-1):")
for R in [10.0,8.0,6.0]:
    p=pr0(R); y=K1n(R,mv,Ev,Jv)+p*K2n(R,mv,Ev,Jv); sgn=np.sign(y); Dv=Dvn(R,mv,Ev,Jv); sD=np.sqrt(Dv)
    kn=Gpr(R,p,mv,Ev,Jv)/Hp(R,p,mv,Ev,Jv)
    ka=Akn(R,mv,Ev,Jv)*sgn/sD**3+ekn(R,mv,Ev,Jv)          # kernel = A_kern*sign(y)/D_v^{3/2} (pure 2nd-kind)
    inn=mv*Hm(R,p,mv,Ev,Jv)/Hp(R,p,mv,Ev,Jv)
    ia=ein(R,mv,Ev,Jv)+Ain(R,mv,Ev,Jv)*sgn/sD             # inner = elem_inn + A_inn*sgn/sqrtDv
    print(f"  r={R}: y={y:+.4f}(sgn={sgn:+.0f})  kernel num={kn:+.6f} split={ka:+.6f} d={abs(kn-ka):.1e}"
          f" | inner num={inn:+.6f} split={ia:+.6f} d={abs(inn-ia):.1e}")
print("\n=> if diffs ~0: kernel & inner both split as ELEM + sign(y)*2nd-kind/sqrt(D_v). Next: 4-block wrap.")
