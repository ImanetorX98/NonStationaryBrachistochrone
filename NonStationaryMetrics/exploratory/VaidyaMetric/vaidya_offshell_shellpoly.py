# -*- coding: utf-8 -*-
# Vaidya tau-branch off-shell, SHELL-POLYNOMIAL machine (no sympy sqrt hell).
# Analytic simplifications on the clean H_tau shell (verified below):
#   G=J/(Delta p_r),  H_pr=(r-2m)DE p_r/(E r^2),  H_m=C0(r)+C2(r) p_r^2,
#   kernel  d_pr G/H_pr = A_V/sqrt(S),  A_V=-J E r^3 DE/Q3,
#   inner   m H_m/H_pr  = P_inner(r)/sqrt(S).
# S=r(r-2m)DE Q3, D0=r(r-2m)DE, Q3=r^2(r-2m)-J^2 DE, DE=(E^2-1)r+2m ; p_r0^2=S/D0^2.
# STEP: build A_V, P_inner symbolically (rational) and VERIFY vs the numerical on-shell values.
import sympy as sp, numpy as np
from scipy.optimize import brentq

r,pr,m,E,J=sp.symbols('r p_r m E J',positive=True)
Dl=r*(r-2*m); f=(r-2*m)/r; DE=(E**2-1)*r+2*m; v=DE/(E**2*r)
D0=sp.expand(r*(r-2*m)*DE); Q3=sp.expand(r**2*(r-2*m)-J**2*DE); S=sp.expand(D0*Q3)
# analytic kernel and H-derivatives
# ingoing shell: p_r0 = -sqrt(S)/D0  (p_r<0). kernel~1/p_r^3, inner~1/p_r both ODD -> sign from ingoing.
A_V=sp.simplify(J*E*r**3*DE/Q3)                      # kernel numerator (ingoing): kernel=A_V/sqrt(S)
C0=sp.diff(Dl*v,m)*f/(2*E*Dl*v)+2/(r*E)              # p_r-independent part of H_m
C2=-Dl*v*E/(r**3*f)                                  # p_r^2 coefficient of H_m
# inner: m H_m/H_pr = mEr^2/((r-2m)DE) (C0/p_r + C2 p_r); ingoing 1/p_r=-D0/sqrtS, p_r=-sqrtS/D0
# => P_inner = -mEr^2 (C0 D0^2 + C2 S)/((r-2m)DE D0)
P_inner=sp.simplify(-m*E*r**2*(C0*D0**2+C2*S)/((r-2*m)*DE*D0))
print("A_V (kernel numerator) =",sp.factor(A_V))
print("P_inner (inner-letter numerator, over sqrt S) =",sp.factor(P_inner))

# ---------- numerical verification against the FULL clean Hamiltonian ----------
Dln=lambda R:R*(R-2*m);
Htau=sp.sqrt(Dl*v/r**2)*sp.sqrt((Dl/r**2)*pr**2+J**2/r**2)-f/E
Hn=sp.lambdify((r,pr,m,E,J),Htau,'numpy')
Hp=sp.lambdify((r,pr,m,E,J),sp.diff(Htau,pr),'numpy')
Hm=sp.lambdify((r,pr,m,E,J),sp.diff(Htau,m),'numpy')
Gexpr=sp.diff(Htau,J)/sp.diff(Htau,pr); Gpr=sp.lambdify((r,pr,m,E,J),sp.diff(Gexpr,pr),'numpy')
Sn=sp.lambdify((r,m,E,J),S,'numpy'); A_Vn=sp.lambdify((r,m,E,J),A_V,'numpy'); P_in=sp.lambdify((r,m,E,J),P_inner,'numpy')
def pr0(R,mv,Ev,Jv):
    g=lambda p:Hn(R,p,mv,Ev,Jv); ps=np.linspace(-60,-1e-4,4000); vals=[g(p) for p in ps]
    s=np.where(np.diff(np.sign(vals)))[0]; return brentq(g,ps[s[0]],ps[s[0]+1])
print("\nverify kernel = A_V/sqrt(S) and inner = P_inner/sqrt(S) (mv=1,E=1.4,J=2.5):")
mv,Ev,Jv=1.0,1.4,2.5
for R in [10.0,8.0,6.0]:
    p=pr0(R,mv,Ev,Jv); sq=np.sqrt(Sn(R,mv,Ev,Jv))
    kern_num=Gpr(R,p,mv,Ev,Jv)/Hp(R,p,mv,Ev,Jv); kern_an=A_Vn(R,mv,Ev,Jv)/sq
    inn_num=mv*Hm(R,p,mv,Ev,Jv)/Hp(R,p,mv,Ev,Jv); inn_an=P_in(R,mv,Ev,Jv)/sq
    print(f"  r={R}: kernel num={kern_num:+.6f} an={kern_an:+.6f} d={abs(kern_num-kern_an):.1e} ;"
          f" inner num={inn_num:+.6f} an={inn_an:+.6f} d={abs(inn_num-inn_an):.1e}")
print("\n=> if diffs ~0: clean rational A_V and P_inner confirmed; proceed to A+B+C reduction/assembly.")
