# -*- coding: utf-8 -*-
# BLOCK 3 assembly (t-branch off-shell, Delta=0 third-kind / genus-2 dilogarithm).
#   Block3 = -int (A_poly/sqrt S) * Pi dr ,  Pi = sum_root rho_root int dr/((r-root) sqrt S)
#   IBP: Block3 = -[Phi_p Pi] + sum_root rho_root sum_j p_j D_{j,root},
#     Phi_p = sum_j p_j U_j ,  D_{j,root} = int U_j/((r-root) sqrt S) dr  (genus-2 dilog at r_root).
# rho_root symbolic in E (verified in offshell_tbranch_thirdkind_Delta.py, rho_E=0).
# Verify assembly == direct block value (math level).
import sympy as sp, numpy as np
from scipy.integrate import quad

r,E=sp.symbols('r E',positive=True)
M,a,J=sp.Rational(1),sp.Rational(9,10),sp.Integer(6); E0=sp.Rational(7,5)
Emu=(E**2-1)*r+2*M; Dl=r**2-2*M*r+a**2
Q2=(2*E**2*J**2*M*r-E**2*J**2*r**2-4*E**2*J*M*a*r+2*E**2*M*a**2*r+E**2*a**2*r**2
    +E**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2-8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
S=sp.expand(r*Emu*Q2); deng=sp.expand(Dl*Emu)
Qp,N3=sp.div(sp.Poly(sp.expand(E**2*J*r**4*Emu),r),sp.Poly(sp.expand(Q2),r))
p_j=[c for c in Qp.all_coeffs()[::-1]]            # p_0,p_1
Qd,R=sp.div(sp.Poly(sp.expand(S),r),sp.Poly(deng,r)); R=R.as_expr()
dengp=sp.diff(deng,r)
roots={'r+':M+sp.sqrt(M**2-a**2),'r-':M-sp.sqrt(M**2-a**2)}
rho={nm:sp.simplify((R/dengp).subs(r,rt)) for nm,rt in roots.items()}

Sn=sp.lambdify(r,S.subs(E,E0),'numpy'); sqn=lambda x:np.sqrt(max(Sn(x),0.0))
dengn=sp.lambdify(r,deng.subs(E,E0),'numpy'); Rn=sp.lambdify(r,R.subs(E,E0),'numpy')
pjn=[float(x.subs(E,E0)) for x in p_j]
rootn={nm:float(rt.subs(E,E0)) for nm,rt in roots.items()}
rhon={nm:float(rho[nm].subs(E,E0)) for nm in roots}
r0=11.0
def Uk(x,k): return quad(lambda t:t**k/sqn(t),r0,x,limit=200)[0]
def Pi_of(x): return quad(lambda t:Rn(t)/(dengn(t)*sqn(t)),r0,x,limit=200)[0]
def Phi_p(x): return sum(pjn[j]*Uk(x,j) for j in range(2))
# direct block 3
def block3_direct(x): return -quad(lambda t:((pjn[0]+pjn[1]*t)/sqn(t))*Pi_of(t),r0,x,limit=150)[0]
# genus-2 dilog letters D_{j,root} = int U_j/((r-root) sqrt S) dr
def D(j,nm,x): return quad(lambda t: Uk(t,j)/((t-rootn[nm])*sqn(t)),r0,x,limit=150)[0]
def block3_asm(x):
    bnd=-Phi_p(x)*Pi_of(x)   # -[Phi_p Pi]_{r0}^{x}, Phi_p(r0)=Pi(r0)=0
    dil=sum(rhon[nm]*sum(pjn[j]*D(j,nm,x) for j in range(2)) for nm in roots)
    return bnd+dil
print("Block 3 (Delta=0 third-kind): direct vs [boundary + genus-2 dilog assembly]")
for x in [10.0,8.0,6.5]:
    d=block3_direct(x); asm=block3_asm(x)
    print(f"  r={x}: direct={d:+.8f}  assembly={asm:+.8f}  diff={abs(d-asm):.1e}")
print("\n=> Block3 = -[Phi_p Pi] + sum_root rho_root(E) sum_j p_j(E) D_{j,root},")
print("   D_{j,root} = genus-2 dilogarithm at the Delta=0 seed-Kerr null point r_root. Symbolic coeff.")
