# -*- coding: utf-8 -*-
# BLOCK 3 of the generic off-shell wrap (t-branch): the THIRD-KIND-at-Delta=0 piece.
# Delta S has a third-kind part  Pi = int R/(deng*sqrt S) dr,  deng = Delta*Emu (poles at
# Delta=0 -> r_pm, and Emu=0 -> r_E). Partial fractions R/(Delta*Emu)=sum_root rho_root/(r-root)
# with SYMBOLIC residues => Pi = sum_root rho_root * Pi_root,  Pi_root=int dr/((r-root) sqrt S)
# the canonical genus-2 THIRD-KIND letters. Block 3 = -int (A_poly/sqrt S) * Pi dr couples the
# 2nd-kind kernel Phi_p to these third-kind letters -> the genus-2 DILOGARITHM at Delta=0.
# STEP 1 here: verify Pi = sum rho_root Pi_root (symbolic residues) numerically.
import sympy as sp, numpy as np
from scipy.integrate import quad

r,E = sp.symbols('r E', positive=True)
M,a,J = sp.Rational(1), sp.Rational(9,10), sp.Integer(6); E0=sp.Rational(7,5)
Emu=(E**2-1)*r+2*M; Dl=r**2-2*M*r+a**2
Q2=(2*E**2*J**2*M*r-E**2*J**2*r**2-4*E**2*J*M*a*r+2*E**2*M*a**2*r+E**2*a**2*r**2
    +E**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2-8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
S=sp.expand(r*Emu*Q2)
deng=sp.expand(Dl*Emu)
Qd,R=sp.div(sp.Poly(sp.expand(S),r),sp.Poly(deng,r)); R=R.as_expr()
print("R =",sp.factor(R)); print("deng = Delta*Emu, roots: r_+,r_- (Delta), r_E (Emu)")

# partial fractions R/(Delta*Emu) over its 3 roots (symbolic in E via apart)
pf=sp.apart(R/deng, r)
print("\npartial fractions R/(Delta*Emu):"); sp.pprint(pf)
# residues at each simple pole: rho_root = R/(d/dr[deng]) at root
dengp=sp.diff(deng,r)
roots={'r+':(M+sp.sqrt(M**2-a**2)), 'r-':(M-sp.sqrt(M**2-a**2)), 'rE':sp.simplify(-2*M/(E**2-1))}
rho={nm:sp.simplify((R/dengp).subs(r,rt)) for nm,rt in roots.items()}
print("\nsymbolic residues rho_root = R/deng'(root):")
for nm in roots: print(f"  rho[{nm}] =",sp.simplify(rho[nm]))

# ---- numeric verification: Pi = sum_root rho_root Pi_root ----
Sn=sp.lambdify(r,S.subs(E,E0),'numpy'); sqn=lambda x:np.sqrt(max(Sn(x),0.0))
dengn=sp.lambdify(r,deng.subs(E,E0),'numpy'); Rn=sp.lambdify(r,R.subs(E,E0),'numpy')
rootn={nm:complex(rt.subs(E,E0)) for nm,rt in roots.items()}
rhon={nm:complex(rho[nm].subs(E,E0)) for nm in roots}
print("\nnumeric roots:",{k:(f"{v.real:.4f}+{v.imag:.4f}i" if abs(v.imag)>1e-9 else f"{v.real:.4f}") for k,v in rootn.items()})
r0=11.0
def Pi_direct(x): return quad(lambda t: Rn(t)/(dengn(t)*sqn(t)), r0, x, limit=200)[0]
def Pi_root(x,rt):  # int dr/((r-root) sqrt S); root may be complex (r+,r- complex here since a<M? a=0.9<1=M so real)
    return quad(lambda t: 1.0/((t-rt)*sqn(t)), r0, x, limit=200)[0]
print("\nverify Pi_direct == sum_root rho_root Pi_root:")
for x in [10.0,8.0,6.5]:
    direct=Pi_direct(x)
    asm=sum(rhon[nm]*Pi_root(x,rootn[nm].real if abs(rootn[nm].imag)<1e-9 else rootn[nm]) for nm in roots)
    print(f"  r={x}: Pi_direct={direct:+.8f}  sum rho*Pi_root={asm.real:+.8f}  diff={abs(direct-asm):.1e}")
print("\n=> Delta S third-kind part Pi = sum_root rho_root(E) Pi_root, symbolic residues. Next: couple to Phi_p.")
