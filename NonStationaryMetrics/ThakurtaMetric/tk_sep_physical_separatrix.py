# -*- coding: utf-8 -*-
# Physical TK-tau separatrix from the GENUINE branch quartic Q2 (kerr/wrap_assembly.py),
# not the simplified r*Delta-J^2*DE. Determine Jc, r_d and whether r_d is ON the
# physical integration path [turning, r0] -> decides if off-shell wrap has an on-path
# power divergence (open counterterm) or an off-path branch pole (closed like on-shell).
import sympy as sp, numpy as np
from scipy.integrate import quad
r,J=sp.symbols('r J')
M=1; a=sp.Rational(9,10); E=sp.Rational(6,5)
Q2=(2*E**2*J**2*M*r - E**2*J**2*r**2 - 4*E**2*J*M*a*r + 2*E**2*M*a**2*r
    + E**2*a**2*r**2 + E**2*r**4 + 4*J**2*M**2 - 4*J**2*M*r + J**2*r**2
    - 8*J*M**2*a + 4*J*M*a*r + 4*M**2*a**2)
DE=(E**2-1)*r+2*M
S=sp.expand(r*Q2*DE)           # the on-shell sextic (r*Emu*Q2)
# separatrix: Q2 has a double root -> disc_r(Q2)=0
res=sp.resultant(Q2,sp.diff(Q2,r),r)
cands=sorted(set(complex(s).real for s in sp.solve(res,J)
             if abs(complex(s).imag)<1e-9 and complex(s).real>0))
print("Jc candidates (Q2 double root, >0):",[f"{c:.4f}" for c in cands])
rplus=float(1+np.sqrt(1-float(a)**2)); print(f"horizon r+={rplus:.4f}")
r0=20.0
for Jc in cands:
    Sn=sp.lambdify(r,S.subs(J,sp.Float(Jc)),'numpy')
    rts=np.roots([float(c) for c in sp.Poly(sp.expand(S.subs(J,sp.Float(Jc))),r).all_coeffs()])
    realpos=sorted(float(np.real(z)) for z in rts if abs(np.imag(z))<1e-6 and np.real(z)>0)
    # double root:
    pp=[(i,j) for i in range(len(rts)) for j in range(i+1,len(rts)) if abs(rts[i]-rts[j])<1e-4]
    rd=float(np.real((rts[pp[0][0]]+rts[pp[0][1]])/2)) if pp else None
    # physical turning: largest root < r0 where S>0 for r>it
    onpath = (rd is not None and rplus < rd < r0)
    print(f" Jc={Jc:.5f} r_d={rd}  realpos={['%.3f'%x for x in realpos]}  r_d on physical path [r+,r0]? {onpath}")
