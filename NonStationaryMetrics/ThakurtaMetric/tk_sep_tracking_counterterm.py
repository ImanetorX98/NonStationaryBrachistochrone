# -*- coding: utf-8 -*-
# TK t-branch separatrix TRACKING counterterm test: does the moving double root dr_d/dlambda cancel
# the off-shell POWER pole 1/(r-r_d)^2, or is it a boundary-layer non-uniformity?
# Physical separatrix (genuine Q2): E=1.2, a=0.9 -> Jc, r_d on-path in ergosphere.
#   double root: Q2(r_d;E,Jc)=0, Q2'(r_d)=0. dilation ray dE/dl=-eps E, dJ/dl=-eps J (E_eff=E/A,J_eff=J/A).
#   dr_d/dl = -eps (E dr_d/dE + J dr_d/dJ). Compare the moving-root contribution to the pole residue.
import sympy as sp, numpy as np
r,E,J=sp.symbols('r E J',positive=True); M=1; a=sp.Rational(9,10)
Q2=(2*E**2*J**2*M*r-E**2*J**2*r**2-4*E**2*J*M*a*r+2*E**2*M*a**2*r+E**2*a**2*r**2
    +E**2*r**4+4*J**2*M**2-4*J**2*M*r+J**2*r**2-8*J*M**2*a+4*J*M*a*r+4*M**2*a**2)
E0=sp.Rational(6,5)
# separatrix Jc, r_d (double root of Q2)
res=sp.resultant(Q2,sp.diff(Q2,r),r)
Jc=sorted(c for c in [complex(s).real for s in sp.solve(res.subs(E,E0),J) if abs(complex(s).imag)<1e-9] if 2<c<5)[0]
Q2c=Q2.subs({E:E0,J:sp.Float(Jc)})
rd=sorted(set(round(float(np.real(z)),6) for z in np.roots([float(c) for c in sp.Poly(sp.expand(Q2c),r).all_coeffs()])
           if abs(np.imag(z))<1e-5 and 1.0<np.real(z)<2.0))[0]
print(f"E=1.2,a=0.9: Jc={Jc:.5f}  r_d={rd:.5f}")
# implicit dr_d/dE, dr_d/dJ from Q2=0 & Q2'=0 (double root moves with (E,J))
Q2r=sp.diff(Q2,r); Q2rr=sp.diff(Q2,r,2)
sub={E:E0,J:sp.Float(Jc),r:sp.Float(rd)}
# d/dE of [Q2=0]: Q2_r r' + Q2_E=0 on the double root; but r_d is double root so Q2_r=0 there.
# Use the pair {Q2=0,Q2_r=0}: differentiate both wrt E along Jc(E). Solve for r_d'(E), Jc'(E).
Q2_E=sp.diff(Q2,E).subs(sub); Q2_J=sp.diff(Q2,J).subs(sub)
Q2r_r=Q2rr.subs(sub); Q2r_E=sp.diff(Q2r,E).subs(sub); Q2r_J=sp.diff(Q2r,J).subs(sub)
# eq1 (Q2=0): Q2_J Jc' + Q2_E = 0  (Q2_r=0)  => Jc' = -Q2_E/Q2_J
Jcp=float(-Q2_E/Q2_J)
# eq2 (Q2_r=0): Q2r_r r' + Q2r_J Jc' + Q2r_E = 0 => r' = -(Q2r_J Jc' + Q2r_E)/Q2r_r
rdp_E=float(-(Q2r_J*Jcp+Q2r_E)/Q2r_r)     # dr_d/dE along the separatrix family
# but we want dr_d/dlambda along the DILATION ray: dE=-eps E, dJ=-eps J INDEPENDENTLY (not along Jc(E))
# partial derivs of r_d holding the OTHER charge: from {Q2=0,Q2_r=0} treat (E,J) both free:
# [Q2_r=0 & Q2=0] define r_d(E,J). d r_d: Q2_r=0 => Q2r_r dr_d + Q2r_E dE + Q2r_J dJ=0 (since Q2=0 auto at double root? no)
# proper: r_d(E,J) solves Q2_r(r_d;E,J)=0 AND Q2(r_d;E,J)=0 simultaneously only on the separatrix.
# Off the separatrix Q2 has no double root. So r_d is defined only WITH Jc(E). The dilation ray leaves
# the separatrix (E/A,J/A generally not on Jc). => the double root DISSOLVES under generic dilation.
print(f"  dJc/dE={Jcp:+.5f}   dr_d/dE (along separatrix family)={rdp_E:+.5f}")
# check: is the dilation ray (E/A,J/A) ON the separatrix? i.e. is Jc(E) homogeneous-compatible?
for A in [1.0,1.05,1.1]:
    Ea=float(E0)/A; Ja=Jc/A
    resA=sp.resultant(Q2,sp.diff(Q2,r),r).subs(E,sp.Float(Ea))
    Jcs=[complex(s).real for s in sp.solve(resA,J) if abs(complex(s).imag)<1e-9 and 2<complex(s).real<5]
    JcA=min((abs(jj-Ja),jj) for jj in Jcs)[1] if Jcs else float('nan')
    print(f"  A={A}: dilation gives J_eff={Ja:.5f}; separatrix at E_eff would be Jc={JcA:.5f}; on-sep? {abs(Ja-JcA)<1e-3}")
print("\n=> se il raggio di dilatazione NON resta sulla separatrice, il doppio root si DISSOLVE:")
print("   l'orbita ATTRAVERSA la separatrice sotto l'evoluzione adiabatica (separatrix crossing).")
print("   => non un controtermine dr_d/dl, ma un problema di strato limite (Neishtadt).")
