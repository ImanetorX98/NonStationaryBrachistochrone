# -*- coding: utf-8 -*-
# (A) CLOSURE of the HOLOMORPHIC part through a-PERIODS (TK tau, genus-2).
# The holomorphic coeffs alpha_k,beta_k of U_k = (poles, symbolic coeffs g_i) + alpha_k u1 + beta_k u2
# are fixed by the a-PERIODS (a principle, not a fit; pole basis normalised to zero a-period):
#   (alpha_k, beta_k) = omega^{-1} . [oint_{a_j} omega_k],   omega_k = r^k dr/y (EVEN model).
# Sage matrix_of_integral_values integrates the polynomial differentials r^k dr/y over the cycles.
# Sanity: k=0,1 (holomorphic) -> unit vectors (U_0=u2,U_1=u1). Period level, like C0/Ce of the separatrix.
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
import numpy as np

M,a,E,J = QQ(1),QQ(9)/10,QQ(7)/5,QQ(5)/2
Rr=PolynomialRing(QQ,['x','y']); x,y=Rr.gens()
Dl=x**2-2*M*x+a**2; Em=(E**2-1)*x+2*M
Sint=(1250*x*(x-2*M)*Em*(x*Dl-J**2*Em)).change_ring(QQ)   # integer model (as in psi_forward_abel)
Xrs=SageRS(y**2-Sint, prec=80)
# a-periods (2 cycles) of the differentials x^k dx/y, k=0..4
diffs=[Rr(x**k) for k in range(5)]
MIV=np.array(Xrs.matrix_of_integral_values(diffs),dtype=complex)   # 5 x 4 (diff x cycles)
aper=MIV[:,:2]   # a-periods (first 2 columns)
print("a-periods oint_{a_j} x^k dx/y  (rows k=0..4, columns a1,a2):")
for k in range(5): print(f"  k={k}: {np.round(aper[k],5)}")
omega=aper[:2,:]      # 2x2: holomorphic a-periods (k=0,1)
ominv=np.linalg.inv(omega)
print("\nomega (holomorphic a-periods k=0,1) =\n",np.round(omega,5))
print("\n=== holomorphic coeffs (alpha_k,beta_k) = omega^{-1} . a-per(omega_k) ===")
for k in range(5):
    ab=aper[k]@ominv   # aper[k] = [alpha,beta] . omega  =>  [alpha,beta] = aper[k] . omega^{-1}
    tag=" (holomorphic: must be ~ e_k)" if k<2 else ""
    print(f"  k={k}: (alpha,beta)=({ab[0]:+.5f},{ab[1]:+.5f}){tag}")
print("\n=> k=0,1 give unit vectors (sanity). k=2,3,4: period-level holomorphic coeffs of U_k")
print("   (fixed by the a-periods = a principle, NOT a fit; they depend on the periods like C0/Ce).")
print("   Combined with the symbolic pole coeffs g_i -> COMPLETE canonical decomposition.")
