# Holomorphic components of omega_a=dE F and omega_b=d(eta), via a-periods (exact).
# b = (omega_A^T)^-1 . (a-periods of omega),  omega_A = a-periods of {1,x}/sqrt(S).
# a-per(omega_a)=sum c_k a-per(x^k/sqrtS);  a-per(omega_b)=a-per(x^3)-2M a-per(x^2).
# det = b^A_0 b^B_1 - b^A_1 b^B_0.  If ~0: no holomorphic x holomorphic weight-2 term.
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
import numpy as np

M = 1.0
Rr = PolynomialRing(QQ, ['x', 'y']); x, y = Rr.gens()
Sint = 1200*x**6-2300*x**5-11428*x**4-5519*x**3+24700*x**2+62500*x   # 1250*S, E=7/5
X = SageRS(y**2 - Sint, prec=80)
P = np.array(X.matrix_of_integral_values([Rr(1), x, x**2, x**3, x**4]), dtype=complex)
print('P shape', P.shape)                       # 5 x 4  (columns: a1,a2,b1,b2)
aper = P[:, :2]                                  # a-periods (5 x 2), rows=x^0..x^4

ck = np.array([-0.5314914946993489, 1.9793812403458497, -0.8122638904451002,
               -0.3604947564454957, 0.18857907778766347])   # analytic c_k at E=7/5
aper_a = (ck @ aper).astype(complex)
aper_b = aper[3] - 2*M*aper[2]                   # a-per(omega_b) = a-per(x^3)-2M a-per(x^2)
omA = aper[:2].T                                 # omega_A: omA[j,i]=a-per(x^i)_j; (a-per=omA^T b)
# aper_a (2,) = omA^T? -> aper_a[j]=sum_i b_i P[i,j]=sum_i b_i aper[i,j]; aper[:2,:] is (2diff,2cyc)
# omega_A[i,j]=aper[i,j] (i=diff 0,1 ; j=cyc). aper_a[j]=sum_i b_i omega_A[i,j] -> aper_a=b @ omega_A
omega_A = aper[:2]                               # (2,2): [i=diff, j=cyc]
bA = np.linalg.solve(omega_A.T.astype(complex), aper_a)
bB = np.linalg.solve(omega_A.T.astype(complex), aper_b.astype(complex))
det = bA[0]*bB[1] - bA[1]*bB[0]
print('b^A (holomorphic component of dE phi0) =', np.round(bA,6))
print('b^B (holomorphic component of the clock) =', np.round(bB,6))
print(f'det = b^A_0 b^B_1 - b^A_1 b^B_0 = {det:.4e}')
print(f'  |det|/(|bA||bB|) = {abs(det)/(np.linalg.norm(bA)*np.linalg.norm(bB)+1e-30):.2e}')
print('=> if ~0: NO holomorphic x holomorphic weight-2 term -> psi_zeta has weight 1 (zeta only).')
print('   if !=0: there is a second weight-2 term (Beilinson abelian integral).')
