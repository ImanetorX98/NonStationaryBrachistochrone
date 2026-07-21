import numpy as np
from scipy.integrate import quad
M=1.0
def setup(E):
    b=E*E-1.0
    Delta=lambda r: r*r-2*M*r
    DE=lambda r: b*r+2*M
    V=lambda r: r*Delta(r)/DE(r)
    # V'(r) exact: 2r(br^2+(3-b)r-4)/(br+2)^2  (M=1)
    Vp=lambda r: 2*r*(b*r*r+(3-b)*r-4)/(b*r+2)**2
    K=lambda r: 1.0/np.sqrt(r*(r-2*M))
    return b,Delta,DE,V,Vp,K

def Phi_direct(rm,E,r0=10.0):
    # original angle integral, tau branch, a=0
    b,Delta,DE,V,Vp,K=setup(E)
    J=np.sqrt(rm*Delta(rm)/DE(rm))
    Vmin=V(rm)
    def integ(s):
        r=rm+s*s
        P=r*Delta(r)-J*J*DE(r)   # = DE*(V-Vmin)
        return 2*J*np.sqrt(r*(r-2*M)*DE(r))/(Delta(r)*np.sqrt(P/(r-rm)))
    val,_=quad(integ,0,np.sqrt(r0-rm),limit=200)
    return val

def Phi_clean(rm,E,r0=10.0):
    b,Delta,DE,V,Vp,K=setup(E)
    Vmin=V(rm); V0=V(r0); Sig=np.sqrt(V0-Vmin)
    # W(V)=K(r)/V'(r) as fn of V; need r(V). invert numerically via V monotone
    from scipy.optimize import brentq
    def rofV(Vt):
        return brentq(lambda r: V(r)-Vt, 2.0+1e-9, 1e4)
    def Wof(Vt):
        r=rofV(Vt); return K(r)/Vp(r)
    val,_=quad(lambda s: Wof(Vmin+s*s),0,Sig,limit=200)
    return 2*np.sqrt(Vmin)*val

for E in [1.2,1.4,2.0]:
    print(f"E={E}")
    for rm in [3.0,4.0,5.0,6.0,7.0,8.0,9.0]:
        pd=Phi_direct(rm,E); pc=Phi_clean(rm,E)
        print(f"  rm={rm}: Phi_direct={pd:.6f}  Phi_clean={pc:.6f}  diff={pd-pc:.2e}")

# ---------------------------------------------------------------------------
# FROZEN-CASE (Vaidya mu=0 / Thakurta a=0) closed reduction, verified above:
#
#   Phi_tau(Vmin) = sqrt(Vmin) int_{Vmin}^{V0} W(V)/sqrt(V-Vmin) dV
#                 = 2 sqrt(Vmin(V0-Vmin)) * I(Vmin),
#   I(Vmin) = int_0^1 W(Vmin+(V0-Vmin)u^2) du,   W(V)=K/V' (fixed, >0, decreasing, convex),
#   V0=V(r0) FIXED (fixed endpoint).  r(V) is an explicit cubic root
#   (r^3-2r^2-b V r-2V=0), so W(V) is elementary-algebraic.
#
# Monotonicity dPhi/dVmin<0  <=>  the crisp inequality
#   (STAR)   -I'(Vmin) > (V0-2Vmin)/(2 Vmin (V0-Vmin)) * I(Vmin),
#            -I' = int_0^1 (-W'(A))(1-u^2) du,   A=Vmin+(V0-Vmin)u^2.
#
# FINDINGS (numerics in schw_ineq.py, check_grazing.py):
#  * (STAR) holds on the whole scattering regime r_min>r_pk, but is TIGHT as
#    r_min->r_pk+ (margin->0): r_pk is DEFINED by dPhi=0, i.e. (STAR) with equality.
#    => the threshold is genuinely transcendental; NO elementary uniform certificate
#    can prove (STAR) for all r_min>r_pk.  A Chebyshev bound (2/3)int(-W')du is too weak.
#  * ELEMENTARY sub-result: for Vmin>=V0/2 (grazing) both pieces of Phi'/Phi are <=0,
#    so Phi strictly decreasing -- but the V0/2 radius (~7.5M for r0=10) is ABOVE the
#    physical turning radii (~4.7-6.4M), so grazing alone does not cover them.
#  * The physical fixed-endpoint turning points sit in Vmin << V0/2, strictly between
#    r_pk and the grazing radius: (STAR) holds there with positive but transcendental margin.
