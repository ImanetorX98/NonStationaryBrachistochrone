"""
Frozen-Schwarzschild fixed-endpoint no-inversion: LARGE-r0 ASYMPTOTIC (transcendental route).

The tight point r_pk (max of Phi_tau) GROWS without bound as r0->inf (numerically
r_pk ~ r0^0.55, V_pk ~ b^-3/2 sqrt(V0)); hence near r_pk we are in the LARGE-V regime
where W(V) has a clean expansion.  This makes the neighbourhood of the max a controlled
perturbation of an EXACTLY solvable model.

Exact identity (verified symbolically):
    2 V W(V) = DE*sqrt(Delta)/N = (br+2) sqrt(r^2-2r) / (br^2+(3-b)r-4),  V=V(r).
Large-r expansion (verified): 2VW = 1 - (1/b)(1/r) + O(1/r^2),  and r ~ sqrt(bV), so
    W(V) = 1/(2V) - 1/(2 b^{3/2} V^{3/2}) + O(V^-2).

MODEL  W0=1/(2V)  gives the Abel integral in closed form:
    Phi_model(x) = sqrt(x) int_x^{V0} (1/(2A))/sqrt(A-x) dA = arctan( sqrt(V0/x - 1) ),
    Phi_model'(x) = -1/(2 sqrt(x(V0-x)))  < 0  everywhere  (monotone: Lemma B holds at
                     leading order -- the geometric-optics angle is a pure arctan).

CORRECTION  W1=-1/(2 b^{3/2} V^{3/2})  also integrates in closed form:
    Phi_corr(x) = -b^{-3/2} sqrt(1/x - 1/V0),
    Phi_corr'(x) = + b^{-3/2} sqrt(V0) / (2 x^{3/2} sqrt(V0-x))  > 0.

Sum:
    Phi'(x) ~ [ 1 / (2 sqrt(x(V0-x))) ] * ( b^{-3/2} sqrt(V0)/x - 1 ).
=> SINGLE sign change at   x_pk ~ b^{-3/2} sqrt(V0),   Phi'<0 for x>x_pk.
   The bracket crosses zero with slope -b^{-3/2}sqrt(V0)/x^2 < 0  =>  Phi''(x_pk) < 0
   (NON-DEGENERATE maximum, with margin).  Hence the single-crossing criterion
   ("every critical point has Phi''<0 => at most one critical point") is NON-TIGHT in
   the asymptotic regime: the tightness that blocked the elementary route disappears.

STATUS: leading order done in closed form; STRUCTURE (single non-degenerate max) verified.
Full rigor for r0>=R_large needs higher-order terms + a rigorous remainder bound to pin
x_pk (leading prediction is quantitatively ~2x too large at moderate r0, converges ~V^-1/2).
Complete theorem = [elementary x>=V0/4 for r0<=R*] + [this asymptotic for r0>=R_large]
+ [validated/interval computation on the compact middle R*<=r0<=R_large].
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import sympy as sp

def verify_symbolic():
    r,b=sp.symbols('r b',positive=True)
    Delta=r*r-2*r; DE=b*r+2; N=b*r*r+(3-b)*r-4
    twoVW=DE*sp.sqrt(Delta)/N
    lead=sp.limit(twoVW,r,sp.oo); c1=sp.simplify(sp.limit((twoVW-1)*r,r,sp.oo))
    print(f"[sym] 2VW -> {lead},  coeff of 1/r = {c1}  (expect 1, -1/b)")

def make_W(E):
    bb=E*E-1.0
    V=lambda rr: rr*(rr*rr-2*rr)/(bb*rr+2)
    Vp=lambda rr: 2*rr*(bb*rr*rr+(3-bb)*rr-4)/(bb*rr+2)**2
    K=lambda rr: 1.0/np.sqrt(rr*(rr-2))
    def rofV(Vt): return brentq(lambda rr: V(rr)-Vt,2+1e-13,1e12)
    def W(Vt): rr=rofV(Vt); return K(rr)/Vp(rr)
    return W,V,bb

def Phi_real(x,V0,W):
    s=V0-x; val,_=quad(lambda u: W(x+s*u*u),0,1,limit=200); return 2*np.sqrt(x*s)*val
def Phi_asym(x,V0,bb):
    return np.arctan(np.sqrt(V0/x-1)) - bb**-1.5*np.sqrt(1/x-1/V0)

if __name__=="__main__":
    verify_symbolic()
    print("\nPhi_real vs Phi_asym (arctan + V^-1/2 correction), E=1.4:")
    E=1.4; W,V,bb=make_W(E)
    for r0 in [40,160,640,2560]:
        V0=V(r0); xpk=bb**-1.5*np.sqrt(V0); x=1.5*xpk
        pr=Phi_real(x,V0,W); pa=Phi_asym(x,V0,bb)
        print(f"  r0={r0:5d} x=1.5 x_pk={x:9.2f}: Phi_real={pr:.5f} Phi_asym={pa:.5f} rel={abs(pr-pa)/pr:.2e}")
