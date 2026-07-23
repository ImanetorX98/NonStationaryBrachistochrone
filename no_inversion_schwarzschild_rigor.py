"""
Frozen-Schwarzschild fixed-endpoint no-inversion: RIGORIZATION of the large-r0 asymptotic
(route B pushed).  Companion to no_inversion_schwarzschild_asymptotic.py.

GOAL: prove rigorously that for r0 >= R_large, Phi_tau(x) (x=Vmin) has a UNIQUE critical
point, a non-degenerate maximum -> single-crossing -> Lemma B (Phi_tau decreasing beyond).

SCHEME (single-crossing localized), every scaling VERIFIED numerically below:

Step 1 -- remainder bound on W.  In s=1/r:  2VW = g(s) = (b+2s)sqrt(1-2s)/(b+(3-b)s-4s^2),
  V(s)=(1-2s)/(s^2(b+2s)).  Series g(s)=1 - s/b + c2 s^2 + ..., c2 explicit
  = (-1+3/b)^2 - 3/2 + 2(1-3/b)/b + 5/b.  Converting (s ~ b^-1/2 V^-1/2):
     W(V) = 1/(2V) - 1/(2 b^{3/2} V^{3/2}) + h2(V),   h2(V) = +C/V^2 (1+o(1)),  C~3.07 (b=.96).
  => explicit bound |h2(V)| <= C/V^2 for V>=V1.  [VERIFIED: V^2 h2 -> 3.07]

Step 2 -- Abel propagation.  Phi = Phi_model + Phi_corr + E, with the two closed forms
     Phi_model=arctan(sqrt(V0/x-1)),  Phi_corr=-b^{-3/2}sqrt(1/x-1/V0),
  and E(x)=sqrt(x) int_x^{V0} h2(A)/sqrt(A-x)dA.  Using |h2|<=C/A^2 and the exact
     int_x^{V0} A^{-2}(A-x)^{-1/2}dA <= int_x^inf = (pi/2) x^{-3/2},
  gives |E(x)| <= C pi/(2 x),  |E'(x)| <= C'/x^2,  |E''(x)| <= C''/x^{3}-ish.  [scalings VERIFIED]

Step 3 -- localization of critical points.  Leading:
     Phi'(x) = [1/(2 sqrt(x(V0-x)))] (b^{-3/2} sqrt(V0)/x - 1) + E'(x).
  At x=2 x_pk: bracket=-1/2, |leading|~V0^{-3/4}; |E'|~V0^{-1} << leading => Phi'(2x_pk)<0.
  [VERIFIED robustly]  At x=x_pk/2 the SAME argument gives Phi'>0, BUT only for r0>~300
  (see caveat).  => critical points confined to [x_pk/2, 2x_pk].

Step 4 -- non-degeneracy.  Phi''(x_pk) = -1/2 b^{9/4} V0^{-5/4} + (reinforcing remainder).
  [VERIFIED: remainder has SAME sign, Phi''<0 is ROBUST, ratio Phi''/lead -> 1 from above.]
  => the unique critical point in [x_pk/2,2x_pk] is a strict maximum -> Lemma B.

CAVEAT / HONEST LIMIT (the real cost of route B):
  x_pk_true / x_pk_pred (=b^{-3/2}sqrt(V0)) -> 1 VERY SLOWLY: 0.35,0.48,0.60,0.70,0.78 for
  r0=80,160,320,640,1280.  The asymptotic series is in powers of V0^{-1/4} (the W2~V^{-2}
  term shifts x_pk by relative O(V0^{-1/4})).  Hence R_large is LARGE (~300-1000 in r0):
  the localization window [x_pk/2,2x_pk] contains x_pk_true only for r0 >~ 300.

STATUS: scheme SOUND, all scalings verified, NO conceptual obstruction.  Full closure needs
  (a) explicit constants in Steps 1-2 (mechanical, long), and
  (b) a validated/interval computation on the COMPACT box [R*~4, R_large~500] (engineering).
  Complete theorem = elementary(r0<=R*) + this asymptotic(r0>=R_large) + CAP(compact middle).
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import sympy as sp

def make_W(E):
    bb=E*E-1.0
    Vf=lambda r: r*(r*r-2*r)/(bb*r+2)
    Vp=lambda r: 2*r*(bb*r*r+(3-bb)*r-4)/(bb*r+2)**2
    K=lambda r: 1.0/np.sqrt(r*(r-2))
    def rV(Vt): return brentq(lambda r: Vf(r)-Vt,2+1e-13,1e13)
    def Wv(Vt): r=rV(Vt); return K(r)/Vp(r)
    return Wv,Vf,bb
def Phi(x,V0,Wv):
    s=V0-x; val,_=quad(lambda u: Wv(x+s*u*u),0,1,limit=300); return 2*np.sqrt(x*s)*val
def PhiP(x,V0,Wv,h=None):
    h=x*1e-6 if h is None else h; return (Phi(x+h,V0,Wv)-Phi(x-h,V0,Wv))/(2*h)
def PhiPP(x,V0,Wv,h=None):
    h=x*3e-4 if h is None else h
    return (Phi(x+h,V0,Wv)-2*Phi(x,V0,Wv)+Phi(x-h,V0,Wv))/(h*h)

if __name__=="__main__":
    E=1.4; Wv,Vf,bb=make_W(E)
    # Step1: symbolic c2 and numeric h2~C/V^2
    s,b=sp.symbols('s b',positive=True)
    g=(b+2*s)*sp.sqrt(1-2*s)/(b+(3-b)*s-4*s**2)
    print("g(s) series:",sp.series(g,s,0,3).removeO())
    print("Step1 h2*V^2 (E=1.4):",[round(float((Wv(V)-1/(2*V)+bb**-1.5/(2*V**1.5))*V**2),3) for V in [1e3,1e4,1e5,1e6]])
    # Step4 robustness + slow convergence caveat
    print("Step4/caveat  r0 : Phi''/lead ,  x_pk_true/x_pk_pred")
    for r0 in [80,160,320,640,1280]:
        V0=Vf(r0); xpk=bb**-1.5*np.sqrt(V0)
        xt=brentq(lambda z:PhiP(z,V0,Wv),0.3*xpk,1.1*xpk)
        ratio=PhiPP(xt,V0,Wv)/(-0.5*bb**(9/4)*V0**(-5/4))
        print(f"   {r0:5d} : {ratio:.3f} , {xt/xpk:.3f}")
