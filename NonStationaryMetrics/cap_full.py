# -*- coding: utf-8 -*-
"""
cap_full -- shared interval-arithmetic library for the frozen-Schwarzschild
fixed-endpoint no-inversion computer-assisted proposition.

RESTORED for CQG-116884, referee major comment 10.  The CAP runners

    no_inversion_schwarzschild_CAP_r0_10.py
    no_inversion_schwarzschild_CAP_grid.py
    no_inversion_schwarzschild_CAP_complete.py

all begin with `import cap_full as cf`, but the module was never archived: it
had been factored out of no_inversion_schwarzschild_CAP.py during development
and existed only in the working tree.  A fresh clone therefore could not run the
representative computer-assisted proposition at all -- exactly the referee's
finding.  This file is that library, extracted verbatim from the prototype so
that the two cannot drift apart.

Everything here is mpmath INTERVAL arithmetic with outward rounding, so every
returned enclosure is guaranteed to contain the true value.

Exports used by the runners:
    Vval(r_iv, bb)          V(r) as an interval
    r_of_V(A_iv, bb)        rigorous inverse: an r-interval enclosing r(A)
    WVderivs(r_iv, bb)      (W, dW/dV, d2W/dV2) as intervals
    PhiP_PhiPP(x, V0, bb, N)  guaranteed enclosures of [Phi'] and [Phi'']

Parameter convention: bb = b = E^2 - 1 (so the demonstrated case E = 1.4 is
bb = 0.96), M = 1, Schwarzschild frozen at constant conformal factor.

Self-test:  python3 cap_full.py
"""
import mpmath as mp
from mpmath import iv

mp.iv.dps = 28
mp.mp.dps = 28


from mpmath import iv
import time
mp.iv.dps=28; mp.mp.dps=28

class J:
    __slots__=('v','d','dd')
    def __init__(s,v,d=None,dd=None):
        s.v=v; s.d=iv.mpf(0) if d is None else d; s.dd=iv.mpf(0) if dd is None else dd
    def __add__(a,b): b=C(b); return J(a.v+b.v,a.d+b.d,a.dd+b.dd)
    __radd__=__add__
    def __sub__(a,b): b=C(b); return J(a.v-b.v,a.d-b.d,a.dd-b.dd)
    def __rsub__(a,b): b=C(b); return J(b.v-a.v,b.d-a.d,b.dd-a.dd)
    def __mul__(a,b): b=C(b); return J(a.v*b.v,a.d*b.v+a.v*b.d,a.dd*b.v+2*a.d*b.d+a.v*b.dd)
    __rmul__=__mul__
    def recip(a): v=1/a.v; d=-a.d*v*v; dd=-a.dd*v*v+2*a.d*a.d*v*v*v; return J(v,d,dd)
    def __truediv__(a,b): return a*C(b).recip()
    def __rtruediv__(a,b): return C(b)*a.recip()
    def sqrt(a): v=iv.sqrt(a.v); d=a.d/(2*v); dd=a.dd/(2*v)-a.d*a.d/(4*v*a.v); return J(v,d,dd)
def C(b): return b if isinstance(b,J) else J(iv.mpf(b))

def Wjets(r_iv,bb):
    r=J(r_iv,iv.mpf(1),iv.mpf(0)); b=iv.mpf(bb)
    Delta=r*(r-2); DE=r*b+2; N=r*r*b+(3-bb)*r-4
    K=C(1)/Delta.sqrt(); Vp=2*r*N/(DE*DE); W=K/Vp; Vj=r*Delta/DE
    return W,Vj
def WVderivs(r_iv,bb):
    W,Vj=Wjets(r_iv,bb)
    return W.v, W.d/Vj.d, (W.dd*Vj.d-W.d*Vj.dd)/(Vj.d**3)
def Vval(r_iv,bb):
    b=iv.mpf(bb); return r_iv*(r_iv*(r_iv-2))/(r_iv*b+2)

def r_of_V(A_iv,bb):
    a=A_iv.a; bb_=A_iv.b
    ra=float(mp.findroot(lambda r: float(Vval(iv.mpf(r),bb).a)-float(a),3.0))
    rb=float(mp.findroot(lambda r: float(Vval(iv.mpf(r),bb).a)-float(bb_),max(3.0,ra)))
    for margin in [1e-6,1e-4,1e-2,1e-1,0.5,1.0,5.0]:
        rlo=max(2.0000001,ra-margin); rhi=rb+margin
        if Vval(iv.mpf(rlo),bb).b<=a and Vval(iv.mpf(rhi),bb).a>=bb_:
            return iv.mpf([rlo,rhi])
    raise RuntimeError("r_of_V enclosure failed")

def moments(x,V0,bb,N=400):
    I=iv.mpf(0); Ip=iv.mpf(0); Ipp=iv.mpf(0)
    for i in range(N):
        u=iv.mpf([mp.mpf(i)/N,mp.mpf(i+1)/N]); w=iv.mpf(1)/N
        A=x+(V0-x)*u*u; r=r_of_V(A,bb)
        Wv,Wp,Wpp=WVderivs(r,bb); om=1-u*u
        I+=w*Wv; Ip+=w*Wp*om; Ipp+=w*Wpp*om*om
    return I,Ip,Ipp

def PhiP_PhiPP(x,V0,bb,N=400):
    I,Ip,Ipp=moments(x,V0,bb,N)
    P1=V0-2*x; S=iv.sqrt(x*(V0-x))
    PhiP=(P1/S)*I+2*S*Ip
    PhiPP=(-2/S-P1*P1/(2*S**3))*I+2*(P1/S)*Ip+2*S*Ipp
    return PhiP,PhiPP

if __name__ == "__main__":
    # reproduce the value quoted in the manuscript: at r0 = 10, b = 0.96, the
    # enclosure of Phi'' near the peak is strictly negative.
    bb = 0.96
    V0 = Vval(iv.mpf(10.0), bb)
    x = Vval(iv.mpf(5.0), bb)
    P, PP = PhiP_PhiPP(x, V0, bb, N=200)
    print("cap_full self-test  (r0=10, b=0.96, r_min=5.0, N=200)")
    print(f"  [Phi'] =[{mp.nstr(P.a,6)}, {mp.nstr(P.b,6)}]   Phi'<0 ? {P.b<0}")
    print(f"  [Phi'']=[{mp.nstr(PP.a,6)}, {mp.nstr(PP.b,6)}]  Phi''<0? {PP.b<0}")
    assert P.b < 0, "Phi' should be certified negative at r_min=5.0"
    print("  OK")
