import mpmath as mp
mp.iv.dps=18; mp.mp.dps=18
import cap_full as cf
from mpmath import iv
from scipy.optimize import brentq
import time

def moments_g(x,V0,bb,N,g=2.0):
    # The panel width must be an ENCLOSURE of the exact difference between the
    # two float endpoints, not a float subtraction promoted to a point interval:
    # the latter rounds once in binary64 and then asserts the rounded value is
    # exact.  Subtracting the two point intervals is outward rounded and is a
    # rigorous enclosure, which can only widen the certificate, never narrow it.
    I=iv.mpf(0);Ip=iv.mpf(0);Ipp=iv.mpf(0);prev=mp.mpf(0)
    for i in range(1,N+1):
        s=mp.mpf(i)/N; un=s**g
        lo=iv.mpf(float(prev));hi=iv.mpf(float(un))
        u=iv.mpf([float(prev),float(un)]);w=hi-lo;prev=un
        A=x+(V0-x)*u*u; r=cf.r_of_V(A,bb); Wv,Wp,Wpp=cf.WVderivs(r,bb); om=1-u*u
        I+=w*Wv;Ip+=w*Wp*om;Ipp+=w*Wpp*om*om
    return I,Ip,Ipp
def cert(r0_iv,a,b,bb,N):
    try:
        V0=cf.Vval(r0_iv,bb); x=iv.mpf([cf.Vval(iv.mpf(a),bb).a,cf.Vval(iv.mpf(b),bb).b])
        I,Ip,Ipp=moments_g(x,V0,bb,N);P1=V0-2*x;S=iv.sqrt(x*(V0-x))
        P=(P1/S)*I+2*S*Ip; PP=(-2/S-P1*P1/(2*S**3))*I+2*(P1/S)*Ip+2*S*Ipp
    except Exception: return None
    if P.b<0 or P.a>0: return 'S'
    if PP.b<0: return 'M'
    return None
def prove(r0f,bb,ra,rb,wthin=0.03,Ncap=(250,600,1500,5000),wmin=0.0015):
    r0_iv=iv.mpf(r0f);stack=[(ra,rb)];nc=0;nS=0;nM=0;t0=time.time()
    while stack:
        a,b=stack.pop(0);w=b-a
        if w>wthin:
            nc+=1
            # cert() returns 'S' or 'M', both truthy: the previous form counted
            # every wide-cell success as 'S', so the S/M split was wrong even
            # though the certified disjunction, and hence PASS/FAIL, was not.
            tag=cert(r0_iv,a,b,bb,250)
            if tag:
                nS+=(tag=='S'); nM+=(tag=='M')
            else:
                m=(a+b)/2;stack.insert(0,(m,b));stack.insert(0,(a,m))
            continue
        tag=None
        for N in Ncap:
            nc+=1;tag=cert(r0_iv,a,b,bb,N)
            if tag:break
        if tag:
            nS+=(tag=='S');nM+=(tag=='M');continue
        if w<wmin: return (False,nc,nS,nM,time.time()-t0)
        m=(a+b)/2;stack.insert(0,(m,b));stack.insert(0,(a,m))
    return (True,nc,nS,nM,time.time()-t0)

# The elementary bound that takes over above the certified window is the
# GRAZING one, V(r_min) >= V(r0)/2: there both terms of the closed form for
# sqrt(x) Phi' are <= 0 (no_inversion_schwarzschild_closedform.py).  An earlier
# version welded at the quarter radius V(r0)/4, where the closed form alone does
# not settle the sign (the boundary term is positive there), so the window now
# runs up to the half-potential radius.
FRAC=2.0
def rhalf(bb,r0):
    Vf=lambda r: r*(r*r-2*r)/(bb*r+2); V0=Vf(r0)
    return brentq(lambda r: Vf(r)-V0/FRAC, 2.001, r0)

def weld(bb,r0,rq,width=1e-10,N=1500):
    """Close the joint between the certified window and the grazing bound.

    prove() certifies r_min in [2.6, rq], and the elementary grazing bound
    covers V(r) >= V(r0)/2.  But rq comes from brentq in binary64, so on
    its own it leaves the true half radius unlocated: if r_{1/2} were a hair
    above rq, the sliver between them would be covered by neither argument.

    Two steps close it, and both are interval arithmetic:
      (a) certify the bridge cell [rq, rq+width] like any other cell;
      (b) enclose V(rq+width) - V(r0)/2 and check the enclosure is strictly
          positive, which places rq+width ABOVE the true half radius, so the
          grazing bound applies from there on.
    Together the two windows overlap and the interval is covered with no gap.

    V is strictly increasing on r >= 2 for every b > 0, so "above in V" is
    "above in r": with r = 2+s, V' = 2rN/(br+2)^2 and
        N(2+s) = b s^2 + 3(b+1) s + 2(b+1),
    every coefficient positive and the constant term 2(b+1) > 0.  That is what
    makes cap_full.r_of_V's two-sided check a genuine enclosure of the preimage
    rather than of one branch of it.
    """
    upper=rq+width
    vgap=cf.Vval(iv.mpf(upper),bb)-cf.Vval(iv.mpf(r0),bb)/int(FRAC)
    return dict(rq=rq,upper=upper,bridge_tag=cert(iv.mpf(r0),rq,upper,bb,N),
                # vgap.a / vgap.b are ivmpf, so str()/nstr() on them still prints
                # '[x, x]'; converting to mp.mpf gives the endpoint itself.
                inequality_enclosure=[mp.nstr(mp.mpf(vgap.a),12),
                                      mp.nstr(mp.mpf(vgap.b),12)],
                above_half=bool(vgap.a>0))

# (E, r0, b=E^2-1).  The r0=10 entry is the configuration of
# no_inversion_schwarzschild_CAP_r0_10.py, at b=0.96 as there; that run stopped at
# r_min=6.05, below the half radius, so it is redone here with the weld.
GRID=[(E,r0,E*E-1.0) for E,r0 in [(1.2,8),(1.6,8),(2.5,8),(1.2,12),(1.6,12),(2.5,12)]]
GRID.append((1.4,10,0.96))
_failed = False
with open('/tmp/cap_grid.log','w',buffering=1) as lg:
    lg.write("E  r0  bb  rhalf  RESULT  nc S M  time\n")
    for E,r0,bb in GRID:
        rq=rhalf(bb,r0)
        # The window is capped at r0-0.1 as well: if rhalf ever exceeded it the
        # certified window would stop short of the weld, so say when that happens
        # instead of letting the joint open silently.
        top=min(rq,r0-0.1)
        ok,nc,nS,nM,dt=prove(r0,bb,2.6,top)
        if not ok: _failed = True
        line=f"E={E} r0={r0} bb={bb:.3f} rhalf={rq:.3f} : {'PASS' if ok else 'FAIL'} (nc={nc} S={nS} M={nM} {dt:.0f}s)"
        print(line,flush=True); lg.write(line+"\n")
        if top<rq:
            _failed=True
            w=f"  WELD GAP: window capped at r0-0.1={top:.6f} below rhalf={rq:.6f}"
            print(w,flush=True); lg.write(w+"\n"); continue
        wd=weld(bb,r0,rq)
        good=(wd['bridge_tag'] is not None) and wd['above_half']
        if not good: _failed = True
        w=(f"  weld [{wd['rq']:.9f},{wd['upper']:.9f}] tag={wd['bridge_tag']} ; "
           f"V(upper)-V(r0)/2 in [{wd['inequality_enclosure'][0]}, "
           f"{wd['inequality_enclosure'][1]}] > 0 ? {wd['above_half']} "
           f": {'OK' if good else 'FAIL'}")
        print(w,flush=True); lg.write(w+"\n")
    print("GRID DONE",flush=True); lg.write("GRID DONE\n")
# A runner that prints FAIL and exits 0 misleads any automation reading only the
# exit code -- and the manifest is read that way.
import sys as _sys
_sys.exit(1 if _failed else 0)
