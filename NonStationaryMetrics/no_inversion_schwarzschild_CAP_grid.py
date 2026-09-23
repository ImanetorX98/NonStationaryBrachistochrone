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

def rquarter(bb,r0):
    Vf=lambda r: r*(r*r-2*r)/(bb*r+2); V0=Vf(r0)
    return brentq(lambda r: Vf(r)-V0/4.0, 2.001, r0)

def weld(bb,r0,rq,width=1e-10,N=1500):
    """Close the joint between the certified window and the grazing bound.

    prove() certifies r_min in [2.6, rq], and the elementary quarter-potential
    bound covers V(r) >= V(r0)/4.  But rq comes from brentq in binary64, so on
    its own it leaves the true quarter radius unlocated: if r_{1/4} were a hair
    above rq, the sliver between them would be covered by neither argument.

    Two steps close it, and both are interval arithmetic:
      (a) certify the bridge cell [rq, rq+width] like any other cell;
      (b) enclose V(rq+width) - V(r0)/4 and check the enclosure is strictly
          positive, which places rq+width ABOVE the true quarter radius, so the
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
    vgap=cf.Vval(iv.mpf(upper),bb)-cf.Vval(iv.mpf(r0),bb)/4
    return dict(rq=rq,upper=upper,bridge_tag=cert(iv.mpf(r0),rq,upper,bb,N),
                # vgap.a / vgap.b are ivmpf, so str()/nstr() on them still prints
                # '[x, x]'; converting to mp.mpf gives the endpoint itself.
                inequality_enclosure=[mp.nstr(mp.mpf(vgap.a),12),
                                      mp.nstr(mp.mpf(vgap.b),12)],
                above_quarter=bool(vgap.a>0))

GRID=[(1.2,8),(1.6,8),(2.5,8),(1.2,12),(1.6,12),(2.5,12)]
_failed = False
with open('/tmp/cap_grid.log','w',buffering=1) as lg:
    lg.write("E  r0  bb  rquarter  RESULT  nc S M  time\n")
    for E,r0 in GRID:
        bb=E*E-1.0; rq=rquarter(bb,r0)
        # The window is capped at r0-0.1 as well: if rquarter ever exceeded it the
        # certified window would stop short of the weld, so say when that happens
        # instead of letting the joint open silently.
        top=min(rq,r0-0.1)
        ok,nc,nS,nM,dt=prove(r0,bb,2.6,top)
        if not ok: _failed = True
        line=f"E={E} r0={r0} bb={bb:.3f} rquarter={rq:.3f} : {'PASS' if ok else 'FAIL'} (nc={nc} S={nS} M={nM} {dt:.0f}s)"
        print(line,flush=True); lg.write(line+"\n")
        if top<rq:
            _failed=True
            w=f"  WELD GAP: window capped at r0-0.1={top:.6f} below rquarter={rq:.6f}"
            print(w,flush=True); lg.write(w+"\n"); continue
        wd=weld(bb,r0,rq)
        good=(wd['bridge_tag'] is not None) and wd['above_quarter']
        if not good: _failed = True
        w=(f"  weld [{wd['rq']:.9f},{wd['upper']:.9f}] tag={wd['bridge_tag']} ; "
           f"V(upper)-V(r0)/4 in [{wd['inequality_enclosure'][0]}, "
           f"{wd['inequality_enclosure'][1]}] > 0 ? {wd['above_quarter']} "
           f": {'OK' if good else 'FAIL'}")
        print(w,flush=True); lg.write(w+"\n")
    print("GRID DONE",flush=True); lg.write("GRID DONE\n")
# A runner that prints FAIL and exits 0 misleads any automation reading only the
# exit code -- and the manifest is read that way.
import sys as _sys
_sys.exit(1 if _failed else 0)
