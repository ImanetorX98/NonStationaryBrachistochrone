import mpmath as mp
mp.iv.dps=18; mp.mp.dps=18
import cap_full as cf
from mpmath import iv
from scipy.optimize import brentq
import time

def moments_g(x,V0,bb,N,g=2.0):
    I=iv.mpf(0);Ip=iv.mpf(0);Ipp=iv.mpf(0);prev=mp.mpf(0)
    for i in range(1,N+1):
        s=mp.mpf(i)/N; un=s**g
        u=iv.mpf([float(prev),float(un)]);w=iv.mpf(float(un)-float(prev));prev=un
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
            if cert(r0_iv,a,b,bb,250): nS+=1
            else: m=(a+b)/2;stack.insert(0,(m,b));stack.insert(0,(a,m));continue
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

GRID=[(1.2,8),(1.6,8),(2.5,8),(1.2,12),(1.6,12),(2.5,12)]
with open('/tmp/cap_grid.log','w',buffering=1) as lg:
    lg.write("E  r0  bb  rquarter  RESULT  nc S M  time\n")
    for E,r0 in GRID:
        bb=E*E-1.0; rq=rquarter(bb,r0)
        ok,nc,nS,nM,dt=prove(r0,bb,2.6,min(rq,r0-0.1))
        line=f"E={E} r0={r0} bb={bb:.3f} rquarter={rq:.3f} : {'PASS' if ok else 'FAIL'} (nc={nc} S={nS} M={nM} {dt:.0f}s)"
        print(line,flush=True); lg.write(line+"\n")
    print("GRID DONE",flush=True); lg.write("GRID DONE\n")
