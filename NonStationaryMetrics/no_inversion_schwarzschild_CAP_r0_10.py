"""Computer-assisted certificate (interval arithmetic) for the fixed-endpoint no-inversion statement in
frozen Schwarzschild: r0 = 10 M, E = 1.4 (b = E^2 - 1 = 0.96), M = 1.

With the interval library cap_full.py (mpmath, 18 digits, outward rounding) the window r_min in [2.6, 6.05] M
is covered by adaptive cells.  On each cell the enclosures, taken in the potential variable V_min = V(r_min)
(V increasing in r), certify the single-crossing condition of the tau half-angle map: either 0 is not in the
enclosure of Phi' (sign certified, tag S) or the enclosure of Phi'' is negative (near the peak, tag M).
Prints COMPLETE, or FAIL with the first cell that cannot be certified; the progress log goes to
$NSM_OUTDIR/cap10_progress.log (default: the system temporary directory).

This run stops at r_min = 6.05 M, below the half-potential radius where the elementary grazing bound takes
over.  The entry (E, r0) = (1.4, 10) of no_inversion_schwarzschild_CAP_grid.py redoes it up to that radius,
closes the joint, and supersedes it.

Run: python3 no_inversion_schwarzschild_CAP_r0_10.py
"""
import mpmath as mp
mp.iv.dps=18; mp.mp.dps=18
import cap_full as cf
from mpmath import iv
import time
from tqdm import tqdm

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
        V0=cf.Vval(r0_iv,bb)
        x=iv.mpf([cf.Vval(iv.mpf(a),bb).a,cf.Vval(iv.mpf(b),bb).b])
        I,Ip,Ipp=moments_g(x,V0,bb,N); P1=V0-2*x; S=iv.sqrt(x*(V0-x))
        P=(P1/S)*I+2*S*Ip; PP=(-2/S-P1*P1/(2*S**3))*I+2*(P1/S)*Ip+2*S*Ipp
    except Exception: return None
    if P.b<0 or P.a>0: return 'S'
    if PP.b<0: return 'M'
    return None

def prove(r0f,bb,ra,rb,wthin=0.03,Ncap=(250,600,1500,5000),wmin=0.0015,logpath=None):
    r0_iv=iv.mpf(r0f); stack=[(ra,rb)]; nc=0;nS=0;nM=0;t0=time.time()
    lg=open(logpath,'w',buffering=1) if logpath else None
    bar=tqdm(total=round(rb-ra,6),desc=f"CAP r0={r0f}",ncols=100,
             bar_format="{l_bar}{bar}| {n:.3f}/{total:.3f} [{elapsed}] {postfix}")
    while stack:
        a,b=stack.pop(0); w=b-a
        if w>wthin:
            nc+=1
            if cert(r0_iv,a,b,bb,250): nS+=1;bar.update(w);bar.set_postfix_str(f"S={nS}M={nM}@{a:.3f}");
            else:
                m=(a+b)/2;stack.insert(0,(m,b));stack.insert(0,(a,m));continue
            if lg:lg.write(f"S wide[{a:.4f},{b:.4f}]\n")
            continue
        tag=None;uN=None
        for N in Ncap:
            nc+=1;tag=cert(r0_iv,a,b,bb,N);uN=N
            if tag:break
        if tag:
            if tag=='S':nS+=1
            else:nM+=1
            bar.update(w);bar.set_postfix_str(f"S={nS}M={nM}@{a:.3f} {tag}N{uN}")
            if lg:lg.write(f"{tag}[{a:.4f},{b:.4f}]N{uN}\n")
            continue
        if w<wmin:
            bar.close();m=f"FAIL[{a:.5f},{b:.5f}]w<{wmin}";print(m,flush=True)
            if lg:lg.write(m+"\n");lg.close()
            return False
        m=(a+b)/2;stack.insert(0,(m,b));stack.insert(0,(a,m))
    bar.close();msg=f"COMPLETE r0={r0f}[{ra},{rb}] nc={nc} S={nS} M={nM} ({time.time()-t0:.0f}s)"
    print(msg,flush=True)
    if lg:lg.write(msg+"\n");lg.close()
    return True

if __name__=="__main__":
    print("=== COMPLETE CAP r0=10 rmin(2.6,6.05) graded dps18 ===",flush=True)
    import os, tempfile
    ok=prove(10.0,0.96,2.6,6.05,logpath=os.path.join(os.environ.get('NSM_OUTDIR', tempfile.gettempdir()), 'cap10_progress.log'))
    print("RESULT:",ok,flush=True)
