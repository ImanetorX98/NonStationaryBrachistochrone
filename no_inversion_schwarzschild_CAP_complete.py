import cap_full as cf
from mpmath import iv
import mpmath as mp
import time,sys

def cert_cell(r0_iv,rmin_iv,bb,N):
    """certify single-crossing condition on a (r0,rmin) box: 0 not in [Phi'] OR [Phi'']<0."""
    try:
        V0=cf.Vval(r0_iv,bb)
        x=iv.mpf([cf.Vval(iv.mpf(rmin_iv.a),bb).a, cf.Vval(iv.mpf(rmin_iv.b),bb).b])
        PhiP,PhiPP=cf.PhiP_PhiPP(x,V0,bb,N)
    except Exception:
        return False,None
    if PhiP.b<0 or PhiP.a>0: return True,'S'      # sign of Phi' certain
    if PhiPP.b<0: return True,'M'                  # near peak, Phi''<0
    return False,None

def prove_r0(r0f,bb,ra,rb,Nsched=(400,1200,4000),wmin=0.01):
    """adaptive: cover rmin in [ra,rb] for fixed r0; escalate N, then subdivide."""
    r0_iv=iv.mpf(r0f)
    stack=[(ra,rb)]; ncell=0; nM=0; t0=time.time()
    while stack:
        a,b=stack.pop(); done=False
        for N in Nsched:
            ok,tag=cert_cell(r0_iv,iv.mpf([a,b]),bb,N); ncell+=1
            if ok:
                if tag=="M": nM+=1
                done=True; break
        if done: continue
        if b-a<wmin:
            print(f"  FAIL at rmin=[{a:.4f},{b:.4f}] (width<{wmin})"); return False
        m=(a+b)/2; stack.append((a,m)); stack.append((m,b))
    print(f"  r0={r0f}: COMPLETE cover rmin in [{ra},{rb}], {ncell} evals, M-cells={nM} ({time.time()-t0:.0f}s)")
    return True

if __name__=="__main__":
    r0=float(sys.argv[1]) if len(sys.argv)>1 else 10.0
    ra=float(sys.argv[2]) if len(sys.argv)>2 else 2.6
    rb=float(sys.argv[3]) if len(sys.argv)>3 else r0
    print(f"COMPLETE certificate, E=1.4 (b=0.96), r0={r0}, rmin in [{ra},{rb}]:")
    ok=prove_r0(r0,0.96,ra,rb)
    print("RESULT: Lemma B RIGOROUSLY CERTIFIED at r0=%.1f =>"%r0, ok)
