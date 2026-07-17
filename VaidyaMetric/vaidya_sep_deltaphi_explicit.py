# MICROSTEP FINALE (numerico robusto): delta phi_tau|_sep ESPLICITO via riduzione POLE-ADAPTED.
# G = P(r)/((r-r_d)^2 sqrtQ) + sum c_k U_k,   sqrtS=(r-r_d)sqrtQ,  Q = S/((r-r_d)^2) (deg4).
# Identita' polinomiale: P'(r-r_d)Q - P(2Q + (1/2)(r-r_d)Q') + (sum c_k r^k)(r-r_d)^2 Q = N_m.
import numpy as np, mpmath as mp, sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq
np.set_printoptions(suppress=True)
E=1.4; m=1.0
# Jc esatto
r,J=sp.symbols('r J'); Es=sp.Rational(7,5)
Ssym=sp.expand(r*(r-2)*((Es**2-1)*r+2)*(r**2*(r-2)-J**2*((Es**2-1)*r+2)))
_Jsol=sp.solve(sp.Eq(sp.resultant(Ssym,sp.diff(Ssym,r),r),0),J)
Jc=float([s for s in _Jsol if s.is_real and float(s)>1][0]);
def P(*c): return np.array(c,dtype=float)              # helper
# --- polinomi (coeff high->low, numpy convention) ---
def poly_S(Jv):    # S(r) = r(r-2)DE(r^2(r-2)-J^2 DE), DE=(E^2-1)r+2
    DE=np.array([E**2-1,2.0]); 
    p=np.polymul(np.polymul([1,0],[1,-2]),DE)               # r(r-2)DE
    br=np.polysub(np.polymul([1,0,0],[1,-2]),Jv**2*np.polymul(DE,[1]))  # r^2(r-2)-J^2 DE
    return np.polymul(p,br)
Sc=poly_S(Jc)
# dm S, dm K:
def poly_S_m(mv,Jv):
    DE=np.array([E**2-1,2*mv]); p=np.polymul(np.polymul([1,0],[1,-2*mv]),DE)
    br=np.polysub(np.polymul([1,0,0],[1,-2*mv]),Jv**2*DE); return np.polymul(p,br)
h=1e-6; dmS=(poly_S_m(m+h,Jc)-poly_S_m(m-h,Jc))/(2*h)
dmS=dmS[-len(Sc):] if len(dmS)>=len(Sc) else np.r_[np.zeros(len(Sc)-len(dmS)),dmS]
K=np.polymul([Jc],[E**2-1,2.0])          # K=J*DE  -> [Jc(E^2-1), 2Jc]
dmK=np.array([2*Jc])                       # dm(J DE)=2J
# N_m = S*dmK - 1/2 K*dmS
Nm=np.polysub(np.polymul(Sc,dmK), 0.5*np.polymul(K,dmS))
# --- r_d (radice doppia) e Q=S/(r-r_d)^2 ---
rts=np.roots(Sc); 
pairs=[(i,j) for i in range(len(rts)) for j in range(i+1,len(rts)) if abs(rts[i]-rts[j])<1e-6]
rd=float(np.real((rts[pairs[0][0]]+rts[pairs[0][1]])/2))
Q,remQ=np.polydiv(Sc, np.polymul([1,-rd],[1,-rd]))
print(f"Jc={Jc:.10f}  r_d={rd:.8f}  |resto polydiv Q|={np.max(np.abs(remQ)):.1e}  deg Q={len(Q)-1}")
Qp=np.polyder(Q)
# --- costruzione sistema lineare identita' ---
pdeg=5; ncoef_c=2
def contrib_P(i):    # coeff di p_i (P=r^i): P'=i r^{i-1}; termine P'(r-rd)Q - P(2Q+0.5(r-rd)Q')
    Pi=np.zeros(pdeg+1); Pi[pdeg-i]=1.0                  # r^i
    Ppi=np.polyder(Pi)
    t1=np.polymul(np.polymul(Ppi,[1,-rd]),Q)
    t2=np.polymul(Pi, np.polyadd(2*Q, 0.5*np.polymul([1,-rd],Qp)))
    return np.polysub(t1,t2)
def contrib_c(k):    # coeff di c_k (r^k):  r^k (r-rd)^2 Q
    rk=np.zeros(k+1); rk[0]=1.0
    return np.polymul(np.polymul(rk,np.polymul([1,-rd],[1,-rd])),Q)
cols=[contrib_P(i) for i in range(pdeg+1)]+[contrib_c(k) for k in range(ncoef_c)]
maxdeg=max(len(c) for c in cols); maxdeg=max(maxdeg,len(Nm))
def padL(a,n): return np.r_[np.zeros(n-len(a)),a]
Amat=np.array([padL(c,maxdeg) for c in cols]).T
bvec=padL(Nm,maxdeg)
xsol,res,rk_,sv=np.linalg.lstsq(Amat,bvec,rcond=None)
resid=Amat@xsol-bvec
pcoef=xsol[:pdeg+1]; ccoef=xsol[pdeg+1:]
print(f"c_k={ccoef}  |residuo identita'|={np.max(np.abs(resid)):.2e}")
Ppoly=pcoef   # P(r)=sum pcoef r^i (high->low)
# --- verifica dG/dr = dm F ---
sQ=lambda x: np.sqrt(np.polyval(Q,x))
def Gpart_alg(x): return np.polyval(Ppoly,x)/((x-rd)**2*sQ(x))
# dm F numerico:
def dmF(x):
    def Fm(mv):
        Sc_=poly_S_m(mv,Jc); return np.polyval(K if False else np.polymul([Jc],[E**2-1,2*mv]),x)/np.sqrt(np.polyval(Sc_,x))
    return (Fm(m+1e-6)-Fm(m-1e-6))/2e-6
# U_k on separatrix:
r0=12.0
def U(x,k): return quad(lambda t:t**k/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def G(x): return Gpart_alg(x)-Gpart_alg(r0)+sum(ccoef[k]*U(x,k) for k in range(ncoef_c))
print("\nverifica dG/dr = dm F (differenza finita di G vs dmF):")
for xt in [11.0,10.0,9.2]:
    dG=(G(xt+1e-5)-G(xt-1e-5))/2e-5
    print(f"  r={xt:4.1f}  dG/dr={dG:+.8f}  dmF={dmF(xt):+.8f}  diff={abs(dG-dmF(xt)):.1e}")
