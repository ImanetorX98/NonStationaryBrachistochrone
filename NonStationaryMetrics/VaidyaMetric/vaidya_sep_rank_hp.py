# Rango di A modulo peso-1 ad ALTA PRECISIONE (mpmath), per confermare = 5 esatto.
# A_perp = A - W (W^+ A);  rango = # autovalori Gram(A_perp) sopra soglia. dps alto -> definitivo.
import numpy as np, mpmath as mp, sympy as sp, sys, logging
from scipy.integrate import quad
try: from tqdm import tqdm
except ImportError:
    def tqdm(it,**k): return it
logging.basicConfig(level=logging.INFO,format='[%(asctime)s] %(message)s',datefmt='%H:%M:%S',stream=sys.stdout)
log=logging.getLogger(); mp.mp.dps=40
E=1.4; m=1.0; r0=12.0
r,J=sp.symbols('r J'); Es=sp.Rational(7,5)
Ssym=sp.expand(r*(r-2)*((Es**2-1)*r+2)*(r**2*(r-2)-J**2*((Es**2-1)*r+2)))
Jc=float([s for s in sp.solve(sp.Eq(sp.resultant(Ssym,sp.diff(Ssym,r),r),0),J) if s.is_real and float(s)>1][0])
def pS(mv,Jv):
    DE=np.array([E**2-1,2*mv]); p=np.polymul(np.polymul([1,0],[1,-2*mv]),DE)
    return np.polymul(p,np.polysub(np.polymul([1,0,0],[1,-2*mv]),Jv**2*DE))
Sc=pS(m,Jc); dmS=(pS(m+1e-6,Jc)-pS(m-1e-6,Jc))/2e-6
K=np.polymul([Jc],[E**2-1,2.0]); Nm=np.polysub(np.polymul(Sc,np.array([2*Jc])),0.5*np.polymul(K,dmS))
rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
rd=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2)); Q,_=np.polydiv(Sc,np.polymul([1,-rd],[1,-rd]))
er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1r,e2r,e3r,e4r=er
def Q4(x): return np.polyval(Q,x)
k2=((e3r-e2r)*(e4r-e1r))/((e4r-e2r)*(e3r-e1r)); pref=2/mp.sqrt((e4r-e2r)*(e3r-e1r))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; qn=mp.exp(mp.pi*1j*tau)
L1=lambda u: mp.jtheta(1,u,qn); L1p=lambda u: mp.jtheta(1,u,qn,1); L1pp=lambda u: mp.jtheta(1,u,qn,2); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,qn,3)/th1p0)
def wsig(z): u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*L1(u)/th1p0
def wzet(z): u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
def wp(z):   u=mp.pi*z/(2*om1); rr=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-rr**2)
def wpp(z):
    u=mp.pi*z/(2*om1); T1=L1p(u)/L1(u); T2=L1pp(u)/L1(u); T3=mp.jtheta(1,u,qn,3)/L1(u)
    return -(mp.pi/(2*om1))**3*(T3-3*T1*T2+2*T1**3)
z_inf=float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4r,mp.inf])))
z_d=z_inf+float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd])))
def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4r,rv,limit=400)[0])
z0=mp.mpf(zr(r0)); zdm=mp.mpf(z_d); zim=mp.mpf(z_inf); iw=1j*mp.mpf(w_im)
ef=[lambda z: mp.mpf(1), lambda z: wzet(z-zdm)-wzet(z+zdm), lambda z: wp(z-zdm)+wp(z+zdm),
    lambda z: wpp(z-zdm)-wpp(z+zdm), lambda z: wp(z), lambda z: wp(z-iw),
    lambda z: wp(z-zim)+wp(z+zim), lambda z: wzet(z-zim)-wzet(z+zim)]
def Pe(i,z):
    if i==0: return z-z0
    if i==1: return (mp.log(wsig(z-zdm))-mp.log(wsig(z+zdm)))-(mp.log(wsig(z0-zdm))-mp.log(wsig(z0+zdm)))
    if i==2: return (-(wzet(z-zdm)+wzet(z+zdm)))-(-(wzet(z0-zdm)+wzet(z0+zdm)))
    if i==3: return ((wp(z-zdm)-wp(z+zdm)))-((wp(z0-zdm)-wp(z0+zdm)))
    if i==4: return (-wzet(z))-(-wzet(z0))
    if i==5: return (-wzet(z-iw))-(-wzet(z0-iw))
    if i==6: return (-(wzet(z-zim)+wzet(z+zim)))-(-(wzet(z0-zim)+wzet(z0+zim)))
    if i==7: return (mp.log(wsig(z-zim))-mp.log(wsig(z+zim)))-(mp.log(wsig(z0-zim))-mp.log(wsig(z0+zim)))
Rall=[0,1,2,3,4,5]; Eall=[0,1,6,7]
pairs=[(i,j) for i in range(8) for j in range(i+1,8) if (i in Rall and j in Eall) or (i in Eall and j in Rall)]
def dA(i,j,z): return ef[i](z)*Pe(j,z)-ef[j](z)*Pe(i,z)
Wf=[]
for i in range(8):
    for j in range(i,8): Wf.append((lambda i,j:(lambda z: ef[i](z)*Pe(j,z)+ef[j](z)*Pe(i,z)))(i,j))
for i in range(8): Wf.append((lambda i:(lambda z: Pe(i,z)+z*ef[i](z)))(i))
Wf.append(lambda z: 2*z); Wf.append(lambda z: mp.mpf(1))
npts=90
log.info("valuto A' (%d) e W' (%d) su %d punti a dps=%d"%(len(pairs),len(Wf),npts,mp.mp.dps))
zs=[mp.mpf(zr(rr)) for rr in np.linspace(11.97,8.83,npts)]
Amat=mp.matrix(npts,len(pairs)); Wmat=mp.matrix(npts,len(Wf))
for a,z in enumerate(tqdm(zs,desc="A',W' hp")):
    for b,(i,j) in enumerate(pairs): Amat[a,b]=dA(i,j,z)
    for b,wf in enumerate(Wf): Wmat[a,b]=wf(z)
log.info("rango via SVD mpmath: rank([A|W]) - rank(W)")
def rankmp(M,thr=mp.mpf('1e-28')):
    U,S,V=mp.svd_r(M); S=[abs(x) for x in S]; s0=max(S)
    return sum(1 for x in S if x/s0>thr), [x/s0 for x in sorted(S,reverse=True)]
# [A|W]
AW=mp.matrix(npts,len(pairs)+len(Wf))
for a in range(npts):
    for b in range(len(pairs)): AW[a,b]=Amat[a,b]
    for b in range(len(Wf)): AW[a,len(pairs)+b]=Wmat[a,b]
rAW,svAW=rankmp(AW); rW,_=rankmp(Wmat)
print("rank([A|W]) =",rAW,"  rank(W) =",rW)
print("=> RANGO(A mod W) =",rAW-rW,"   relazioni di Fay =",len(pairs)-(rAW-rW))
print("primi valori singolari [A|W] (norm):")
for k in range(min(30,len(svAW))): print("   sv[%2d]=%s"%(k,mp.nstr(svAW[k],4)))
