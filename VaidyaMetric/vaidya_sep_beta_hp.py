# beta ad ALTA PRECISIONE (mpmath dps=40): se delta phi_w2 sta nel span dei 5 dilog-base,
# il residuo del sistema -> ~1e-38 e beta ha 38 cifre -> PSLQ per forma chiusa.
import numpy as np, mpmath as mp, sympy as sp, sys, logging
from scipy.integrate import quad
try: from tqdm import tqdm
except ImportError:
    def tqdm(it,**k): return it
logging.basicConfig(level=logging.INFO,format='[%(asctime)s] %(message)s',datefmt='%H:%M:%S',stream=sys.stdout)
log=logging.getLogger(); mp.mp.dps=45
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
Qp=np.polyder(Q); Qpp=np.polyder(Qp); Nmp=np.polyder(Nm); Nmpp=np.polyder(Nmp)
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
z_inf=float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4r,mp.inf]))); sa=float(mp.sqrt(a4))
c_r=float(mp.re(e4r-(2/sa)*wzet(z_inf))); z_d=z_inf+float(mp.re(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,rd])))
def zr(rv): return float(quad(lambda x:1/np.sqrt(Q4(x)),e4r,rv,limit=400)[0])
z0=mp.mpf(zr(r0)); zdm=mp.mpf(z_d); zim=mp.mpf(z_inf); iw=1j*mp.mpf(w_im)
s=np.sqrt(Q4(rd)); a1=np.polyval(Qp,rd)/(4*s); a2=np.polyval(Qpp,rd)/12
F=lambda x: np.polyval(Nm,x)/np.polyval(Q,x)
Fp=lambda x:(np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x))/np.polyval(Q,x)**2
Fpp=lambda x:((np.polyval(Nmpp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qpp,x))*np.polyval(Q,x)-2*np.polyval(Qp,x)*(np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x)))/np.polyval(Q,x)**3
h0=F(rd); h1=Fp(rd)*s; h2=0.5*(Fpp(rd)*s**2+Fp(rd)*(np.polyval(Qp,rd)/2))
b1zd=mp.mpf(h2-3*a1*h1+(6*a1**2-3*a2)*h0)/mp.mpf(s)**3; b2zd=mp.mpf(h1-3*a1*h0)/mp.mpf(s)**3; b3zd=mp.mpf(h0)/mp.mpf(s)**3
def b2h(ei): return mp.mpf(np.polyval(Nm,ei))/mp.mpf((ei-rd)**3)*(4/mp.mpf(np.polyval(Qp,ei))**2)
e1_zd=mp.mpf(rd**3-2*rd**2)/mp.mpf(s); B=c_r+(1/sa)*mp.re(wzet(2*z_inf)); Aa=mp.mpf(-1)/sa
e1_zi=Aa*(2*B+rd-2); e2_zi=Aa**2
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
cvec=[mp.mpf(0),b1zd,b2zd,-b3zd/2,b2h(e4r),b2h(e3r),mp.mpf(0),mp.mpf(0)]
dvec=[mp.mpf(0),e1_zd,mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(0),e2_zi,e1_zi]
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
def Rex(z): rr=r_of_z(z); return np.polyval(Nm,complex(rr))/((complex(rr)-rd)**3*np.polyval(Q,complex(rr)))
def Eex(z): rr=r_of_z(z); return (rr**3-2*rr**2)/(rr-rd)
zt=mp.mpf('0.19'); cvec[0]=Rex(zt)-sum(cvec[i]*ef[i](zt) for i in range(1,8)); dvec[0]=Eex(zt)-sum(dvec[i]*ef[i](zt) for i in range(1,8))
Rall=[0,1,2,3,4,5]; Eall=[0,1,6,7]
pairs=[(i,j) for i in range(8) for j in range(i+1,8) if (i in Rall and j in Eall) or (i in Eall and j in Rall)]
enam=['1','Z_zd','P_zd','Pp_zd','wp0','wpiw','P_zi','Z_zi']
wk=[(cvec[i]*dvec[j]-cvec[j]*dvec[i])/2 for (i,j) in pairs]
def dA(i,j,z): return ef[i](z)*Pe(j,z)-ef[j](z)*Pe(i,z)
basis_names=[('wp0','P_zi'),('wp0','Z_zi'),('Z_zd','wp0'),('1','wp0'),('Pp_zd','P_zi')]
basis=[pairs.index((min(enam.index(a),enam.index(b)),max(enam.index(a),enam.index(b)))) for a,b in basis_names]
Wf=[]
for i in range(8):
    for j in range(i,8): Wf.append((lambda i,j:(lambda z: ef[i](z)*Pe(j,z)+ef[j](z)*Pe(i,z)))(i,j))
for i in range(8): Wf.append((lambda i:(lambda z: Pe(i,z)+z*ef[i](z)))(i))
Wf.append(lambda z: 2*z); Wf.append(lambda z: mp.mpf(1))
ncol=5+len(Wf); npts=ncol+15
log.info("mpmath dps=%d: sistema %d punti x %d colonne (5 dilog-base + %d peso-1)"%(mp.mp.dps,npts,ncol,len(Wf)))
zs=[mp.mpf(zr(rr)) for rr in np.linspace(11.97,8.83,npts)]
Amat=mp.matrix(npts,ncol); bvec=mp.matrix(npts,1)
for a,z in enumerate(tqdm(zs,desc="build hp system")):
    for c,bidx in enumerate(basis): i,j=pairs[bidx]; Amat[a,c]=dA(i,j,z)
    for c,wf in enumerate(Wf): Amat[a,5+c]=wf(z)
    bvec[a]=sum(wk[k]*dA(pairs[k][0],pairs[k][1],z) for k in range(len(pairs)))
log.info("design = [5 A-base] + [colonne W indipendenti]")
import scipy.linalg as sla
Wnp=np.array([[complex(Amat[a,5+c]) for c in range(len(Wf))] for a in range(npts)])
_,_,Pw=sla.qr(Wnp,pivoting=True)
Rw=np.abs(np.diag(sla.qr(Wnp[:,Pw],mode='r')[0])); rkW=int(np.sum(Rw/Rw[0]>1e-12))
wkeep=list(Pw[:rkW])
cols=list(range(5))+[5+c for c in wkeep]     # 5 A-base + W indip
log.info("W indipendenti = %d; design totale %d colonne"%(rkW,len(cols)))
Ar=mp.matrix(npts,len(cols))
for a in range(npts):
    for c,cc in enumerate(cols): Ar[a,c]=Amat[a,cc]
sol_r=mp.lu_solve(Ar.T*Ar, Ar.T*bvec)
resid=max(abs(x) for x in (Ar*sol_r-bvec))
sol=mp.matrix(ncol,1)
for c,cc in enumerate(cols): sol[cc]=sol_r[c]
print("residuo sistema (mpmath, base indip.) =",mp.nstr(resid,4))
print(">>> se ~1e-40: delta_phi_w2 sta ESATTAMENTE nei 5 dilog -> beta esatti; se ~1e-9: rango>5")
print("beta (5 dilog-base) ad alta precisione:")
for b in range(5): print("   beta[%d] (%s,%s) = %s"%(b,basis_names[b][0],basis_names[b][1],mp.nstr(sol[b],35)))
