# (1) ESTRAZIONE forma esplicita a 5 termini di delta phi|_sep via riduzione di Fay.
# delta phi = 1/2 G~ eta + sum_{b=1}^5 beta_b A_base_b + [peso-1].
# beta da ALGEBRA LINEARE ESATTA (residuo ~1e-15, non fit): risolvo su derivate (forma chiusa).
import numpy as np, mpmath as mp, sympy as sp, sys, logging
from scipy.integrate import quad
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(it,**k): return it
logging.basicConfig(level=logging.INFO,format='[%(asctime)s] %(message)s',datefmt='%H:%M:%S',stream=sys.stdout)
log=logging.getLogger(); mp.mp.dps=30
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
# residui (coeff R e eta')
s=np.sqrt(Q4(rd)); a1=np.polyval(Qp,rd)/(4*s); a2=np.polyval(Qpp,rd)/12
F=lambda x: np.polyval(Nm,x)/np.polyval(Q,x)
Fp=lambda x:(np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x))/np.polyval(Q,x)**2
Fpp=lambda x:((np.polyval(Nmpp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qpp,x))*np.polyval(Q,x)
   -2*np.polyval(Qp,x)*(np.polyval(Nmp,x)*np.polyval(Q,x)-np.polyval(Nm,x)*np.polyval(Qp,x)))/np.polyval(Q,x)**3
h0=F(rd); h1=Fp(rd)*s; h2=0.5*(Fpp(rd)*s**2+Fp(rd)*(np.polyval(Qp,rd)/2))
b1zd=(h2-3*a1*h1+(6*a1**2-3*a2)*h0)/s**3; b2zd=(h1-3*a1*h0)/s**3; b3zd=h0/s**3
def b2h(ei): return np.polyval(Nm,ei)/((ei-rd)**3)*(4/np.polyval(Qp,ei)**2)
e1_zd=(rd**3-2*rd**2)/s; B=c_r+(1/sa)*float(mp.re(wzet(2*z_inf))); Aa=-1/sa
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
cvec=[0,b1zd,b2zd,-b3zd/2,b2h(e4r),b2h(e3r),0.0,0.0]
dvec=[0,e1_zd,0.0,0.0,0.0,0.0,e2_zi,e1_zi]
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
def Rex(z): rr=r_of_z(z); return np.polyval(Nm,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
def Eex(z): rr=r_of_z(z); return (rr**3-2*rr**2)/(rr-rd)
zt=mp.mpf('0.19'); cvec[0]=complex(Rex(zt)-sum(cvec[i]*ef[i](zt) for i in range(1,8)))
dvec[0]=complex(Eex(zt)-sum(dvec[i]*ef[i](zt) for i in range(1,8)))
Rall=[0,1,2,3,4,5]; Eall=[0,1,6,7]
pairs=[(i,j) for i in range(8) for j in range(i+1,8) if (i in Rall and j in Eall) or (i in Eall and j in Rall)]
enam=['1','Z_zd','P_zd','Pp_zd','wp0','wpiw','P_zi','Z_zi']
wk=[0.5*(cvec[i]*dvec[j]-cvec[j]*dvec[i]) for (i,j) in pairs]     # coeff di A_k in delta phi_w2
def dApair(i,j,z): return ef[i](z)*Pe(j,z)-ef[j](z)*Pe(i,z)      # A'_{ij}
# base peso-1 V e derivate V'
Vp=[]
for i in range(8):
    for j in range(i,8): Vp.append((lambda i,j:(lambda z: ef[i](z)*Pe(j,z)+ef[j](z)*Pe(i,z)))(i,j))
for i in range(8): Vp.append((lambda i:(lambda z: Pe(i,z)+z*ef[i](z)))(i))
Vp.append(lambda z: 2*z); Vp.append(lambda z: mp.mpf(1))
# QR-pivot per scegliere 5 pair-base
log.info("setup curva: Jc=%.6f r_d=%.6f e4=%.4f  |  costruisco matrici derivate (forma chiusa)"%(Jc,rd,e4r))
zsQ=[mp.mpf(zr(rr)) for rr in np.linspace(11.9,8.9,60)]
Amat=np.array([[complex(dApair(i,j,z)) for (i,j) in pairs] for z in tqdm(zsQ,desc="QR-pivot sampling")])
Vmat=np.array([[complex(v(z)) for v in Vp] for z in zsQ])
Aperp=Amat-Vmat@np.linalg.lstsq(Vmat,Amat,rcond=None)[0]
import scipy.linalg as sla
_,_,P=sla.qr(Aperp,pivoting=True)
basis=list(P[:5]); print("5 dilog-base scelti:",[f"A[{enam[pairs[b][0]]},{enam[pairs[b][1]]}]" for b in basis])
# risolvi delta phi_w2' = sum_b beta_b A'_base + sum_l gamma_l V'_l   (derivate, forma chiusa)
log.info("risolvo sistema lineare esatto su 160 punti (derivate)")
zs=[mp.mpf(zr(rr)) for rr in np.linspace(11.95,8.85,160)]
LHS=np.array([complex(sum(wk[k]*dApair(pairs[k][0],pairs[k][1],z) for k in range(len(pairs)))) for z in tqdm(zs,desc="LHS delta_phi_w2'")])
Des=np.array([[complex(dApair(pairs[b][0],pairs[b][1],z)) for b in basis]+[complex(v(z)) for v in Vp] for z in tqdm(zs,desc="design matrix")])
coef,res_,rk_,sv_=np.linalg.lstsq(Des,LHS,rcond=None)
resid=np.max(np.abs(Des@coef-LHS))
beta=coef[:5]
print(f"residuo sistema (deve essere ~1e-13): {resid:.1e}")
print("beta (coeff dei 5 dilog-base):",np.array2string(beta.real,precision=6))
np.save('/tmp/fay_beta.npy',{'basis':basis,'beta':beta,'gamma':coef[5:],'pairs':pairs,'enam':enam})

# ===== VERIFICA END-TO-END: delta phi = 1/2 G~ eta + sum_5 beta A_base + peso-1  vs diretto =====
gamma=coef[5:]
# primitive integrate del peso-1 (V_int): {Pe_i Pe_j}, {z Pe_i}, {z^2}, {z}
def Vint(l,z):
    idx=0
    for i in range(8):
        for j in range(i,8):
            if idx==l: return Pe(i,z)*Pe(j,z); 
            idx+=1
    for i in range(8):
        if idx==l: return z*Pe(i,z)
        idx+=1
    if idx==l: return z*z
    idx+=1
    if idx==l: return z
def Abase_true(b,z):
    i,j=pairs[basis[b]]
    return mp.quad(lambda t: ef[i](t)*Pe(j,t)-ef[j](t)*Pe(i,t),[z0,z],maxdegree=7)
def Gt(z): return sum(cvec[i]*Pe(i,z) for i in range(8))
def etaf(z): return sum(dvec[i]*Pe(i,z) for i in range(8))
def dphi_5term(z):
    w2=sum(beta[b]*complex(Abase_true(b,z)) for b in range(5))
    w2+=sum(gamma[l]*complex(Vint(l,z)-Vint(l,z0)) for l in range(len(gamma)))
    return 0.5*complex(Gt(z))*complex(etaf(z))+w2
def sQ(x): return np.sqrt(Q4(x))
def U(x,kk): return quad(lambda t:t**kk/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def eta_dir(x): return U(x,3)-2*U(x,2)
def dmF(x):
    def Fm(mv): return np.polyval(np.polymul([Jc],[E**2-1,2*mv]),x)/np.sqrt(np.polyval(pS(mv,Jc),x))
    return (Fm(m+1e-6)-Fm(m-1e-6))/2e-6
def dphi_direct(x): return quad(lambda t: dmF(t)*eta_dir(t), r0, x, limit=200)[0]
log.info("verifica end-to-end su integrali annidati veri (5 dilog-base)")
print("\n=== delta phi = 1/2 G~eta + 5 dilog-base + peso-1  vs DIRETTO ===")
for rr in tqdm([11.0,10.0,9.2],desc="verifica nested"):
    z=mp.mpf(zr(rr)); d5=dphi_5term(z).real; dd=dphi_direct(rr)
    print(f"  r={rr:4.1f}  5-term={d5:+.8f}  diretto={dd:+.8f}  diff={abs(d5-dd):.1e}")
