# TK tau separatrice: BLOCK ASSEMBLY (decompone int G~ eta' in dilog nominati).
# delta phi ~ int R_tau eta dz = 1/2 G~ eta + 1/2 sum(c_i d_j - c_j d_i) A[e_i,e_j].
# c_i = residui di R_tau (N_tau); d_j = residui di eta'; A = dilog ellittici antisimmetrici.
import numpy as np, mpmath as mp, sympy as sp, sys, logging
from scipy.integrate import quad
try: from tqdm import tqdm
except ImportError:
    def tqdm(it,**k): return it
logging.basicConfig(level=logging.INFO,format='[%(asctime)s] %(message)s',datefmt='%H:%M:%S',stream=sys.stdout)
log=logging.getLogger(); mp.mp.dps=30
M=1.0; a=0.9; E=1.2; r0=20.0
r,J=sp.symbols('r J'); Ms=1; asq=sp.Rational(81,100); Es=sp.Rational(6,5)
DE=(Es**2-1)*r+2*Ms; Delta=r**2-2*Ms*r+asq; g=sp.expand(r*Delta-J**2*DE)
Jc=float([complex(s).real for s in sp.solve(sp.Eq(sp.resultant(g,sp.diff(g,r),r),0),J) if abs(complex(s).imag)<1e-9 and complex(s).real>5][0])
log.info("Jc=%.6f, costruisco curva e sorgente"%Jc)
Sc=np.array([float(c) for c in sp.Poly(sp.expand(r*(r-2*Ms)*DE*(r*Delta-J**2*DE)).subs(J,Jc),r).all_coeffs()])
Ntau=np.array([float(c) for c in sp.Poly(sp.expand(Es*J*r**4*(r-2*Ms)**2*DE).subs(J,Jc),r).all_coeffs()])
rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
rd=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2)); Q,_=np.polydiv(Sc,np.polymul([1,-rd],[1,-rd]))
er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1r,e2r,e3r,e4r=er; Qp=np.polyder(Q); Qpp=np.polyder(Qp)
Ntaup=np.polyder(Ntau); Ntaupp=np.polyder(Ntaup)
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
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
# --- residui R_tau (analitici, stesse formule Vaidya con N_tau) ---
s=np.sqrt(Q4(rd)); a1=np.polyval(Qp,rd)/(4*s); a2=np.polyval(Qpp,rd)/12
F=lambda x: np.polyval(Ntau,x)/np.polyval(Q,x)
Fp=lambda x:(np.polyval(Ntaup,x)*np.polyval(Q,x)-np.polyval(Ntau,x)*np.polyval(Qp,x))/np.polyval(Q,x)**2
Fpp=lambda x:((np.polyval(Ntaupp,x)*np.polyval(Q,x)-np.polyval(Ntau,x)*np.polyval(Qpp,x))*np.polyval(Q,x)-2*np.polyval(Qp,x)*(np.polyval(Ntaup,x)*np.polyval(Q,x)-np.polyval(Ntau,x)*np.polyval(Qp,x)))/np.polyval(Q,x)**3
h0=F(rd); h1=Fp(rd)*s; h2=0.5*(Fpp(rd)*s**2+Fp(rd)*(np.polyval(Qp,rd)/2))
b1zd=(h2-3*a1*h1+(6*a1**2-3*a2)*h0)/s**3; b2zd=(h1-3*a1*h0)/s**3; b3zd=h0/s**3
b2h=lambda ei: np.polyval(Ntau,ei)/((ei-rd)**3)*(4/np.polyval(Qp,ei)**2)
# eta' = (r^3-2M r^2)/(r-r_d): residui
e1_zd=(rd**3-2*M*rd**2)/s
B=c_r+(1/sa)*float(mp.re(wzet(2*z_inf))); Aq=-1/sa
# eta'-residui a z_inf: eta'=(r^3-2M r^2)/(r-r_d); r->inf ~ r^2 -> stessa struttura Vaidya (P3=r^3-2M r^2)
# uso formule Vaidya con 2m->2M:
e1_zi=Aq*(2*B*0+0)  # placeholder; calcolo via contorno sotto
# blocchi pari
Zf=lambda A_: (lambda z: wzet(z-A_)-wzet(z+A_)); Pf=lambda A_: (lambda z: wp(z-A_)+wp(z+A_)); Ppf=lambda A_: (lambda z: wpp(z-A_)-wpp(z+A_))
ef=[lambda z: mp.mpf(1), Zf(zdm), Pf(zdm), Ppf(zdm), (lambda z: wp(z)), (lambda z: wp(z-iw)), Pf(zim), Zf(zim)]
def Pe(i,z):
    if i==0: return z-z0
    if i==1: return (mp.log(wsig(z-zdm))-mp.log(wsig(z+zdm)))-(mp.log(wsig(z0-zdm))-mp.log(wsig(z0+zdm)))
    if i==2: return (-(wzet(z-zdm)+wzet(z+zdm)))-(-(wzet(z0-zdm)+wzet(z0+zdm)))
    if i==3: return ((wp(z-zdm)-wp(z+zdm)))-((wp(z0-zdm)-wp(z0+zdm)))
    if i==4: return (-wzet(z))-(-wzet(z0))
    if i==5: return (-wzet(z-iw))-(-wzet(z0-iw))
    if i==6: return (-(wzet(z-zim)+wzet(z+zim)))-(-(wzet(z0-zim)+wzet(z0+zim)))
    if i==7: return (mp.log(wsig(z-zim))-mp.log(wsig(z+zim)))-(mp.log(wsig(z0-zim))-mp.log(wsig(z0+zim)))
# coeff R (c) e eta' (d) sui blocchi. Fisso C0,Ce e i residui z_inf per match numerico a un punto.
def Rex(z): rr=r_of_z(z); return np.polyval(Ntau,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
def Epex(z): rr=r_of_z(z); return (rr**3-2*M*rr**2)/(rr-rd)
# eta'-residui z_inf via principio: e2_zi, e1_zi da parti principali (contorno, 1 volta)
def lauf(fun,A_,order,eps=1e-3):
    return {n:complex(mp.quad(lambda th: fun(A_+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi)) for n in range(1,order+1)}
ez=lauf(Epex,zim,2); e2_zi=ez[2].real; e1_zi=ez[1].real
cvec=[0,b1zd,b2zd,-b3zd/2,b2h(e4r),b2h(e3r),0.0,0.0]
dvec=[0,e1_zd,0.0,0.0,0.0,0.0,e2_zi,e1_zi]
zt=mp.mpf('0.05')
cvec[0]=complex(Rex(zt)-sum(cvec[i]*ef[i](zt) for i in range(1,8)))
dvec[0]=complex(Epex(zt)-sum(dvec[i]*ef[i](zt) for i in range(1,8)))
log.info("check R,eta' ricostruiti: R=%.1e eta'=%.1e"%(
    max(abs(complex(Rex(mp.mpf(t))-sum(cvec[i]*ef[i](mp.mpf(t)) for i in range(8)))) for t in [0.03,0.07]),
    max(abs(complex(Epex(mp.mpf(t))-sum(dvec[i]*ef[i](mp.mpf(t)) for i in range(8)))) for t in [0.03,0.07])))
Rall=[0,1,2,3,4,5]; Eall=[0,1,6,7]
pairs=[(i,j) for i in range(8) for j in range(i+1,8) if (i in Rall and j in Eall) or (i in Eall and j in Rall)]
_pc={}
def Penum(i,z):
    k=(i,float(z)); 
    if k in _pc: return _pc[k]
    v=mp.quad(lambda t: ef[i](t),[z0,z]); _pc[k]=v; return v
def A_(i,j,z): return mp.quad(lambda t: ef[i](t)*Penum(j,t)-ef[j](t)*Penum(i,t),[z0,z],maxdegree=8)
def Gt(z): return sum(cvec[i]*Penum(i,z) for i in range(8))
def etaf(z): return sum(dvec[i]*Penum(i,z) for i in range(8))
pairsnz=[(i,j) for i in range(8) for j in range(i+1,8) if (cvec[i]*dvec[j]-cvec[j]*dvec[i])!=0 and abs(cvec[i]*dvec[j]-cvec[j]*dvec[i])>1e-9]
def dphi_block(z):
    tot=0.5*Gt(z)*etaf(z)
    for i,j in pairsnz: tot+=0.5*(cvec[i]*dvec[j]-cvec[j]*dvec[i])*A_(i,j,z)
    return tot
def sQ(x): return np.sqrt(Q4(x))
def U(x,k): return quad(lambda t:t**k/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def eta_dir(x): return U(x,3)-2*M*U(x,2)
def dEF(x): return np.polyval(Ntau,x)/np.polyval(Sc,x)**1.5
def dphi_direct(x): return quad(lambda t: dEF(t)*eta_dir(t), r0, x, limit=200)[0]
log.info("block assembly (%d atomi A) vs diretto"%len(pairsnz))
print("=== TK tau BLOCK ASSEMBLY: 1/2 G~eta + 1/2 sum(cd-cd)A[e_i,e_j] vs DIRETTO ===")
for rr in tqdm([19.0,18.0,17.0],desc="TK-block"):
    z=mp.mpf(zr(rr)); db=complex(dphi_block(z)).real; dd=dphi_direct(rr)
    print(f"  r={rr:4.1f}  block={db:+.8e}  dir={dd:+.8e}  diff={abs(db-dd):.1e}")
print("atomi A[e_i,e_j] (dilog ellittici) = %d ; coeff = prodotti di residui chiusi (N_tau)."%len(pairsnz))
