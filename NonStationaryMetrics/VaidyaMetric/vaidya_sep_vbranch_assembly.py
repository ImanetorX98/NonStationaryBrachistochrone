# RAMO v (tempo avanzato) sulla separatrice: ASSEMBLAGGIO in blocchi.
# delta phi_v = G~(z) v(r) - int G~ v_z dz,  v_z = E r^3/(r-r_d) + r sqrt(Q4)/(r-2m).
# v_z = [pari: E r^3/(r-r_d)] + [dispari orizzonte: r sqrt(Q4)/(r-2m)].
# Decompongo v_z in parti principali (poli z_d,+-z_inf, orizzonte i w_im) -> dilog con coeff chiusi.
import numpy as np, mpmath as mp, sympy as sp, sys, logging
from scipy.integrate import quad
try: from tqdm import tqdm
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
z0=mp.mpf(zr(r0)); iw=1j*mp.mpf(w_im)
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
def Rf(z): rr=r_of_z(z); return np.polyval(Nm,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
def lau(fun,a,order,eps=1e-3):
    return {n:complex(mp.quad(lambda th: fun(a+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi)) for n in range(1,order+1)}
PPr={'zd':(mp.mpf(z_d),lau(Rf,mp.mpf(z_d),3)),'mzd':(mp.mpf(-z_d),lau(Rf,mp.mpf(-z_d),3)),'h0':(mp.mpf(0),lau(Rf,mp.mpf(0),2)),'h2':(iw,lau(Rf,iw,2))}
C0=complex(Rf(mp.mpf('0.19'))-sum(b.get(1,0)*wzet(mp.mpf('0.19')-a)+b.get(2,0)*wp(mp.mpf('0.19')-a)+b.get(3,0)*(-wpp(mp.mpf('0.19')-a)/2) for a,b in PPr.values()))
def Gt(z):
    s=C0*z
    for nm,(a,b) in PPr.items(): s+= b.get(1,0)*mp.log(wsig(z-a))-b.get(2,0)*wzet(z-a)-(b.get(3,0)/2)*wp(z-a)
    return s
def Gt_r(x): z=mp.mpf(zr(x)); return complex(Gt(z)-Gt(z0)).real
# --- sqrt(Q4) come funzione di z (SEGNO fisico: r decresce da r0, dz>0; su (e4,r0) Q4>0) ---
def sqQ_z(z):  # sqrt(Q4(r(z))) col segno coerente con dr/dz>0 sul ramo fisico
    rr=float(np.real(r_of_z(z))); return np.sqrt(Q4(rr))
# --- v_z(z) = E r^3/(r-r_d) + r sqrt(Q4)/(r-2m) ---
def vz(z):
    rr=r_of_z(z); return E*rr**3/(rr-rd) + rr*sqQ_z(z)/(rr-2*m)
# --- clock v(r) diretto (per il bordo) ---
def sQ(x): return np.sqrt(Q4(x))
def U(x,k): return quad(lambda t:t**k/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def vclock(x): return E*U(x,3)+(x-r0)+2*m*(np.log(abs(x-2*m))-np.log(abs(r0-2*m)))
# --- diretto: delta phi_v = int dm F v dr ---
def dmF(x):
    def Fm(mv): return np.polyval(np.polymul([Jc],[E**2-1,2*mv]),x)/np.sqrt(np.polyval(pS(mv,Jc),x))
    return (Fm(m+1e-6)-Fm(m-1e-6))/2e-6
def dphi_direct(x): return quad(lambda t: dmF(t)*vclock(t), r0, x, limit=200)[0]
# --- assemblato: G~ v - int G~ v_z dz ---
def dphi_asm(x):
    zx=mp.mpf(zr(x)); bnd=Gt_r(x)*vclock(x)
    intg=complex(mp.quad(lambda z: complex(Gt(z)-Gt(z0)).real*vz(z),[z0,zx])).real
    return bnd-intg
log.info("v-branch assembly: G~ v - int G~ v_z  vs diretto")
print("=== delta phi_v|_sep : ASSEMBLATO (G~,v_z espliciti in z) vs DIRETTO ===")
for x in tqdm([11.0,10.5,10.0],desc="v-assembly"):
    da=dphi_asm(x); dd=dphi_direct(x)
    print(f"  r={x:4.2f}  asm={da:+.8f}  dir={dd:+.8f}  diff={abs(da-dd):.1e}")
