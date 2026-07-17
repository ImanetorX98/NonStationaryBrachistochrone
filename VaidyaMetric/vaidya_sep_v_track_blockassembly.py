# VAIDYA ramo v: BLOCK ASSEMBLY del weight-2. delta phi_v = G~ v - int G~ v_z dz.
# v_z = E r^3/(r-r_d) + r sqrt(Q4)/(r-2m)  [2o termine DISPARI, orizzonte 2-torsione].
# Decompongo v_z in parti principali (contorno) -> int G~ v_z = sum coeff * dilog. Verifico=diretto.
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
# --- TRACKING v: N_m -> N_tot = N_m + (Jc/m) N_J,  N_J = S dJK - 1/2 K dJS = S*DE + J^2 r(r-2m) DE^3 ---
Msy=sp.symbols('M'); DEs=(Es**2-1)*r+2*Msy; Ss=r*(r-2*Msy)*DEs*(r**2*(r-2*Msy)-J**2*DEs); Ks=J*DEs
N_J_sym=sp.expand((Ss*sp.diff(Ks,J)-sp.Rational(1,2)*Ks*sp.diff(Ss,J)).subs(Msy,m).subs(J,Jc))
NJ=np.array([float(c) for c in sp.Poly(N_J_sym,r).all_coeffs()])
dJcdm=Jc/m   # = j(E), scaling lineare Vaidya
_L=max(len(Nm),len(NJ)); Nm=np.r_[np.zeros(_L-len(Nm)),Nm]+dJcdm*np.r_[np.zeros(_L-len(NJ)),NJ]  # Nm = N_tot
log.info('TRACKING v: dJc/dm=Jc/m=%.5f, uso N_tot'%dJcdm)
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
def lau(fun,aa,order,eps=1e-3):
    return {n:complex(mp.quad(lambda th: fun(aa+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi)) for n in range(1,order+1)}
PPr={'zd':(mp.mpf(z_d),lau(Rf,mp.mpf(z_d),3)),'mzd':(mp.mpf(-z_d),lau(Rf,mp.mpf(-z_d),3)),'h0':(mp.mpf(0),lau(Rf,mp.mpf(0),2)),'h2':(iw,lau(Rf,iw,2))}
C0=complex(Rf(mp.mpf('0.19'))-sum(b.get(1,0)*wzet(mp.mpf('0.19')-aa)+b.get(2,0)*wp(mp.mpf('0.19')-aa)+b.get(3,0)*(-wpp(mp.mpf('0.19')-aa)/2) for aa,b in PPr.values()))
def Gt(z):
    s=C0*z
    for nm,(aa,b) in PPr.items(): s+= b.get(1,0)*mp.log(wsig(z-aa))-b.get(2,0)*wzet(z-aa)-(b.get(3,0)/2)*wp(z-aa)
    return s
def Gt_r(x): z=mp.mpf(zr(x)); return complex(Gt(z)-Gt(z0)).real
def sqQ_z(z): return (1/sa)*(wp(z-z_inf)-wp(z+z_inf))   # = dr/dz = sqrt(Q4), meromorfa (ok z complesso)
def vz(z):    rr=r_of_z(z); return E*rr**3/(rr-rd)+rr*sqQ_z(z)/(rr-2*m)   # esplicito (dispari, sqQ=dr/dz)
log.info("decompongo v_z in parti principali (contorno)")
# v_z poli: z_d,-z_d (da r^3/(r-r_d)); z_inf,-z_inf (da r^3 e sqrt(Q4)); iw (orizzonte, termine dispari)
def vz_c(z): return complex(vz(z))
PPv={'zd':(mp.mpf(z_d),lau(vz_c,mp.mpf(z_d),1)),'mzd':(mp.mpf(-z_d),lau(vz_c,mp.mpf(-z_d),1)),
     'zi':(mp.mpf(z_inf),lau(vz_c,mp.mpf(z_inf),2)),'mzi':(mp.mpf(-z_inf),lau(vz_c,mp.mpf(-z_inf),2)),
     'iw':(iw,lau(vz_c,iw,1))}
ztv=mp.mpf('0.12')
def prcv(z): return sum(b.get(1,0)*wzet(z-aa)+b.get(2,0)*wp(z-aa) for aa,b in PPv.values())
Cv=complex(vz_c(ztv)-prcv(ztv))
# verifica ricostruzione v_z
recerr=max(abs(complex(vz_c(mp.mpf(t))-(Cv+prcv(mp.mpf(t))))) for t in [0.10,0.15,0.20])
log.info("check v_z ricostruito da parti principali: %.1e  (poli: +-z_d,+-z_inf,iw)"%recerr)
for nm,(aa,b) in PPv.items(): log.info("  v_z polo %s: res=%s  P-coeff=%s"%(nm,mp.nstr(b.get(1,0),4),mp.nstr(b.get(2,0),4)))

# --- BLOCK ASSEMBLY: int G~ v_z dz = somma atomi (G~-termini x v_z-parti-principali) ---
# G~ termini: (coeff,kind,polo); v_z: Cv + sum e1 zeta(z-a) + e2 P(z-a)
Gterms=[(C0,'z',mp.mpf(0))]
for nm,(aa,b) in PPr.items():
    Gterms+=[(b.get(1,0),'lns',aa),(-b.get(2,0),'zet',aa),(-(b.get(3,0))/2,'wp',aa)]
Vterms=[(Cv,'one',mp.mpf(0))]
for nm,(aa,b) in PPv.items():
    if abs(b.get(1,0))>1e-12: Vterms.append((b.get(1,0),'zet',aa))
    if abs(b.get(2,0))>1e-12: Vterms.append((b.get(2,0),'wp',aa))
def gf(kind,p,z):
    if kind=='z': return z
    if kind=='one': return mp.mpf(1)
    if kind=='lns': return mp.log(wsig(z-p))
    if kind=='zet': return wzet(z-p)
    if kind=='wp': return wp(z-p)
_ic={}
def inner(k,p,z):
    key=(k,float(p.real),float(p.imag),float(z))
    if key in _ic: return _ic[key]
    v=mp.quad(lambda t: gf(k,p,t),[z0,z]); _ic[key]=v; return v
def atomP(ki,pi,kj,pj,z): return mp.quad(lambda t: gf(ki,pi,t)*gf(kj,pj,t),[z0,z],maxdegree=8)  # PRODOTTO
def intGvz_atoms(z):  # int G~ v_z dz = sum_{i,j} ci cj int gf_i gf_j dz
    return sum(ci*cj*atomP(ki,pi,kj,pj,z) for ci,ki,pi in Gterms for cj,kj,pj in Vterms)
def sQ(x): return np.sqrt(Q4(x))
def U(x,k): return quad(lambda t:t**k/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def vclock(x): return E*U(x,3)+(x-r0)+2*m*(np.log(abs(x-2*m))-np.log(abs(r0-2*m)))
def src_track(x): return np.polyval(Nm,x)/np.polyval(Sc,x)**1.5   # Nm = N_tot (tracking)
def dphi_direct(x): return quad(lambda t: src_track(t)*vclock(t), r0, x, limit=200)[0]
def dphi_block(x):
    z=mp.mpf(zr(x)); return complex(Gt(z)).real*vclock(x)-complex(intGvz_atoms(z)).real
natoms=len(Gterms)*len(Vterms)
log.info("block assembly ramo v: %d prodotti (G~ %d x v_z %d)"%(natoms,len(Gterms),len(Vterms)))
print("=== VAIDYA v TRACKING BLOCK: delta phi_v = G~ v - sum atomi(G~ x v_z)  vs DIRETTO ===")
for x in tqdm([11.0,10.5,10.0],desc="v-block"):
    db=dphi_block(x); dd=dphi_direct(x)
    print(f"  r={x:4.2f}  block={db:+.8f}  dir={dd:+.8f}  diff={abs(db-dd):.1e}")
# conta atomi dilog (lns x zet)
ndl=sum(1 for ci,ki,pi in Gterms for cj,kj,pj in Vterms if ki=='lns' and kj=='zet' and abs(ci*cj)>1e-9)
print("atomi dilog (lnσ×ζ) = %d ; totale prodotti = %d ; coeff = residui G~ x v_z (chiusi)."%(ndl,natoms))
