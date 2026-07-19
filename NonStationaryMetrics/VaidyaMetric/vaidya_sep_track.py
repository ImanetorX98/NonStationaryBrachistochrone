# SEPARATRICE-TRACKING (J=Jc(m)): correzione adiabatica con Jc mobile.
# sorgente = dm F + (Jc/m) dJ F = N_tot/S^{3/2},  N_tot = N_m + (Jc/m) N_J  (polinomiale).
# STESSA macchina della tau (pole-adapted G, residui chiusi, dilog): cambia solo il numeratore.
# Verifica: delta phi_track assemblato (G esplicito) = diretto int[dmF+(Jc/m)dJF] eta dr.
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
jE=Jc/m   # dJc/dm = j(E) = Jc/m
# --- numeratori simbolici N_m, N_J, N_tot ---
Msym=sp.symbols('M',positive=True)
DEm=(Es**2-1)*r+2*Msym; Sm=r*(r-2*Msym)*DEm*(r**2*(r-2*Msym)-J**2*DEm); Km=J*DEm
N_m_sym=sp.expand((Sm*sp.diff(Km,Msym)-sp.Rational(1,2)*Km*sp.diff(Sm,Msym)).subs(Msym,m))
N_J_sym=sp.expand((Sm.subs(Msym,m))*sp.diff((Km.subs(Msym,m)),J)-sp.Rational(1,2)*(Km.subs(Msym,m))*sp.diff((Sm.subs(Msym,m)),J))
def poly_np(expr): return np.array([float(c) for c in sp.Poly(expr.subs(J,Jc),r).all_coeffs()])
Nm=poly_np(N_m_sym); NJ=poly_np(N_J_sym)
# allinea gradi e somma N_tot = N_m + jE*N_J
L=max(len(Nm),len(NJ)); Nm_p=np.r_[np.zeros(L-len(Nm)),Nm]; NJ_p=np.r_[np.zeros(L-len(NJ)),NJ]
Ntot=Nm_p+jE*NJ_p
def pS(mv,Jv):
    DE=np.array([E**2-1,2*mv]); p=np.polymul(np.polymul([1,0],[1,-2*mv]),DE)
    return np.polymul(p,np.polysub(np.polymul([1,0,0],[1,-2*mv]),Jv**2*DE))
Sc=pS(m,Jc); rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
rd=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2)); Q,_=np.polydiv(Sc,np.polymul([1,-rd],[1,-rd]))
er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1r,e2r,e3r,e4r=er
def Q4(x): return np.polyval(Q,x)
log.info("Jc=%.6f j(E)=Jc/m=%.6f  N_tot grado=%d"%(Jc,jE,len(Ntot)-1))
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
def Rt(z): rr=r_of_z(z); return np.polyval(Ntot,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))   # R_track
def lau(fun,a,order,eps=1e-3):
    return {n:complex(mp.quad(lambda th: fun(a+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi)) for n in range(1,order+1)}
log.info("parti principali di R_track (contorno)")
PPr={'zd':(mp.mpf(z_d),lau(Rt,mp.mpf(z_d),3)),'mzd':(mp.mpf(-z_d),lau(Rt,mp.mpf(-z_d),3)),
     'h0':(mp.mpf(0),lau(Rt,mp.mpf(0),2)),'h2':(iw,lau(Rt,iw,2))}
def const_of(fun,PP,zt): return complex(fun(zt)-sum(b.get(1,0)*wzet(zt-a)+b.get(2,0)*wp(zt-a)+b.get(3,0)*(-wpp(zt-a)/2) for a,b in PP.values()))
C0=const_of(Rt,PPr,mp.mpf('0.19'))
def Gt(z):
    s=C0*z
    for nm,(a,b) in PPr.items(): s+= b.get(1,0)*mp.log(wsig(z-a))-b.get(2,0)*wzet(z-a)-(b.get(3,0)/2)*wp(z-a)
    return s
# verifica dG/dz = R_track
d=max(abs(complex((Gt(mp.mpf(t)+mp.mpf('1e-7'))-Gt(mp.mpf(t)-mp.mpf('1e-7')))/mp.mpf('2e-7')-Rt(mp.mpf(t)))) for t in [0.15,0.22])
log.info("check dG_track/dz = R_track : %.1e"%d)
def Gt_r(x): z=mp.mpf(zr(x)); return complex(Gt(z)-Gt(z0)).real
# clock tau: eta=U_3-2U_2
def sQ(x): return np.sqrt(Q4(x))
def U(x,k): return quad(lambda t:t**k/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def eta(x): return U(x,3)-2*m*U(x,2)
def etap(x): return (x**3-2*m*x**2)/((x-rd)*sQ(x))   # eta' forma dz
# sorgente diretta: dm F + (Jc/m) dJ F = N_tot/S^{3/2}
Sn=lambda x: np.polyval(Sc,x)
def src(x): return np.polyval(Ntot,x)/Sn(x)**1.5
def dphi_direct(x): return quad(lambda t: src(t)*eta(t), r0, x, limit=200)[0]
def dphi_asm(x):    # IBP: G~ eta - int G~ eta'  (delta phi = int src*eta = int (dG/dr) eta)
    return Gt_r(x)*eta(x) - quad(lambda t: Gt_r(t)*etap(t), r0, x, limit=120)[0]
log.info("delta phi_track (assemblato G esplicito) vs diretto")
print("=== delta phi_track|_sep (tau clock, J=Jc(m) tracked) : ASSEMBLATO vs DIRETTO ===")
for x in tqdm([11.0,10.0,9.2],desc="track"):
    print(f"  r={x:4.1f}  asm={dphi_asm(x):+.8f}  dir={dphi_direct(x):+.8f}  diff={abs(dphi_asm(x)-dphi_direct(x)):.1e}")
print("\nStessa macchina della tau (N_m -> N_tot): residui b_n^a, C0, C_e in forma chiusa identica.")
