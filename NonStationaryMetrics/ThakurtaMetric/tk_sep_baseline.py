# SEPARATRICE Thakurta-Kerr (tau, tempo proprio). Stesso pipeline Vaidya, sorgente E-deriv.
# delta phi_TK = -Ehat (A'/A) int eta * dE F dr,  dE F = N_tau/S^{3/2}, N_tau=E J r^4(r-2M)^2 DE.
# clock eta=U_3-2M U_2. Curva S=r(r-2M)DE(rDelta-J^2 DE), separatrice Jc (radice doppia cubico).
import numpy as np, mpmath as mp, sympy as sp, sys, logging
from scipy.integrate import quad
try: from tqdm import tqdm
except ImportError:
    def tqdm(it,**k): return it
logging.basicConfig(level=logging.INFO,format='[%(asctime)s] %(message)s',datefmt='%H:%M:%S',stream=sys.stdout)
log=logging.getLogger(); mp.mp.dps=30
M=1.0; a=0.9; E=1.2; r0=20.0
r,J=sp.symbols('r J'); Ms=1; asq=sp.Rational(81,100); Es=sp.Rational(6,5)
DE=(Es**2-1)*r+2*Ms; Delta=r**2-2*Ms*r+asq
g=sp.expand(r*Delta-J**2*DE)   # cubico
Jc=float([complex(s).real for s in sp.solve(sp.Eq(sp.resultant(g,sp.diff(g,r),r),0),J) if abs(complex(s).imag)<1e-9 and complex(s).real>5][0])
log.info("Jc(TK)=%.6f"%Jc)
S_sym=sp.expand(r*(r-2*Ms)*DE*(r*Delta-J**2*DE)).subs(J,Jc)
Sc=np.array([float(c) for c in sp.Poly(S_sym,r).all_coeffs()])
# N_tau = E J r^4 (r-2M)^2 DE  (sorgente dE F)
Ntau_sym=sp.expand(Es*J*r**4*(r-2*Ms)**2*DE).subs(J,Jc)
Ntau=np.array([float(c) for c in sp.Poly(Ntau_sym,r).all_coeffs()])
rts=np.roots(Sc); pr=[(i,j) for i in range(6) for j in range(i+1,6) if abs(rts[i]-rts[j])<1e-6]
rd=float(np.real((rts[pr[0][0]]+rts[pr[0][1]])/2)); Q,_=np.polydiv(Sc,np.polymul([1,-rd],[1,-rd]))
er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1r,e2r,e3r,e4r=er
def Q4(x): return np.polyval(Q,x)
log.info("r_d=%.5f  radici Q4=%s  a4=%.4f  turning e4=%.4f"%(rd,np.round(er,4),a4,e4r))
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
def Rf(z): rr=r_of_z(z); return np.polyval(Ntau,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
# check r(z)
log.info("check r(z): err=%.1e"%max(abs(float(np.real(r_of_z(mp.mpf(zr(rv)))))-rv) for rv in [19.,18.,17.]))
def lau(fun,aa,order,eps=1e-3):
    return {n:complex(mp.quad(lambda th: fun(aa+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi)) for n in range(1,order+1)}
PPr={'zd':(mp.mpf(z_d),lau(Rf,mp.mpf(z_d),3)),'mzd':(mp.mpf(-z_d),lau(Rf,mp.mpf(-z_d),3)),'h0':(mp.mpf(0),lau(Rf,mp.mpf(0),2)),'h2':(iw,lau(Rf,iw,2))}
C0=complex(Rf(mp.mpf('0.05'))-sum(b.get(1,0)*wzet(mp.mpf('0.05')-aa)+b.get(2,0)*wp(mp.mpf('0.05')-aa)+b.get(3,0)*(-wpp(mp.mpf('0.05')-aa)/2) for aa,b in PPr.values()))
def Gt(z):
    s=C0*z
    for nm,(aa,b) in PPr.items(): s+= b.get(1,0)*mp.log(wsig(z-aa))-b.get(2,0)*wzet(z-aa)-(b.get(3,0)/2)*wp(z-aa)
    return s
dchk=max(abs(complex((Gt(mp.mpf(t)+mp.mpf('1e-7'))-Gt(mp.mpf(t)-mp.mpf('1e-7')))/mp.mpf('2e-7')-Rf(mp.mpf(t)))) for t in [0.05,0.09])
log.info("check dG/dz=R_tau: %.1e"%dchk)
def Gt_r(x): z=mp.mpf(zr(x)); return complex(Gt(z)-Gt(z0)).real
def sQ(x): return np.sqrt(Q4(x))
def U(x,k): return quad(lambda t:t**k/((t-rd)*sQ(t)),r0,x,limit=200)[0]
def eta(x): return U(x,3)-2*M*U(x,2)
def etap(x): return (x**3-2*M*x**2)/((x-rd)*sQ(x))
def Sn(x): return np.polyval(Sc,x)
def dEF(x): return np.polyval(Ntau,x)/Sn(x)**1.5
def dphi_direct(x): return quad(lambda t: dEF(t)*eta(t), r0, x, limit=200)[0]
def dphi_asm(x): return Gt_r(x)*eta(x)-quad(lambda t: Gt_r(t)*etap(t), r0, x, limit=120)[0]
print("\n=== TK tau separatrice: I = int eta dEF dr  (delta phi = -Ehat A'/A * I) : ASM vs DIR ===")
for x in tqdm([19.0,18.0,17.0],desc="TK-tau"):
    print(f"  r={x:4.1f}  asm={dphi_asm(x):+.8e}  dir={dphi_direct(x):+.8e}  diff={abs(dphi_asm(x)-dphi_direct(x)):.1e}")
