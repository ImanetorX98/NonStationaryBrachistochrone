# SEPARATRICE Thakurta-Kerr ramo t (tempo COORDINATA). Curva R6=r Q2 DE; clock 2a+3a specie.
# delta phi_t = -Ehat(A'/A) int eta_t * dE F dr,  dE F = N_t/R6^{3/2}.
# clock t = eta2 + eta3: eta2'=P3/sqrt R6 (2a specie), eta3'=R_Delta/(Delta sqrt R6) (3a, orizzonti).
# In z (separatrice sqrt R6=(r-r_d)sqrt Q4): eta_t'(dz) = rho_t/(r-r_d), rho_t=P3+R_Delta/Delta.
import numpy as np, mpmath as mp, sympy as sp, sys, logging
from scipy.integrate import quad
try: from tqdm import tqdm
except ImportError:
    def tqdm(it,**k): return it
logging.basicConfig(level=logging.INFO,format='[%(asctime)s] %(message)s',datefmt='%H:%M:%S',stream=sys.stdout)
log=logging.getLogger(); mp.mp.dps=30
M=1.0; a=0.9; E=1.2; r0=20.0; SIGN='prograda'; Jc=19.088436  # Jc+ (r_d=-6.62, turning 11.40)
r,J,Esy,Msy,asy=sp.symbols('r J E M a'); asub={Msy:1,asy:sp.Rational(9,10)}
DE=(Esy**2-1)*r+2*Msy; Delta=r**2-2*Msy*r+asy**2
Q2=(2*Esy**2*J**2*Msy*r-Esy**2*J**2*r**2-4*Esy**2*J*Msy*asy*r+2*Esy**2*Msy*asy**2*r+Esy**2*asy**2*r**2+Esy**2*r**4
    +4*J**2*Msy**2-4*J**2*Msy*r+J**2*r**2-8*J*Msy**2*asy+4*J*Msy*asy*r+4*Msy**2*asy**2)
R6f=sp.expand(r*Q2*DE); Ktf=r*DE*(J*(r-2*Msy)+2*Msy*asy)/Delta
Ntf=sp.expand(sp.simplify(sp.diff(Ktf/sp.sqrt(R6f),Esy)*R6f**sp.Rational(3,2)))
rho_tf=sp.cancel(sp.together((Esy**2*r**3-2*Msy*asy*Ktf/r)/((r-2*Msy)/r)))
# sostituisco M,a numerici (E,J restano simbolici per ora)
R6=sp.expand(R6f.subs(asub)); Nt=sp.expand(Ntf.subs(asub)); rho_t=sp.cancel(rho_tf.subs(asub)); DeltaN=Delta.subs(asub)
# FIX: dividi per il denominatore VERO (=2500*Delta), non Delta (bug: P3 era 2500x)
P3poly,_=sp.div(sp.Poly(sp.numer(rho_t),r),sp.Poly(sp.denom(rho_t),r))
P3=sp.expand(P3poly.as_expr()); RD=sp.expand(sp.simplify((rho_t-P3)*DeltaN))   # rho_t = P3 + RD/Delta
Delta=DeltaN
log.info("simbolico OK. Jc=%.5f (%s), raffino (Jc,r_d) ad alta prec"%(Jc,SIGN))
# raffina: Q2(rd,Jc)=0 e dQ2/dr(rd,Jc)=0  (sistema doppia radice), near (r_d0,Jc)
Q2E=Q2.subs(asub).subs(Esy,sp.Rational(6,5))
Q2n=sp.lambdify((r,J),Q2E,'mpmath'); Q2rn=sp.lambdify((r,J),sp.diff(Q2E,r),'mpmath')
sol=mp.findroot(lambda rd_,Jc_:[Q2n(rd_,Jc_),Q2rn(rd_,Jc_)], (mp.mpf('-6.62'),mp.mpf(str(Jc))))
rd_hp=float(sol[0]); Jc=float(sol[1]); rd0=rd_hp
log.info("raffinato: Jc=%.10f  r_d=%.10f"%(Jc,rd_hp))
# numerico a Jc
R6c=np.array([float(c) for c in sp.Poly(R6.subs(Esy,sp.Rational(6,5)).subs(J,Jc),r).all_coeffs()])
Ntc=np.array([float(c) for c in sp.Poly(sp.expand(Nt.subs(Esy,sp.Rational(6,5))).subs(J,Jc),r).all_coeffs()])
P3n=sp.lambdify(r,P3.subs(Esy,sp.Rational(6,5)).subs(J,Jc),'numpy'); RDn=sp.lambdify(r,RD.subs(Esy,sp.Rational(6,5)).subs(J,Jc),'numpy'); Dn=sp.lambdify(r,Delta,'numpy')
rd=rd_hp; Q,rem=np.polydiv(R6c,np.polymul([1,-rd],[1,-rd]))
log.info("resto polydiv R6/(r-r_d)^2 = %.1e (deve ~0)"%np.max(np.abs(rem)))
er=np.sort(np.real(np.roots(Q))); a4=Q[0]; e1r,e2r,e3r,e4r=er
def Q4(x): return np.polyval(Q,x)
log.info("r_d=%.4f  radici Q4=%s  turning e4=%.4f  orizzonti r_pm=%.4f,%.4f"%(rd,np.round(er,3),e4r,1+np.sqrt(1-a*a),1-np.sqrt(1-a*a)))
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
def Rf(z): rr=r_of_z(z); return np.polyval(Ntc,rr)/((rr-rd)**3*np.polyval(Q,complex(rr)))
log.info("check r(z): %.1e"%max(abs(float(np.real(r_of_z(mp.mpf(zr(rv)))))-rv) for rv in [19.,18.]))
def lau(fun,aa,order,eps=1e-3):
    return {n:complex(mp.quad(lambda th: fun(aa+eps*mp.exp(1j*th))*(eps*mp.exp(1j*th))**(n-1)*(1j*eps*mp.exp(1j*th)),[0,2*mp.pi])/(2j*mp.pi)) for n in range(1,order+1)}
PPr={'zd':(mp.mpf(z_d),lau(Rf,mp.mpf(z_d),3)),'mzd':(mp.mpf(-z_d),lau(Rf,mp.mpf(-z_d),3)),'h0':(mp.mpf(0),lau(Rf,mp.mpf(0),2)),'h2':(iw,lau(Rf,iw,2))}
C0=complex(Rf(mp.mpf('0.05'))-sum(b.get(1,0)*wzet(mp.mpf('0.05')-aa)+b.get(2,0)*wp(mp.mpf('0.05')-aa)+b.get(3,0)*(-wpp(mp.mpf('0.05')-aa)/2) for aa,b in PPr.values()))
def Gt(z):
    s=C0*z
    for nm,(aa,b) in PPr.items(): s+= b.get(1,0)*mp.log(wsig(z-aa))-b.get(2,0)*wzet(z-aa)-(b.get(3,0)/2)*wp(z-aa)
    return s
log.info("check dG/dz=R_t: %.1e"%max(abs(complex((Gt(mp.mpf(t)+mp.mpf('1e-7'))-Gt(mp.mpf(t)-mp.mpf('1e-7')))/mp.mpf('2e-7')-Rf(mp.mpf(t)))) for t in [0.05,0.09]))
def Gt_r(x): z=mp.mpf(zr(x)); return complex(Gt(z)-Gt(z0)).real
def sQ(x): return np.sqrt(Q4(x))
# clock t: eta_t' = rho_t/sqrt(R6) = [P3 + RD/Delta]/((r-r_d) sqrt Q4);  eta_t = int
def etatp(x): return (P3n(x)+RDn(x)/Dn(x))/((x-rd)*sQ(x))   # dt/dr
def etat(x): return quad(lambda t: etatp(t), r0, x, limit=200)[0]
def dEF(x): return np.polyval(Ntc,x)/np.polyval(R6c,x)**1.5
def dphi_direct(x): return quad(lambda t: dEF(t)*etat(t), r0, x, limit=200)[0]
# eta_t' in z-form per IBP: rho_t/(r-r_d)  (dz = dr/sqrtQ4)
def etatp_z(x): return (P3n(x)+RDn(x)/Dn(x))/(x-rd)
def dphi_asm(x): return Gt_r(x)*etat(x)-quad(lambda t: Gt_r(t)*etatp_z(mp.mpf(zr(t)) and t)/sQ(t)*sQ(t),r0,x,limit=120)[0] if False else Gt_r(x)*etat(x)-quad(lambda t: Gt_r(t)*etatp(t),r0,x,limit=120)[0]
print("\n=== TK t-branch prograda: delta phi ~ int eta_t dEF : ASM vs DIR ===")
for x in tqdm([19.0,18.0,17.0],desc="TK-t"):
    print(f"  r={x:4.1f}  asm={dphi_asm(x):+.8e}  dir={dphi_direct(x):+.8e}  diff={abs(dphi_asm(x)-dphi_direct(x)):.1e}")
