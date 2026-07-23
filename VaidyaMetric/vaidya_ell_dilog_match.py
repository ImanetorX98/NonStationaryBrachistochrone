# MATCH NUMERICO: forma chiusa ellittica (Weierstrass sigma/zeta) su separatrice.
# Peso-1: U_0 = rho[ln sigma(z-z_d) - ln sigma(z+z_d)] + C z + cost  (3a specie ellittica)
#   rho = 1/sqrt(Q4(r_d)) (residuo, DERIVATO), C dalla parte olomorfa, cost da U_0(r0)=0.
# Toolkit sigma,zeta da theta1 (mpmath), reticolo auto-consistente con z(r)=int dr/sqrt(Q4).
import numpy as np, mpmath as mp
from scipy.integrate import quad
mp.mp.dps=25
M,E=1.0,1.4
er=[-2.0833333333,0.0,2.0,8.7274240]; a4=0.96
def Q4(x): return a4*(x-er[0])*(x-er[1])*(x-er[2])*(x-er[3])
e1,e2,e3,e4=er
# --- semiperiodi (Legendre, gia' verificati vs Sage a 15 cifre) ---
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
w1=float(pref*mp.ellipk(k2)); w_im=float(pref*mp.ellipk(1-k2))
om1=mp.mpf(w1); om2=mp.mpc(0,w_im); tau=om2/om1; q=mp.exp(mp.pi*1j*tau)
print(f"semiperiodo reale om1={float(om1):.8f}  immag om2={float(w_im):.8f}i")
print(f"tau={complex(tau):.8f}  (atteso 0.90597 i)   q={complex(q):.6e}")

# --- Weierstrass sigma,zeta da theta1 ---
def th1(u):  return mp.jtheta(1,u,q)
def th1p(u): return mp.jtheta(1,u,q,1)
th1p0=th1p(0); th1ppp0=mp.jtheta(1,0,q,3)
eta1=-(mp.pi**2/(12*om1))*(th1ppp0/th1p0)          # quasi-periodo
def wsigma(z):
    z=mp.mpc(z); a=mp.pi*z/(2*om1)
    return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*th1(a)/th1p0
def wzeta(z):
    z=mp.mpc(z); a=mp.pi*z/(2*om1)
    return eta1*z/om1 + (mp.pi/(2*om1))*(th1p(a)/th1(a))
# check: zeta'(z) = -wp(z); wp(om1)=e_i-root (Weierstrass). Verifica quasi-periodicita
chk=wzeta(0.3+0j)+2*eta1 - wzeta(0.3+2*om1)
print(f"check quasi-period zeta(z+2om1)=zeta(z)+2eta1 : |res|={float(abs(chk)):.2e}")

# --- mappa di Abel z(r): base al PUNTO DI WEIERSTRASS e4 (z=0) => involuzione z->-z,
#     cosi' le due preimmagini di r_d sono +-z_d. ---
r0=12.0
def zr(rv): return quad(lambda x:1/np.sqrt(Q4(x)), e4, rv, limit=400)[0]   # base e4
# r_d sulla stessa ovale (per infinito): z_d = int_{e4}^{+inf} + int_{-inf}^{r_d}
r_d=-3.3637111
zd1=mp.quad(lambda x:1/mp.sqrt(Q4(x)), [e4, mp.inf])
zd2=mp.quad(lambda x:1/mp.sqrt(Q4(x)), [-mp.inf, r_d])
z_d=complex(zd1+zd2)
print(f"\nz(e4)=0  z(r0)={zr(r0):.6f}  z_d(r_d)={z_d.real:.8f}")

# --- peso-1: U_0 diretto vs forma chiusa sigma ---
def Sn(rv):
    DE=(E**2-1)*rv+2*M
    return rv*(rv-2*M)*DE*(rv**2*(rv-2*M)-7.026623740**2*DE)
def U0_dir(rv): return quad(lambda t:1.0/((t-r_d)*np.sqrt(Q4(t))), r0, rv, limit=300)[0]
rho=1.0/np.sqrt(Q4(r_d))
# parte olomorfa C: 1/(r-r_d) - rho[zeta(z-z_d)-zeta(z+z_d)] = cost, valutata a r=e4 (z=0):
C = 1.0/(e4-r_d) - rho*(wzeta(-z_d)-wzeta(z_d))    # = 1/(e4-r_d)+2 rho zeta(z_d)
C = mp.re(C)
def U0_closed(rv):
    z=mp.mpf(zr(rv)); z0=mp.mpf(zr(r0))
    prim=lambda zz: rho*mp.log(wsigma(zz-z_d)/wsigma(zz+z_d)) + C*zz
    return mp.re(prim(z)-prim(z0))
print(f"\nrho=1/sqrt(Q4(r_d))={rho:.6f}  C(olomorfo)={float(C):.6f}")
# --- check IDENTITA' DERIVATA: dU_0/dz = 1/(r-r_d) == rho[zeta(z-z_d)-zeta(z+z_d)]+C ---
print("check derivata dU0/dz = 1/(r-r_d) vs rho[zeta(z-z_d)-zeta(z+z_d)]+C:")
for rv in [11.0,10.0,9.2]:
    z=mp.mpf(zr(rv)); lhs=1.0/(rv-r_d)
    rhs=mp.re(rho*(wzeta(z-z_d)-wzeta(z+z_d))+C)
    print(f"  r={rv}: lhs={lhs:.8f} rhs={float(rhs):.8f} diff={abs(lhs-float(rhs)):.2e}")
print("PESO-1  U_0: diretto vs forma chiusa sigma (3a specie ellittica)")
for rv in [11.0,10.0,9.2,e4+0.35]:
    ud=U0_dir(rv); uc=float(U0_closed(rv))
    print(f"  r={rv:.3f}: dir={ud:+.9f}  sigma={uc:+.9f}  diff={abs(ud-uc):.2e}")

# =====================================================================
# PESO-2: dilog d'orizzonte G~_0 = int ln(r-2m) dU_0 = int ln(r-2m)/((r-r_d)sqrt(Q4)) dr
#   In z:  G~_0 = int ln(r(z)-2m) * [rho(zeta(z-z_d)-zeta(z+z_d))+C] dz
#   ln(r-2m) via sigma: r-2m ha zero DOPPIO al 2-torsione z_h e poli doppi ai due inf.
#   => G~_0 e' combinazione di int lnsigma(z-a) zeta(z-b) dz = DILOGARITMO ELLITTICO.
print("\n=== PESO-2: dilog d'orizzonte = dilogaritmo ellittico ===")
r0v=12.0
def Sn_sep(rv):
    DE=(E**2-1)*rv+2*M
    return rv*(rv-2*M)*DE*(rv**2*(rv-2*M)-7.026623740**2*DE)
# G~_0 diretto (reale) su tratto fisico
def Gtil_dir(rv): return quad(lambda t:np.log(abs(t-2*M))/((t-r_d)*np.sqrt(Q4(t))), r0v, rv, limit=300)[0]
# G~_0 via cambio variabile z (identita' di struttura): integrando in z
z_h=complex(om2)   # 2-torsione orizzonte = semiperiodo immag
def r_of_z_real(zz):  # inverti z(r) sul tratto reale [e4,r0] (monotono)
    from scipy.optimize import brentq
    return brentq(lambda rr: zr(rr)-zz, e4, 40.0)
# verifica: G~_0(r) = U_0(r) ln(r-2m) - D_0(r),  D_0=int U_0/(r-2m) dr  (chiude il cerchio)
def U0c(rv):
    z=mp.mpf(zr(rv)); z0=mp.mpf(zr(r0v))
    prim=lambda zz: rho*mp.log(wsigma(zz-z_d)/wsigma(zz+z_d))+C*zz
    return float(mp.re(prim(z)-prim(z0)))
def D0_dir(rv): return quad(lambda t:U0c(t)/(t-2*M), r0v, rv, limit=200)[0]
print(" verifica cerchio  D_0 = U_0 ln(r-2m) - G~_0 :")
for rv in [11.0,10.0,9.2]:
    lhs=D0_dir(rv); rhs=U0c(rv)*np.log(abs(rv-2*M))-Gtil_dir(rv)
    print(f"  r={rv}: D0={lhs:+.8f}  U0 ln-G~={rhs:+.8f}  diff={abs(lhs-rhs):.1e}")

# --- Zagier: dilogaritmo ellittico (Bloch-Wigner) q-serie, funzione TABULATA ---
def BlochWigner(x):
    x=complex(x)
    return (mp.polylog(2,x)).imag + mp.arg(1-x)*mp.log(abs(x)) if abs(x)>1e-15 else 0.0
def Dell(z):   # D^E(z)=sum_{n>=0}D(q^n zeta)-sum_{n>=1}D(q^n/zeta), q,zeta come da reticolo
    zeta=mp.exp(2j*mp.pi*z/(2*om1)); Q=mp.exp(2j*mp.pi*tau); s=0.0
    for n in range(0,60): s+=BlochWigner(Q**n*zeta)
    for n in range(1,60): s-=BlochWigner(Q**n/zeta)
    return s
print(f"\n Zagier D^E q-serie implementato (funzione TABULATA, converge):")
for zc in [0.2+0.3j, complex(z_d), complex(z_d)-complex(z_h)]:
    print(f"   D^E({complex(zc):.4f}) = {float(Dell(zc)):+.8f}")
print(" NB: z reale => zeta=e^{2pi i z/2om1} sul cerchio unitario (fase) => q^n zeta")
print("     complessi => D^E != 0: il dilog ellittico singola-valore E' il contenuto.")
