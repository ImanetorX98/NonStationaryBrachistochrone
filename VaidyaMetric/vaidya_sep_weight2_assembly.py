# SEPARATRICE Vaidya - assemblaggio WEIGHT-2 (W_kj) con dilog ellittico.
# r(z) esplicito in zeta; f_k=U_k'=r^k/(r-r_d); atomi E2 = polilog ellittico Gamma-tilde(1,1)
# (= iterated_integral GiNaC, verificato 1e-18). Ogni passo verificato vs quadratura.
import numpy as np, mpmath as mp
from scipy.integrate import quad
mp.mp.dps=30
er=[-2.0833333333,0.0,2.0,8.7274240]; a4=0.96; e1,e2,e3,e4=er; r_d=-3.3637111; r0=12.0
def Q4(x): return a4*(x-er[0])*(x-er[1])*(x-er[2])*(x-er[3])
# --- reticolo + sigma,zeta,P ESATTI da theta1 (periodi fisici 2om1,2om_im) ---
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; q=mp.exp(mp.pi*1j*tau)
L1=lambda u: mp.jtheta(1,u,q); L1p=lambda u: mp.jtheta(1,u,q,1); L1pp=lambda u: mp.jtheta(1,u,q,2); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,q,3)/th1p0)
def wsig(z): z=mp.mpc(z); u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*L1(u)/th1p0
def wzet(z): z=mp.mpc(z); u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
def wp(z):   z=mp.mpc(z); u=mp.pi*z/(2*om1); r=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-r**2)
def zr(rv): return quad(lambda x:1/np.sqrt(Q4(x)), e4, rv, limit=400)[0]   # z(r) fisico
sa=np.sqrt(a4)
z_inf=float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4,mp.inf]))       # z_infty
c_r=float(mp.re(e4-(2/sa)*wzet(z_inf)))
z_d=z_inf+float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,r_d]))
rho=1/np.sqrt(Q4(r_d)); C0=float(mp.re(1.0/(e4-r_d)-rho*(wzet(-z_d)-wzet(z_d))))

# --- r(z) come funzione ellittica ESPLICITA in zeta ---
def r_of_z(z): return c_r-(1/sa)*(wzet(z-z_inf)-wzet(z+z_inf))
print("=== PASSO 1: r(z) esplicito in zeta vs r(z) numerico (inversione) ===")
print("  om1=%.6f w_im=%.6f tau=%.6fi z_inf=%.5f z_d=%.5f c_r=%.5f rho=%.5f C0=%.5f"
      %(float(om1),w_im,float(tau.imag),z_inf,z_d,c_r,rho,C0))
for rv in [11.5,11.0,10.0,9.2]:
    z=mp.mpf(zr(rv)); err=float(abs(mp.re(r_of_z(z))-rv))
    print(f"  r={rv:5.1f}  r(z)_expl={float(mp.re(r_of_z(z))):9.5f}  err={err:.2e}")

# --- f_k(z)=r^k/(r-r_d) esplicito; verifica U_k'=f_k ---
def f_k(z,k): rr=r_of_z(z); return rr**k/(rr-r_d)
print("\n=== PASSO 2: f_k=r^k/(r-r_d) esplicito vs integrando diretto ===")
for k in [0,1,2,3]:
    for rv in [11.0,9.5]:
        z=mp.mpf(zr(rv)); direct=rv**k/((rv-r_d))   # note: integrand in dz is r^k/(r-r_d)
        err=float(abs(mp.re(f_k(z,k))-direct))
        print(f"  k={k} r={rv:4.1f}  f_k={float(mp.re(f_k(z,k))):10.5f} vs {direct:10.5f}  err={err:.1e}")

# ============ PASSO 3: weight-2  W_kj = 2 int U_k dU_j - U_k U_j ============
# Primitive ESPLICITE U_k(z) in sigma,zeta (da vaidya_separatrix_explicit_Uk, verif 1e-9).
za2=wzet(2*z_inf)
Cid2=float(mp.re((wzet(mp.mpf('0.11')-z_inf)-wzet(mp.mpf('0.11')+z_inf))**2
      -(-2*za2*(wzet(mp.mpf('0.11')-z_inf)-wzet(mp.mpf('0.11')+z_inf))
        +wp(mp.mpf('0.11')-z_inf)+wp(mp.mpf('0.11')+z_inf))))
LS =lambda z: mp.log(wsig(z-z_inf)/wsig(z+z_inf))
U0p=lambda z: rho*mp.log(wsig(z-z_d)/wsig(z+z_d))+C0*z
V1p=lambda z: c_r*z-(1/sa)*LS(z)
V2p=lambda z: c_r**2*z-(2*c_r/sa)*LS(z)+(1/a4)*(-2*za2*LS(z)-wzet(z-z_inf)-wzet(z+z_inf)+Cid2*z)
def Q4f(x): return a4*(x-er[0])*(x-er[1])*(x-er[2])*(x-er[3])
sQ=lambda x: np.sqrt(Q4f(x))
a4c,b3,b2,b1,b0=[float(c) for c in np.poly(er)*a4]
z0=mp.mpf(zr(r0))
def prim0(pf,z): return mp.re(pf(z)-pf(z0))               # primitiva con base z0
def V0v(z): return z-z0
def V1v(z): return prim0(V1p,z)
def V2v(z): return prim0(V2p,z)
# V3,V4,V5 via Hermite (in r): serve sqrt(Q4) alla quota r(z). Uso r fisico dell'estremo.
def hermiteV(x):
    rr=x; s=sQ(rr); s0=sQ(r0)
    V0=zr(rr)-zr(r0); V1=float(V1v(mp.mpf(zr(rr)))); V2=float(V2v(mp.mpf(zr(rr))))
    V3=(2*(s-s0)         -3*b3*V2-2*b2*V1-1*b1*V0)/(4*a4c)
    V4=(2*(rr*s-r0*s0)    -5*b3*V3-4*b2*V2-3*b1*V1-2*b0*V0)/(6*a4c)
    V5=(2*(rr**2*s-r0**2*s0)-7*b3*V4-6*b2*V3-5*b1*V2-4*b0*V1)/(8*a4c)
    return [V0,V1,V2,V3,V4,V5]
def U0v(z): return prim0(U0p,z)
def Uk_r(x,k):   # U_k come funzione di r (estremo), via primitive esplicite
    Vs=hermiteV(x); return sum(r_d**(k-1-i)*Vs[i] for i in range(k))+r_d**k*float(U0v(mp.mpf(zr(x))))
# --- DIRETTO: U_k, W_kj per quadratura in r (sqrt S=(r-r_d)sqrt Q4) ---
def Ukd(x,k): return quad(lambda t:t**k/((t-r_d)*sQ(t)), r0, x, limit=300)[0]
def Wd(x,k,j): return quad(lambda t:(Ukd(t,k)*t**j-Ukd(t,j)*t**k)/((t-r_d)*sQ(t)), r0, x, limit=200)[0]
# --- ASSEMBLATO: 2 int U_k dU_j - U_k U_j, con U_k ESPLICITO (special functions) ---
def IntUkdUj(x,k,j):  # int_{r0}^x U_k(r) f_j(r) dz,  f_j dz = r^j/((r-r_d)sqrtQ4) dr
    return quad(lambda t: Uk_r(t,k)*t**j/((t-r_d)*sQ(t)), r0, x, limit=200)[0]
def Wasm(x,k,j): return 2*IntUkdUj(x,k,j)-Uk_r(x,k)*Uk_r(x,j)

print("\n=== PASSO 3: W_kj = 2 int U_k dU_j - U_k U_j  (esplicito vs diretto) ===")
print("  (U_k = funzioni speciali esplicite sigma,zeta,+Hermite; int U_k dU_j = polilog ell. weight-2)")
for (k,j) in [(2,3),(0,2),(0,3),(1,3),(1,2)]:
    for x in [11.0,10.0]:
        wa=Wasm(x,k,j); wd=Wd(x,k,j); print(f"  W_{k}{j}(r={x:4.1f})  asm={wa:+.6f}  dir={wd:+.6f}  diff={abs(wa-wd):.1e}")
