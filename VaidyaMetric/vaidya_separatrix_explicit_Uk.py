# SEPARATRICE Vaidya: mattoni U_k SCRITTI ESPLICITAMENTE in Weierstrass sigma,zeta,P.
# (funzioni speciali standard, non integrali in quadratura). Coefficienti da RESIDUI +
# RIDUZIONE DI HERMITE sistematica (esatta), verificati vs quadratura ~1e-9.
#
# Curva ellittica E: w^2=Q4(r), sqrt(S)=(r-r_d)sqrt(Q4), coord di Abel z=int dr/sqrt(Q4).
# Split algebrico: r^k/(r-r_d)=poly_{k-1}(r)+r_d^k/(r-r_d) => U_k = sum r_d^{k-1-i} V_i + r_d^k U_0.
#   V_i = int r^i dz = int r^i/sqrt(Q4) dr.
#   V_0=z; V_1 (3a specie inf), V_2 (2a specie) espliciti in sigma,zeta.
#   V_{>=3}: RICORSIONE di Hermite da d/dr(r^k sqrt Q4):
#     (2k+4)a4 V_{k+3}+(2k+3)b3 V_{k+2}+(2k+2)b2 V_{k+1}+(2k+1)b1 V_k+2k b0 V_{k-1} = 2 r^k sqrt(Q4)
import numpy as np, mpmath as mp
from scipy.integrate import quad
mp.mp.dps=30
r0=12.0; er=[-2.0833333333,0.0,2.0,8.7274240]; a4=0.96; e1,e2,e3,e4=er; r_d=-3.3637111
def Q4(x): return a4*(x-er[0])*(x-er[1])*(x-er[2])*(x-er[3])
a4c,b3,b2,b1,b0=[float(c) for c in np.poly(er)*a4]      # coeff Q4 (high->low)
# --- reticolo + sigma,zeta,P ESATTI da theta1 ---
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; q=mp.exp(mp.pi*1j*tau)
L1=lambda u: mp.jtheta(1,u,q); L1p=lambda u: mp.jtheta(1,u,q,1); L1pp=lambda u: mp.jtheta(1,u,q,2); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,q,3)/th1p0)
def wsig(z): z=mp.mpc(z); u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*L1(u)/th1p0
def wzet(z): z=mp.mpc(z); u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
def wp(z):   z=mp.mpc(z); u=mp.pi*z/(2*om1); r=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-r**2)
def zr(rv): return quad(lambda x:1/np.sqrt(Q4(x)), e4, rv, limit=400)[0]
sa=np.sqrt(a4); sQ=lambda x: np.sqrt(Q4(x))
a=float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4,mp.inf]))         # z_inf (immagine di infinito)
c_r=float(mp.re(e4-(2/sa)*wzet(a))); za=wzet(2*a)
Cid2=float(mp.re((wzet(mp.mpf('0.11')-a)-wzet(mp.mpf('0.11')+a))**2
      -(-2*za*(wzet(mp.mpf('0.11')-a)-wzet(mp.mpf('0.11')+a))+wp(mp.mpf('0.11')-a)+wp(mp.mpf('0.11')+a))))
z_d=a+float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,r_d]))
rho=1/np.sqrt(Q4(r_d)); C0=float(mp.re(1.0/(e4-r_d)-rho*(wzet(-z_d)-wzet(z_d))))
# --- primitive esplicite (funzioni di z) ---
LS =lambda z: mp.log(wsig(z-a)/wsig(z+a))
U0p=lambda z: rho*mp.log(wsig(z-z_d)/wsig(z+z_d))+C0*z
V1p=lambda z: c_r*z-(1/sa)*LS(z)
V2p=lambda z: c_r**2*z-(2*c_r/sa)*LS(z)+(1/a4)*(-2*za*LS(z)-wzet(z-a)-wzet(z+a)+Cid2*z)
def pd(pf,x): z=mp.mpf(zr(x)); z0=mp.mpf(zr(r0)); return float(mp.re(pf(z)-pf(z0)))
def V0(x): return zr(x)-zr(r0)
def V1(x): return pd(V1p,x)
def V2(x): return pd(V2p,x)
def V3(x): return (2*(sQ(x)-sQ(r0))       -3*b3*V2(x)-2*b2*V1(x)-1*b1*V0(x))/(4*a4c)
def V4(x): return (2*(x*sQ(x)-r0*sQ(r0))   -5*b3*V3(x)-4*b2*V2(x)-3*b1*V1(x)-2*b0*V0(x))/(6*a4c)
def V5(x): return (2*(x**2*sQ(x)-r0**2*sQ(r0))-7*b3*V4(x)-6*b2*V3(x)-5*b1*V2(x)-4*b0*V1(x))/(8*a4c)
def U0(x): return pd(U0p,x)
Vs=[V0,V1,V2,V3,V4,V5]
def Uk(x,k): return sum(r_d**(k-1-i)*Vs[i](x) for i in range(k)) + r_d**k*U0(x)
def Udir(x,k): return quad(lambda t:t**k/((t-r_d)*np.sqrt(Q4(t))), r0, x, limit=300)[0]

print("Separatrice Vaidya: tau=%.6f i, z_inf=%.5f, z_d=%.5f, c_r=%.5f, rho=%.5f"%(float(tau.imag),a,z_d,c_r,rho))
print("Mattoni U_k = funzioni speciali ESPLICITE (Weierstrass sigma,zeta,P + algebrico) vs quadratura:")
print("  r      diff(U0..U5)")
for x in [11.5,11.0,10.0,9.5,9.0]:
    dmax=max(abs(Uk(x,k)-Udir(x,k)) for k in range(6))
    print(f"  {x:4.1f}   {dmax:.1e}")
print("\nWeight-1 COMPLETO (U_0..U_5). V_{>=3} da ricorsione di Hermite esatta.")
print("Restano: L_2m (2a specie, stessa macchina) e weight-2 W_kj,D_k -> Bloch-Wigner.")
