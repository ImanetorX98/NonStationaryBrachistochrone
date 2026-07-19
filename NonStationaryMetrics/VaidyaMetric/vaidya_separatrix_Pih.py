# Building block d'orizzonte 2a specie sulla separatrice, ESPLICITO in Weierstrass zeta.
# Pi_h = int dr/((r-2m) sqrt(Q4)),  r=2m branch point di Q4 (2-torsione z_h).
#   Pi_h = beta*zeta(z-z_h) + gamma*z,  beta=-4/Q4'(2m) (dal polo doppio, analitico),
#   z_h = mezzo-periodo immaginario, gamma costante da identita' derivata. Verif 1e-14.
import numpy as np, mpmath as mp
from scipy.integrate import quad
mp.mp.dps=30
r0=12.0; er=[-2.0833333333,0.0,2.0,8.7274240]; a4=0.96; e1,e2,e3,e4=er
def Q4(x): return a4*(x-er[0])*(x-er[1])*(x-er[2])*(x-er[3])
def Q4p(x):
    s=0.0
    for i in range(4):
        p=a4
        for j in range(4):
            if j!=i: p*=(x-er[j])
        s+=p
    return s
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; q=mp.exp(mp.pi*1j*tau)
L1=lambda u: mp.jtheta(1,u,q); L1p=lambda u: mp.jtheta(1,u,q,1); L1pp=lambda u: mp.jtheta(1,u,q,2); th1p0=L1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,q,3)/th1p0)
def wzet(z): z=mp.mpc(z); u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(L1p(u)/L1(u))
def wp(z): z=mp.mpc(z); u=mp.pi*z/(2*om1); r=L1p(u)/L1(u); return -eta1/om1-(mp.pi/(2*om1))**2*(L1pp(u)/L1(u)-r**2)
def zr(rv): return quad(lambda x:1/np.sqrt(Q4(x)), e4, rv, limit=400)[0]
z_h=mp.mpc(0, quad(lambda x:1/np.sqrt(-Q4(x)),e3,e4,limit=300)[0])   # 2-torsione orizzonte
beta=-4.0/Q4p(2.0)
gamma=1.0/(11.0-2.0)-float(mp.re(-beta*wp(mp.mpf(zr(11.0))-z_h)))    # da Pi_h'=1/(r-2m)=-beta P+gamma
def Pi_h(x):
    z=mp.mpf(zr(x)); z0=mp.mpf(zr(r0)); p=lambda zz: beta*wzet(zz-z_h)+gamma*zz
    return float(mp.re(p(z)-p(z0)))
def Pih_dir(x): return quad(lambda t:1/((t-2.0)*np.sqrt(Q4(t))), r0, x, limit=300)[0]
print("Pi_h = beta*zeta(z-z_h)+gamma*z,  beta=%.6f, z_h=%.6fi, gamma=%.6f"%(beta,float(z_h.imag),gamma))
for x in [11.5,11.0,10.0,9.5,9.0]:
    print(f"  r={x}: dir={Pih_dir(x):+.8f} expl={Pi_h(x):+.8f} diff={abs(Pih_dir(x)-Pi_h(x)):.1e}")
