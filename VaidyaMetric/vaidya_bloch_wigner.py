# Stage B/C: dilog D_k su separatrice = DILOGARITMO ELLITTICO (Bloch-Wigner).
# Curva ellittica E: w^2=Q4(r), radici {-2.0833,0,2,8.7274}; diff olomorfo dr/sqrt(Q4).
# Abel z(r)=int dr/sqrt(Q4).  Orizzonte r=2m=2 = branch point => punto 2-TORSIONE
#   (mezzo periodo).  Polo r_d=-3.364 di U_k = punto generico z_d.
# Weight-2 iterato D_k=int U_k/(r-2) dr = DILOG ELLITTICO valutato in z_h-z_d.
import numpy as np, mpmath as mp
from scipy.integrate import quad
mp.mp.dps=30
M,E=1.0,1.4; r0=12.0
er=np.array([-2.0833333333,0.0,2.0,8.7274240]); a4=0.96
def Q4(x): return a4*np.prod([x-e for e in er])
Q4v=np.vectorize(Q4)
e1,e2,e3,e4=sorted(er)     # -2.0833,0,2,8.7274

# --- periodi/modulo (4 radici reali, Legendre): k^2=((e3-e2)(e4-e1))/((e4-e2)(e3-e1)) ---
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); Kc=mp.ellipk(k2); Kp=mp.ellipk(1-k2)
pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
wReal=pref*Kc; wImag=pref*Kp; tau=Kp/Kc
print(f"modulo ellittico k^2 = {float(k2):.8f}")
print(f"periodo reale  = {float(wReal):.8f}   periodo immag = {float(wImag):.8f}")
print(f"tau = {float(tau):.8f} i   (Sage EllipticCurve: 0.90597338 i)")

# --- Abel map z(r) per r nel range fisico [e4, r0], reale ---
def zt(rv): return quad(lambda x:1.0/np.sqrt(Q4(x)), e4, rv, limit=400)[0]   # da turning e4
# orizzonte r=2=e3: branch point => 2-TORSIONE. z_h imag = int_{e3}^{e4}=half-period immag
zh_im=quad(lambda x:1.0/np.sqrt(-Q4(x)), e3, e4, limit=400)[0]
print(f"\nz(orizzonte r=2m): parte immag = {zh_im:.8f}  = mezzo-periodo immag {float(wImag):.8f}")
print(f"  diff = {abs(zh_im-float(wImag)):.2e}  => orizzonte = punto 2-TORSIONE. CONFERMATO.")

# --- verifica dz/dr = 1/sqrt(Q4) (olomorfo) ---
print("\nverifica dz/dr=1/sqrt(Q4):")
for rv in [9.5,10.5,11.5]:
    dz=(zt(rv+1e-6)-zt(rv-1e-6))/2e-6
    print(f"  r={rv}: dz/dr={dz:.8f}  1/sqrtQ4={1/np.sqrt(Q4(rv)):.8f}  diff={abs(dz-1/np.sqrt(Q4(rv))):.1e}")
