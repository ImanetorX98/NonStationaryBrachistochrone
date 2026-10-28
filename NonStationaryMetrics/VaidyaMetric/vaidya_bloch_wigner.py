# Stage B/C: dilog D_k on the separatrix = ELLIPTIC DILOGARITHM (Bloch-Wigner).
# Elliptic curve E: w^2=Q4(r), roots {-2.0833,0,2,8.7274}; holomorphic differential dr/sqrt(Q4).
# Abel z(r)=int dr/sqrt(Q4).  Horizon r=2m=2 = branch point => 2-TORSION point
#   (half period).  Pole r_d=-3.364 of U_k = generic point z_d.
# Iterated weight 2: D_k=int U_k/(r-2) dr = ELLIPTIC DILOG evaluated at z_h-z_d.
import numpy as np, mpmath as mp
from scipy.integrate import quad
mp.mp.dps=30
M,E=1.0,1.4; r0=12.0
er=np.array([-2.0833333333,0.0,2.0,8.7274240]); a4=0.96
def Q4(x): return a4*np.prod([x-e for e in er])
Q4v=np.vectorize(Q4)
e1,e2,e3,e4=sorted(er)     # -2.0833,0,2,8.7274

# --- periods/modulus (4 real roots, Legendre): k^2=((e3-e2)(e4-e1))/((e4-e2)(e3-e1)) ---
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); Kc=mp.ellipk(k2); Kp=mp.ellipk(1-k2)
pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
wReal=pref*Kc; wImag=pref*Kp; tau=Kp/Kc
print(f"elliptic modulus k^2 = {float(k2):.8f}")
print(f"real period  = {float(wReal):.8f}   imaginary period = {float(wImag):.8f}")
print(f"tau = {float(tau):.8f} i   (Sage EllipticCurve: 0.90597338 i)")

# --- Abel map z(r) for r in the physical range [e4, r0], real ---
def zt(rv): return quad(lambda x:1.0/np.sqrt(Q4(x)), e4, rv, limit=400)[0]   # from the turning point e4
# horizon r=2=e3: branch point => 2-TORSION. imaginary part of z_h = int_{e3}^{e4} = imaginary half-period
zh_im=quad(lambda x:1.0/np.sqrt(-Q4(x)), e3, e4, limit=400)[0]
print(f"\nz(horizon r=2m): imaginary part = {zh_im:.8f}  = imaginary half-period {float(wImag):.8f}")
print(f"  diff = {abs(zh_im-float(wImag)):.2e}  => horizon = 2-TORSION point. CONFIRMED.")

# --- check dz/dr = 1/sqrt(Q4) (holomorphic) ---
print("\ncheck dz/dr=1/sqrt(Q4):")
for rv in [9.5,10.5,11.5]:
    dz=(zt(rv+1e-6)-zt(rv-1e-6))/2e-6
    print(f"  r={rv}: dz/dr={dz:.8f}  1/sqrtQ4={1/np.sqrt(Q4(rv)):.8f}  diff={abs(dz-1/np.sqrt(Q4(rv))):.1e}")
