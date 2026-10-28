# SEPARATRIX |J|=Jc: ELLIPTIC reduction (genus 1) of U_k and of the dilog D_k.
# Stage A: su Jc, sqrt(S)=(r-r_d) sqrt(Q4).  U_k = int r^k/[(r-r_d) sqrt(Q4)] dr.
#   r^k/(r-r_d) = poly_{k-1}(r) + r_d^k/(r-r_d)  =>
#   U_k = (elliptic, first/second kind) + r_d^k * Pi(r_d)   [elliptic third kind, pole r_d].
# Check: elliptic U_k == direct U_k (genus-2 quadrature at Jc)  to ~1e-10.
import numpy as np, sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq

M,E=1.0,1.4; r0=12.0
# --- exact Jc (bracketing the double root) ---
r=sp.symbols('r'); Js=sp.symbols('J',positive=True)
B=r**2*(r-2*M)-Js**2*((E**2-1)*r+2*M)
Jc=float(min([s for s in sp.solve(sp.Eq(sp.discriminant(sp.Poly(B,r),r),0),Js) if s.is_real and s>0],key=abs))
# r_d = double root
Bc=np.array([float(c) for c in sp.Poly(B.subs(Js,Jc),r).all_coeffs()])
rts=np.roots(Bc); r_d=float(sorted(rts,key=lambda z:abs(z-rts[np.argmin([abs(z1-z2) for i,z1 in enumerate(rts) for z2 in rts if not (z1 is z2)][:1])]))[0].real) if False else None
# simple: r_d = the repeated root
rr=sorted([z.real for z in rts if abs(z.imag)<1e-6])
# find the closest pair
pairs=[(i,j,abs(rr[i]-rr[j])) for i in range(len(rr)) for j in range(i+1,len(rr))]
i,j,_=min(pairs,key=lambda t:t[2]); r_d=0.5*(rr[i]+rr[j])
print(f'Jc={Jc:.8f}  r_d(double)={r_d:.8f}')

def DE(r): return (E**2-1)*r+2*M
def Sn(r,Jv):  return r*(r-2*M)*DE(r)*(r**2*(r-2*M)-Jv**2*DE(r))
# Q4 = S/(r-r_d)^2  by NUMERICAL polynomial division (r_d float, approximate double root)
Ssep=sp.expand(r*(r-2*M)*((E**2-1)*r+2*M)*(r**2*(r-2*M)-Jc**2*((E**2-1)*r+2*M)))
Scoef=[float(c) for c in sp.Poly(Ssep,r).all_coeffs()]        # degree 6, decreasing powers
den=np.polynomial.polynomial.polypow([-r_d,1],2)              # (r-r_d)^2, increasing powers
q,rem=np.polynomial.polynomial.polydiv(Scoef[::-1],den)       # increasing powers
print('division remainder (must be ~0):',np.max(np.abs(rem)))
Q4c=q[::-1]                                                   # Q4 coefficients, decreasing powers
Q4=lambda x:np.polyval(Q4c,x)
Q4roots=sorted([z.real for z in np.roots(Q4c) if abs(z.imag)<1e-6])
print('Q4 roots:',[round(z,5) for z in Q4roots])

# --- direct genus 2 at Jc ---
sqS=lambda x:np.sqrt(Sn(x,Jc))
def Uk_dir(x,k): return quad(lambda t:t**k/sqS(t),r0,x,limit=300)[0]
# --- elliptic: sqrt(S)=(r-r_d) sqrt(Q4), same sign on [xf,r0] (r>r_d>0? here r_d<0) ---
sgn=np.sign((r0-r_d))  # r_d<0 => r-r_d>0 along the whole trajectory
def Uk_ell(x,k): return quad(lambda t:t**k/((t-r_d)*np.sqrt(Q4(t))),r0,x,limit=300)[0]

# inner endpoint (physical turning point on the separatrix = double root? no, turning point = 8.727)
# the physical brachistochrone at Jc turns at the root 8.727 (>r0=12? no, 8.727<12) -> it starts from r0=12 inwards
rt=max([z for z in Q4roots if z<r0])
xf=rt+0.3
print(f'physical turning point (root of Q4) = {rt:.5f}  xf={xf:.5f}')
print('\nStage A  U_k: elliptic against direct genus 2 at Jc')
for k in range(5):
    ud=Uk_dir(xf,k); ue=Uk_ell(xf,k)
    print(f'  k={k}: dir={ud:+.8f}  ell={ue:+.8f}  diff={abs(ud-ue):.2e}')
