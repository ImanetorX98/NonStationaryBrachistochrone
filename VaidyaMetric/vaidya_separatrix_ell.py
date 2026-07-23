# SEPARATRICE |J|=Jc: riduzione ELLITTICA (genere 1) di U_k e del dilog D_k.
# Stage A: su Jc, sqrt(S)=(r-r_d) sqrt(Q4).  U_k = int r^k/[(r-r_d) sqrt(Q4)] dr.
#   r^k/(r-r_d) = poly_{k-1}(r) + r_d^k/(r-r_d)  =>
#   U_k = (ellittico 1a/2a specie) + r_d^k * Pi(r_d)   [3a specie ellittica, polo r_d].
# Verifica: U_k ellittico == U_k diretto (genus-2 quad @Jc)  a ~1e-10.
import numpy as np, sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq

M,E=1.0,1.4; r0=12.0
# --- Jc esatto (radice doppia bracket) ---
r=sp.symbols('r'); Js=sp.symbols('J',positive=True)
B=r**2*(r-2*M)-Js**2*((E**2-1)*r+2*M)
Jc=float(min([s for s in sp.solve(sp.Eq(sp.discriminant(sp.Poly(B,r),r),0),Js) if s.is_real and s>0],key=abs))
# r_d = radice doppia
Bc=np.array([float(c) for c in sp.Poly(B.subs(Js,Jc),r).all_coeffs()])
rts=np.roots(Bc); r_d=float(sorted(rts,key=lambda z:abs(z-rts[np.argmin([abs(z1-z2) for i,z1 in enumerate(rts) for z2 in rts if not (z1 is z2)][:1])]))[0].real) if False else None
# semplice: r_d = la radice ripetuta
rr=sorted([z.real for z in rts if abs(z.imag)<1e-6])
# trova la coppia piu vicina
pairs=[(i,j,abs(rr[i]-rr[j])) for i in range(len(rr)) for j in range(i+1,len(rr))]
i,j,_=min(pairs,key=lambda t:t[2]); r_d=0.5*(rr[i]+rr[j])
print(f'Jc={Jc:.8f}  r_d(doppia)={r_d:.8f}')

def DE(r): return (E**2-1)*r+2*M
def Sn(r,Jv):  return r*(r-2*M)*DE(r)*(r**2*(r-2*M)-Jv**2*DE(r))
# Q4 = S/(r-r_d)^2  via divisione polinomiale NUMERICA (r_d float, radice doppia approx)
Ssep=sp.expand(r*(r-2*M)*((E**2-1)*r+2*M)*(r**2*(r-2*M)-Jc**2*((E**2-1)*r+2*M)))
Scoef=[float(c) for c in sp.Poly(Ssep,r).all_coeffs()]        # deg6, decrescente
den=np.polynomial.polynomial.polypow([-r_d,1],2)              # (r-r_d)^2, crescente
q,rem=np.polynomial.polynomial.polydiv(Scoef[::-1],den)       # crescente
print('resto divisione (deve ~0):',np.max(np.abs(rem)))
Q4c=q[::-1]                                                   # Q4 coeff decrescente
Q4=lambda x:np.polyval(Q4c,x)
Q4roots=sorted([z.real for z in np.roots(Q4c) if abs(z.imag)<1e-6])
print('Q4 radici:',[round(z,5) for z in Q4roots])

# --- diretto genus-2 @Jc ---
sqS=lambda x:np.sqrt(Sn(x,Jc))
def Uk_dir(x,k): return quad(lambda t:t**k/sqS(t),r0,x,limit=300)[0]
# --- ellittico: sqrt(S)=(r-r_d) sqrt(Q4), stesso segno su [xf,r0] (r>r_d>0? r_d<0 qui) ---
sgn=np.sign((r0-r_d))  # r_d<0 => r-r_d>0 su tutta la traiettoria
def Uk_ell(x,k): return quad(lambda t:t**k/((t-r_d)*np.sqrt(Q4(t))),r0,x,limit=300)[0]

# estremo dentro (turning fisico su separatrice = radice doppia? no, turning = 8.727)
# la brachistocrona fisica su Jc ha turning al root 8.727 (>r0=12? no 8.727<12) -> parte da r0=12 verso dentro
rt=max([z for z in Q4roots if z<r0])
xf=rt+0.3
print(f'turning fisico (root Q4) = {rt:.5f}  xf={xf:.5f}')
print('\nStage A  U_k: ellittico vs diretto genus-2 @Jc')
for k in range(5):
    ud=Uk_dir(xf,k); ue=Uk_ell(xf,k)
    print(f'  k={k}: dir={ud:+.8f}  ell={ue:+.8f}  diff={abs(ud-ue):.2e}')
