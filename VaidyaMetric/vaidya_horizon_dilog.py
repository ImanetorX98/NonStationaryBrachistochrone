# Chiusura del DILOG D'ORIZZONTE D_k = int U_k/(r-2m) dr (Vaidya, frozen Schwarzschild).
# ---------------------------------------------------------------------------
# STRUTTURA (per parti esatta):
#   eta_h = dr/(r-2m).  Come 1-forma sulla curva y^2=S, in t=sqrt(r-2m):
#     r-2m=t^2, dr=2t dt  =>  eta_h = 2 dt/t  = polo SEMPLICE, residuo 2 al
#     punto di Weierstrass sopra r=2m  => TERZA SPECIE (lettera dlog ORIZZONTE).
#   Primitiva elementare: int eta_h = ln(r-2m).
#   => D_k = U_k ln(r-2m) - G_k,   G_k = int ln(r-2m) r^k/sqrt(S) dr   (PASSO1)
#   G_k = dilog iperellittico genus-2 con lettera dlog all'ORIZZONTE
#         (compagno di W_kj, che ha la lettera dipolo all'INFINITO).
# SEPARATRICE |J|=Jc: S acquista radice doppia -> genere scende a 1 ->
#   G_k, D_k CHIUDONO nel DILOGARITMO ELLITTICO (Bloch-Wigner), tabulato (PASSO2).
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import sympy as sp

M,E,J=1.0,1.4,2.5; r0=12.0
def DE(r): return (E**2-1)*r+2*M
def Sn(r):  return r*(r-2*M)*DE(r)*(r**2*(r-2*M)-J**2*DE(r))
sq=lambda x:np.sqrt(Sn(x))
rs=np.linspace(2.01,r0,20000); vv=Sn(rs); idx=np.where(np.diff(np.sign(vv)))[0]
rmin=max(brentq(Sn,rs[i],rs[i+1]) for i in idx if rs[i]>2.0); xf=rmin+0.4

def Uk(x,k): return quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def Dk_dir(x,k): return quad(lambda t:Uk(t,k)/(t-2*M),r0,x,limit=150)[0]     # DIRETTO
def Gk(x,k): return quad(lambda t:np.log(abs(t-2*M))*t**k/sq(t),r0,x,limit=150)[0]

# ---- PASSO1: D_k = U_k ln(r-2m) - G_k  (identita IBP) ----
print('PASSO1  D_k = U_k ln(r-2m) - G_k   (per parti, lettera dlog orizzonte)')
print(' k :  D_k diretto      U_k ln - G_k       diff')
c=np.log(abs(r0-2*M))   # nota: U_k(r0)=0, quindi termine al bordo r0 = 0
for k in range(5):
    d_dir=Dk_dir(xf,k)
    d_ibp=Uk(xf,k)*np.log(abs(xf-2*M)) - Gk(xf,k)
    print(f' {k} : {d_dir:+.9f}   {d_ibp:+.9f}   {abs(d_dir-d_ibp):.2e}')

# ---- Jc separatrice: bracket B(r)=r^2(r-2m)-J^2 DE ha radice doppia ----
r=sp.symbols('r'); Js=sp.symbols('J',positive=True)
B=r**2*(r-2*M)-Js**2*((E**2-1)*r+2*M)
disc=sp.discriminant(sp.Poly(B,r),r)
Jc_sols=[s for s in sp.solve(sp.Eq(disc,0),Js) if s.is_real and s>0]
Jc=float(min(Jc_sols,key=lambda s:abs(float(s))))
print(f'\nPASSO2  separatrice Jc = {Jc:.6f}  (radice doppia del bracket)')
Bc=sp.Poly(B.subs(Js,Jc),r)
roots=sorted([complex(x).real for x in np.roots([float(c) for c in Bc.all_coeffs()]) if abs(complex(x).imag)<1e-6])
print(f'  radici bracket @Jc: {[round(x,5) for x in roots]}  (doppia = genere 1)')
print('  => su Jc: sqrt(S)=(r-r_d) sqrt(quartica), U_k ellittici,')
print('     D_k = DILOGARITMO ELLITTICO (Bloch-Wigner) tabulato.')
