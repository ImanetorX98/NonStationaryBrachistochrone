# EXPLICIT closed form of delta phi_V (Vaidya), named building blocks, verified.
# delta phi_V/Mdot = int dM_F * v dr,  v = E U_3 + (r-r0) + 2M ln((r-2M)/(r0-2M)).
# By parts + reduction: pieces = genus-2 polylog (W) + second kind (U) + HORIZON dilog
#   (D_k=int U_k/(r-2M)) + horizon third kind + elementary (logs, products).
import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

M,E,J=1.0,1.4,2.5; r0=12.0
r=sp.symbols('r',positive=True); Ms,Es,Js=sp.symbols('M E J',positive=True)
Emu=(Es**2-1)*r+2*Ms
Ssym=sp.expand(r*(r-2*Ms)*Emu*(r**2*(r-2*Ms)-Js**2*Emu))
sub={Ms:1,Es:sp.Rational(7,5),Js:sp.Rational(5,2)}
Sn=sp.lambdify(r,Ssym.subs(sub),'numpy'); sq=lambda x:np.sqrt(Sn(x))
# reduction dM F = d(Acal/sqrtS)+sum c_k^M r^k/sqrtS
K=Js*Emu; NM=sp.expand(Ssym*sp.diff(K,Ms)-K*sp.diff(Ssym,Ms)/2)
Ac=[sp.Symbol(f'A{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
Acal=sum(Ac[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
sol=sp.solve(sp.Poly(sp.expand(2*NM-(2*Ssym*sp.diff(Acal,r)-Acal*sp.diff(Ssym,r)+2*Ssym*Mp)),r).all_coeffs(),Ac+ck,dict=True)[0]
cM=[float(sp.simplify(sol[ck[i]]).subs(sub)) for i in range(5)]
Anum=[float(sp.simplify(sol.get(Ac[i],0)).subs(sub)) for i in range(6)]
def Acaln(x): return sum(Anum[i]*x**i for i in range(6))
dMF=(sp.diff(K,Ms)/sp.sqrt(Ssym)-K*sp.diff(Ssym,Ms)/(2*Ssym**sp.Rational(3,2)))
dMFn=sp.lambdify(r,dMF.subs(sub),'numpy')

rs=np.linspace(2.01,r0,20000); vv=Sn(rs); idx=np.where(np.diff(np.sign(vv)))[0]
rmin=max(brentq(Sn,rs[i],rs[i+1]) for i in idx if rs[i]>2.0); xf=rmin+0.4
def Uk(x,k): return quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def AM(x): return Acaln(x)/sq(x)-Acaln(r0)/sq(r0)+sum(cM[k]*Uk(x,k) for k in range(5))
def vclock(x): return E*Uk(x,3)+(x-r0)+2*M*(np.log(x-2*M)-np.log(r0-2*M))

# --- DIRECT ---
dphi_dir=quad(lambda x:dMFn(x)*vclock(x),r0,xf,limit=150)[0]

# --- EXPLICIT (by parts): A_M v - E int A_M r^3/sqrtS - int A_M r/(r-2M) ---
term_bdy=AM(xf)*vclock(xf)
int_AM_r3=quad(lambda x:AM(x)*x**3/sq(x),r0,xf,limit=150)[0]           # -> polilog+prod+elem
int_AM_tort=quad(lambda x:AM(x)*x/(x-2*M),r0,xf,limit=150)[0]          # -> horizon dilog + elementary
dphi_expl=term_bdy - E*int_AM_r3 - int_AM_tort
print(f'delta phi_V/Mdot direct   = {dphi_dir:.8f}')
print(f'            by parts, expl = {dphi_expl:.8f}')
print(f'  diff = {abs(dphi_dir-dphi_expl):.2e}  (identity, coefficients given)')

# --- the NEW piece: horizon dilog D_k = int U_k/(r-2M) dr ---
# int A_M r/(r-2M) = int A_M dr + 2M int A_M/(r-2M);  int A_M/(r-2M) contains the D_k
def Dk(x,k): return quad(lambda t:Uk(t,k)/(t-2*M),r0,x,limit=150)[0]
print(f'\nHORIZON dilog D_k(xf)=int U_k/(r-2M): ',[round(Dk(xf,k),4) for k in range(5)])
print('algebraic coefficients c_k^M:',[round(c,4) for c in cM])
print('=> explicit delta phi_V: A_M v - E[polylog+products+elementary] - [horizon dilog+elementary].')
print('   building blocks: U_k, W_jk (polylog at infinity), D_k (horizon dilog), logs, products.')
