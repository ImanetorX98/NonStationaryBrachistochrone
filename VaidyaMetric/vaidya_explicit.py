# Forma ESPLICITA chiusa di delta phi_V (Vaidya), building block nominati, verificata.
# delta phi_V/Mdot = int dM_F * v dr,  v = E U_3 + (r-r0) + 2M ln((r-2M)/(r0-2M)).
# Per parti + riduzione: pezzi = polilog genus-2 (W) + 2a specie (U) + dilog ORIZZONTE
#   (D_k=int U_k/(r-2M)) + 3a specie orizzonte + elementari (log, prodotti).
import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

M,E,J=1.0,1.4,2.5; r0=12.0
r=sp.symbols('r',positive=True); Ms,Es,Js=sp.symbols('M E J',positive=True)
Emu=(Es**2-1)*r+2*Ms
Ssym=sp.expand(r*(r-2*Ms)*Emu*(r**2*(r-2*Ms)-Js**2*Emu))
sub={Ms:1,Es:sp.Rational(7,5),Js:sp.Rational(5,2)}
Sn=sp.lambdify(r,Ssym.subs(sub),'numpy'); sq=lambda x:np.sqrt(Sn(x))
# riduzione dM F = d(Acal/sqrtS)+sum c_k^M r^k/sqrtS
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

# --- DIRETTO ---
dphi_dir=quad(lambda x:dMFn(x)*vclock(x),r0,xf,limit=150)[0]

# --- ESPLICITO (per parti): A_M v - E int A_M r^3/sqrtS - int A_M r/(r-2M) ---
term_bdy=AM(xf)*vclock(xf)
int_AM_r3=quad(lambda x:AM(x)*x**3/sq(x),r0,xf,limit=150)[0]           # -> polilog+prod+elem
int_AM_tort=quad(lambda x:AM(x)*x/(x-2*M),r0,xf,limit=150)[0]          # -> dilog orizzonte+elem
dphi_expl=term_bdy - E*int_AM_r3 - int_AM_tort
print(f'delta phi_V/Mdot diretto  = {dphi_dir:.8f}')
print(f'            per-parti espl = {dphi_expl:.8f}')
print(f'  diff = {abs(dphi_dir-dphi_expl):.2e}  (identita, coeff dati)')

# --- il pezzo NUOVO: dilog orizzonte D_k = int U_k/(r-2M) dr ---
# int A_M r/(r-2M) = int A_M dr + 2M int A_M/(r-2M);  int A_M/(r-2M) contiene i D_k
def Dk(x,k): return quad(lambda t:Uk(t,k)/(t-2*M),r0,x,limit=150)[0]
print(f'\ndilog ORIZZONTE D_k(xf)=int U_k/(r-2M): ',[round(Dk(xf,k),4) for k in range(5)])
print('coeff algebrici c_k^M:',[round(c,4) for c in cM])
print('=> delta phi_V esplicito: A_M v - E[polilog+prod+elem] - [dilog-orizzonte+elem].')
print('   building block: U_k, W_jk (polilog inf), D_k (dilog orizzonte), log, prodotti.')
