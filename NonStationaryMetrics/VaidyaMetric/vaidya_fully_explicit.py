# FORMA COMPLETAMENTE ESPLICITA di delta phi_V/Mdot (Vaidya adiabatico).
# Ogni integrale espanso in building block NOMINATI:
#   U_k      = int r^k/sqrt(S) dr                     (2a specie, polilog genus-2)
#   W_3k     = int (U_3 dU_k - U_k dU_3)              (polilog peso-2, genus-2)
#   D_k      = int U_k/(r-2M) dr                      (dilog ORIZZONTE, peso-2)
#   L_2M     = int Acal^M r/((r-2M) sqrt(S)) dr       (3a specie, polo orizzonte)   [nota sotto]
#   I_poly   = int Acal^M r^3/S dr                    (ELEMENTARE: log/fratti razionali)
# Coefficienti c_k^M, a_j^M ALGEBRICI (da riduzione dM F = d(Acal^M/sqrtS)+sum c_k^M r^k/sqrtS).
#
# ASSEMBLAGGIO (per parti + riduzione, esatto):
#   delta phi_V/Mdot = A_M v
#       - E*I_poly - (E/2) sum_k c_k^M (U_k U_3 - W_3k)
#       - sum_j a_j^M U_j - sum_k c_k^M (r U_k - U_{k+1}) - 2M*Lh - 2M sum_k c_k^M D_k
# con A_M = Acal^M/sqrtS + sum_k c_k^M U_k,  v = E U_3 + (r-r0) + 2M ln((r-2M)/(r0-2M)).
import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

M,E,J=1.0,1.4,2.5; r0=12.0
r=sp.symbols('r',positive=True); Ms,Es,Js=sp.symbols('M E J',positive=True)
Emu=(Es**2-1)*r+2*Ms
Ssym=sp.expand(r*(r-2*Ms)*Emu*(r**2*(r-2*Ms)-Js**2*Emu))
sub={Ms:1,Es:sp.Rational(7,5),Js:sp.Rational(5,2)}
Sn=sp.lambdify(r,Ssym.subs(sub),'numpy'); sq=lambda x:np.sqrt(Sn(x))

# ---- riduzione: dM F = d(Acal/sqrtS) + sum c_k r^k/sqrtS ----
K=Js*Emu; NM=sp.expand(Ssym*sp.diff(K,Ms)-K*sp.diff(Ssym,Ms)/2)
Ac=[sp.Symbol(f'A{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
Acal=sum(Ac[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
sol=sp.solve(sp.Poly(sp.expand(2*NM-(2*Ssym*sp.diff(Acal,r)-Acal*sp.diff(Ssym,r)+2*Ssym*Mp)),r).all_coeffs(),Ac+ck,dict=True)[0]
cM=[float(sp.simplify(sol[ck[i]]).subs(sub)) for i in range(5)]
aM=[float(sp.simplify(sol.get(Ac[i],0)).subs(sub)) for i in range(6)]  # coeff di Acal^M, deg5
def Acaln(x): return sum(aM[i]*x**i for i in range(6))
dMF=(sp.diff(K,Ms)/sp.sqrt(Ssym)-K*sp.diff(Ssym,Ms)/(2*Ssym**sp.Rational(3,2)))
dMFn=sp.lambdify(r,dMF.subs(sub),'numpy')
print('c_k^M =',[round(c,5) for c in cM])
print('a_j^M =',[round(a,5) for a in aM])

# ---- estremo di integrazione (dentro, oltre il turning) ----
rs=np.linspace(2.01,r0,20000); vv=Sn(rs); idx=np.where(np.diff(np.sign(vv)))[0]
rmin=max(brentq(Sn,rs[i],rs[i+1]) for i in idx if rs[i]>2.0); xf=rmin+0.4

# ---- building block ----
def Uk(x,k): return quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def W3k(x,k): # int (U_3 dU_k - U_k dU_3) = int (U_3 r^k - U_k r^3)/sqrtS dr
    return quad(lambda t:(Uk(t,3)*t**k-Uk(t,k)*t**3)/sq(t),r0,x,limit=150)[0]
def Dk(x,k): return quad(lambda t:Uk(t,k)/(t-2*M),r0,x,limit=150)[0]
def Lh(x): return quad(lambda t:Acaln(t)*t/((t-2*M)*sq(t)),r0,x,limit=150)[0]
def Ipoly(x): return quad(lambda t:Acaln(t)*t**3/Sn(t),r0,x,limit=150)[0]
def AM(x): return Acaln(x)/sq(x)-Acaln(r0)/sq(r0)+sum(cM[k]*Uk(x,k) for k in range(5))
def vclock(x): return E*Uk(x,3)+(x-r0)+2*M*(np.log(x-2*M)-np.log(r0-2*M))

# ---- DIRETTO ----
dphi_dir=quad(lambda x:dMFn(x)*vclock(x),r0,xf,limit=150)[0]

# ---- ASSEMBLAGGIO COMPLETAMENTE ESPLICITO ----
# A_M = Acal/sqrtS - alpha0 + sum c_k U_k,  alpha0 = Acal(r0)/sqrtS(r0) costante.
x=xf
a0 = Acaln(r0)/np.sqrt(Sn(r0))
bdy = AM(x)*vclock(x)
# E * int A_M r^3/sqrtS = E[ I_poly - a0 U_3 + sum c_k (1/2)(U_k U_3 - W_3k) ]
p2  = E*( Ipoly(x) - a0*Uk(x,3) + 0.5*sum(cM[k]*(Uk(x,k)*Uk(x,3)-W3k(x,k)) for k in range(5)) )
# int A_M r/(r-2M) = Lh - a0 (v - E U_3) + sum c_k[(r U_k - U_{k+1}) + 2M D_k]
p3  = ( Lh(x) - a0*(vclock(x)-E*Uk(x,3))
      + sum(cM[k]*(x*Uk(x,k)-Uk(x,k+1)) for k in range(5))
      + 2*M*sum(cM[k]*Dk(x,k) for k in range(5)) )
dphi_expl = bdy - p2 - p3

print(f'\ndelta phi_V/Mdot  DIRETTO        = {dphi_dir:.10f}')
print(f'                  ESPLICITO full = {dphi_expl:.10f}')
print(f'  diff = {abs(dphi_dir-dphi_expl):.2e}   <-- identita (coeff algebrici dati)')

# ---- valori dei building block (per il paper) ----
print('\nBuilding block @ xf=%.4f:'%xf)
print('  U_k   =',[round(Uk(x,k),4) for k in range(6)])
print('  W_3k  =',[round(W3k(x,k),4) for k in range(5)])
print('  D_k   =',[round(Dk(x,k),4) for k in range(5)])
print('  L_2M  =',round(Lh(x),4),'   I_poly =',round(Ipoly(x),4))
print('  A_M   =',round(AM(x),4),'   v =',round(vclock(x),4))
