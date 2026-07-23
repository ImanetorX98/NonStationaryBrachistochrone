# ASIMMETRIA accrescimento/evaporazione della correzione adiabatica delta phi_V.
# Accrescimento (mdot>0): Vaidya INGOING, clock = tempo AVANZATO v = E U_3 + r_*.
# Evaporazione (mdot<0): Vaidya OUTGOING, clock = tempo RITARDATO u = E U_3 - r_*.
#   (r_* = tortoise = (r-r0) + 2m ln((r-2m)/(r0-2m)),  v-u = 2 r_*)
# Stessa riduzione dM F (dipende solo da S congelata); cambia solo il clock.
# => delta phi/mdot|_accr = A + B ,  delta phi/mdot|_evap = A - B ,
#    A = int dM F * E U_3 dr   (parte SIMMETRICA, polilog all'INFINITO)
#    B = int dM F * r_*   dr   (parte ANTISIMMETRICA, dilog d'ORIZZONTE)
# ASIMMETRIA NETTA (a |mdot| uguale) = (accr)-(evap) = 2B: portata SOLO dall'orizzonte.
import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

M,E,J=1.0,1.4,2.5; r0=12.0
r=sp.symbols('r',positive=True); Ms,Es,Js=sp.symbols('M E J',positive=True)
Emu=(Es**2-1)*r+2*Ms
Ssym=sp.expand(r*(r-2*Ms)*Emu*(r**2*(r-2*Ms)-Js**2*Emu))
sub={Ms:1,Es:sp.Rational(7,5),Js:sp.Rational(5,2)}
Sn=sp.lambdify(r,Ssym.subs(sub),'numpy'); sq=lambda x:np.sqrt(Sn(x))
K=Js*Emu
dMF=(sp.diff(K,Ms)/sp.sqrt(Ssym)-K*sp.diff(Ssym,Ms)/(2*Ssym**sp.Rational(3,2)))
dMFn=sp.lambdify(r,dMF.subs(sub),'numpy')

rs=np.linspace(2.01,r0,20000); vv=Sn(rs); idx=np.where(np.diff(np.sign(vv)))[0]
rmin=max(brentq(Sn,rs[i],rs[i+1]) for i in idx if rs[i]>2.0); xf=rmin+0.4
def U3(x): return quad(lambda t:t**3/sq(t),r0,x,limit=200)[0]
def rstar(x): return (x-r0)+2*M*(np.log(x-2*M)-np.log(r0-2*M))
def v_adv(x): return E*U3(x)+rstar(x)     # accrescimento
def u_ret(x): return E*U3(x)-rstar(x)     # evaporazione

# --- parti A (simmetrica) e B (antisimmetrica=orizzonte) ---
A=quad(lambda x:dMFn(x)*E*U3(x), r0,xf,limit=150)[0]
B=quad(lambda x:dMFn(x)*rstar(x),r0,xf,limit=150)[0]
# --- diretti coi due clock ---
accr=quad(lambda x:dMFn(x)*v_adv(x),r0,xf,limit=150)[0]
evap=quad(lambda x:dMFn(x)*u_ret(x),r0,xf,limit=150)[0]

print("=== delta phi_V / mdot  (fattore comune mdot) ===")
print(f"  A (simmetrica, polilog INF)  = {A:+.8f}")
print(f"  B (antisimm., dilog ORIZZ)   = {B:+.8f}")
print(f"  accr = A+B (clock v avanzato)= {A+B:+.8f}   diretto={accr:+.8f}  diff={abs(A+B-accr):.1e}")
print(f"  evap = A-B (clock u ritard.) = {A-B:+.8f}   diretto={evap:+.8f}  diff={abs(A-B-evap):.1e}")
print(f"\nASIMMETRIA NETTA (accr-evap, |mdot| uguale) = 2B = {2*B:+.8f}")
print("  => portata INTERAMENTE dal termine d'ORIZZONTE (tortoise r_*); la parte")
print("     E U_3 (polilog all'infinito) e' COMUNE e si cancella nell'asimmetria.")

# --- segno fisico: delta phi_accr = mdot(A+B)>0? , delta phi_evap = mdot(A-B) con mdot<0 ---
print("\n=== segno fisico (delta phi = mdot * [.]) ===")
print(f"  accrescimento mdot>0: delta phi = +|mdot|*({A+B:+.4f})")
print(f"  evaporazione  mdot<0: delta phi = -|mdot|*({A-B:+.4f})")
print(f"  somma (a |mdot| uguale) = |mdot|*[(A+B)-(A-B)] = |mdot|*2B = {2*B:+.4f}|mdot|")

# --- struttura di B: bulk (2a specie) + settore d'ORIZZONTE (dilog D_k) ---
# r_* = (r-r0)  +  2m ln((r-2m)/(r0-2m))
#        \___bulk___/     \____log d'orizzonte____/
B_bulk=quad(lambda x:dMFn(x)*(x-r0), r0,xf,limit=150)[0]
B_hor =quad(lambda x:dMFn(x)*2*M*(np.log(x-2*M)-np.log(r0-2*M)), r0,xf,limit=150)[0]
print("\n=== struttura dell'asimmetria B = B_bulk + B_hor ===")
print(f"  B_bulk (2a specie, r-r0)           = {B_bulk:+.8f}")
print(f"  B_hor  (log d'orizzonte -> dilog D_k)= {B_hor:+.8f}")
print(f"  somma = {B_bulk+B_hor:+.8f}   B diretto = {B:+.8f}   diff={abs(B_bulk+B_hor-B):.1e}")
# B_hor per parti: 2m[ A_m ln(r-2m) ]_{r0}^{xf} - 2m int A_m/(r-2m) dr,  int A_m/(r-2m) = sum c_k D_k
print("  => B_hor = 2m[A_m ln(r-2m)] - 2m*int A_m/(r-2m);  int A_m/(r-2m) = dilog d'ORIZZONTE D_k.")
print("     L'asimmetria accr/evap trascendente (log-orizzonte) E' il dilog D_k gia' derivato.")
