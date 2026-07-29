# ASIMMETRIA accrescimento/evaporazione della correzione adiabatica delta phi_V.
# Accrescimento (mdot>0): Vaidya INGOING  ds^2=-f dv^2+2 dv dr, clock AVANZATO  v = t + r_*.
# Evaporazione (mdot<0): Vaidya OUTGOING ds^2=-f du^2-2 du dr, clock RITARDATO u = t - r_*.
# DERIVAZIONE del clock uscente (referee Issue 4, non un flip ad-hoc): u,v sono definiti dal MEDESIMO
#   tempo di Killing t via u=t-r_*, v=t+r_* (r_*=tortoise). Quindi lungo QUALSIASI traiettoria
#   v - u = 2 r_*  =>  du/dr = dv/dr - 2 dr_*/dr. Con la brachistocrona congelata dv/dr=E r^3/sqrt(S)+r/(r-2m)
#   (v-branch, derivato dalla metrica ingoing) e dr_*/dr=r/(r-2m):
#     du/dr = E r^3/sqrt(S) - r/(r-2m)   =>   u = E U_3 - r_*   (DERIVATO dalla metrica outgoing).
#   (r_* = (r-r0) + 2m ln((r-2m)/(r0-2m)),  dr_*/dr=r/(r-2m).)
# Stessa riduzione dM F (dipende solo da S congelata: ingoing/outgoing Schwarzschild = stessa geometria);
# cambiano SOLO il clock (v vs u) e il segno di mdot (accrescimento m su' / evaporazione m giu').
# => delta phi/mdot|_accr = A + B ,  delta phi/mdot|_evap = A - B  (COEFFICIENTI, /mdot con segno),
#    A = int dM F * E U_3 dr   (parte all'INFINITO, polilog)
#    B = int dM F * r_*   dr   (parte d'ORIZZONTE, dilog tortoise)
# PHYSICS (delta phi = mdot*[.]; mdot>0 accr, mdot<0 evap): delta phi_accr=+|mdot|(A+B),
#    delta phi_evap=-|mdot|(A-B). Quindi (a |mdot| uguale):
#    DIFFERENZA fisica = accr - evap = |mdot|*2A   (portata da A -- NON da B!)
#    SOMMA fisica      = accr + evap = |mdot|*2B   = FALLIMENTO DELL'ANTISIMMETRIA (portata da B, orizzonte)
# (referee Issue 4 FIX: B controlla il fallimento dell'antisimmetria, non la differenza fisica.)
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
def v_adv(x): return E*U3(x)+rstar(x)     # accrescimento (ingoing, derivato metrica ingoing)
def u_ret(x): return E*U3(x)-rstar(x)     # evaporazione (outgoing, DERIVATO: u=v-2r_*)

# --- verifica derivazione clock outgoing: du/dr = dv/dr - 2 dr_*/dr (referee Issue 4) ---
def dvdr(x): return E*x**3/sq(x)+x/(x-2*M)      # v-branch (metrica ingoing)
def dudr(x): return E*x**3/sq(x)-x/(x-2*M)      # u-branch (derivato)
def drstardr(x): return x/(x-2*M)
print("verifica clock outgoing derivato: du/dr == dv/dr - 2 dr_*/dr :",
      f"{max(abs(dudr(xx)-(dvdr(xx)-2*drstardr(xx))) for xx in [6.0,8.0,10.0]):.1e}")

# --- parti A (simmetrica) e B (antisimmetrica=orizzonte) ---
A=quad(lambda x:dMFn(x)*E*U3(x), r0,xf,limit=150)[0]
B=quad(lambda x:dMFn(x)*rstar(x),r0,xf,limit=150)[0]
# --- diretti coi due clock ---
accr=quad(lambda x:dMFn(x)*v_adv(x),r0,xf,limit=150)[0]
evap=quad(lambda x:dMFn(x)*u_ret(x),r0,xf,limit=150)[0]

print("=== delta phi_V : coefficienti (/mdot con segno) e quantita' FISICHE ===")
print(f"  A (parte INFINITO, polilog)  = {A:+.8f}")
print(f"  B (parte ORIZZONTE, dilog)   = {B:+.8f}")
print(f"  coeff accr = A+B (clock v)   = {A+B:+.8f}   diretto={accr:+.8f}  diff={abs(A+B-accr):.1e}")
print(f"  coeff evap = A-B (clock u)   = {A-B:+.8f}   diretto={evap:+.8f}  diff={abs(A-B-evap):.1e}")
# --- FISICA: delta phi = mdot*[.], mdot>0 accr, mdot<0 evap => segni ---
dphi_accr=(A+B)     # /|mdot|, mdot>0
dphi_evap=-(A-B)    # /|mdot|, mdot<0 (delta phi_evap = mdot(A-B) = -|mdot|(A-B))
print("\n=== FISICA (delta phi / |mdot|) ===")
print(f"  delta phi_accr/|mdot| = +(A+B) = {dphi_accr:+.6f}")
print(f"  delta phi_evap/|mdot| = -(A-B) = {dphi_evap:+.6f}")
print(f"  DIFFERENZA  accr-evap = 2A = {dphi_accr-dphi_evap:+.6f}   (portata da A, l'INFINITO -- NON da B!)")
print(f"  SOMMA       accr+evap = 2B = {dphi_accr+dphi_evap:+.6f}   = FALLIMENTO DELL'ANTISIMMETRIA (B, orizzonte)")
print("  (referee Issue 4 FIX: la vecchia dicitura 'asimmetria=2B' confondeva la differenza dei coeff")
print("   (A+B)-(A-B)=2B con la differenza FISICA accr-evap=2A; B misura la deviazione dall'antisimmetria.)")

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
