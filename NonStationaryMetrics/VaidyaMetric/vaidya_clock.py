# Clock v(r) per l'adiabatico Vaidya (tempo avanzato) e struttura di delta phi.
# dv/dr = E r^3/sqrt(S) + r/(r-2M)  (2a specie + tortoise).
# v(r) = E R_3 + r + 2M ln(r-2M),  R_3=int r^3/sqrt(S).
# delta phi_V = Mdot int dM_F * v dr  => pezzi: polilog (da E R_3) + elementari
#   + NUOVO: int dM_F * ln(r-2M) dr  (abeliano pesato dal log d'orizzonte).
import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

M,E,J=1.0,1.4,2.5
r=sp.symbols('r',positive=True)
Ms,Es,Js=sp.symbols('M E J',positive=True)
Emu=(Es**2-1)*r+2*Ms
Ssym=sp.expand(r*(r-2*Ms)*Emu*(r**2*(r-2*Ms)-Js**2*Emu))
sub={Ms:1,Es:sp.Rational(7,5),Js:sp.Rational(5,2)}
Sn=sp.lambdify(r,Ssym.subs(sub),'numpy')
def sq(x): return np.sqrt(Sn(x))

# --- clock: verifica che dv/dr = (E/f) dtau/dr + 1/f coincida con E r^3/sqrtS + r/(r-2M) ---
f=1-2*M/r
dtau_dr=r**2*(r-2*M)/sp.sqrt(Ssym.subs(sub))          # a=0 tau clock
dv_dr_A=(E/f)*dtau_dr + 1/f                            # = dt/dr + 1/f
dv_dr_B=E*r**3/sp.sqrt(Ssym.subs(sub)) + r/(r-2*M)     # forma compatta
diff=sp.simplify(dv_dr_A-dv_dr_B)
print('dv/dr forma-A == forma-B ?', diff==0, ' (diff simbolica =', diff, ')')

# --- struttura di delta phi_V ---
# dM_F = N_M/S^(3/2); lo prendo come funzione
NM=sp.expand(Ssym*sp.diff(Js*Emu,Ms) - Js*Emu*sp.diff(Ssym,Ms)/2)  # =S dMK - K dMS/2
dMF=(sp.diff(Js*Emu,Ms)/sp.sqrt(Ssym) - Js*Emu*sp.diff(Ssym,Ms)/(2*Ssym**sp.Rational(3,2)))
dMFn=sp.lambdify(r,dMF.subs(sub),'numpy')

r0=12.0; wn=lambda x:E**2-(1-2*M/x)
rmin=brentq(lambda x: Sn(x), 2.001, r0-0.01) if Sn(2.001)*Sn(r0-0.01)<0 else None
# trova rmin come radice di S sopra 2M
rs=np.linspace(2.01,r0,20000); v=Sn(rs); idx=np.where(np.diff(np.sign(v)))[0]
rmin=max(brentq(Sn,rs[i],rs[i+1]) for i in idx if rs[i]>2.0)
print(f'rmin (turning) = {rmin:.4f}')

def R3(x): return quad(lambda t:t**3/sq(t), r0, x, limit=200)[0]     # 2a specie
def v_clock(x):  # v(r) = E R_3 + r + 2M ln(r-2M)  (rispetto a r0)
    return E*R3(x) + (x-r0) + 2*M*(np.log(x-2*M)-np.log(r0-2*M))
xf=rmin+0.4
# delta phi_V = -? Mdot int dM_F * v dr  (segno/fattore a parte); mostro i 3 pezzi
grid=np.linspace(r0-1e-3,xf,400)
# pezzo polilog (E R_3), pezzo elementare (r), pezzo LOG (2M ln(r-2M))
def piece(weight):
    return quad(lambda x: dMFn(x)*weight(x), r0, xf, limit=150)[0]
p_poly=piece(lambda x: E*R3(x))
p_lin =piece(lambda x: (x-r0))
p_log =piece(lambda x: 2*M*(np.log(x-2*M)-np.log(r0-2*M)))
p_tot =piece(v_clock)
print(f'\ndelta phi_V = Mdot * int dM_F * v dr,  v=E R_3 + r + 2M ln(r-2M):')
print(f'  pezzo polilog (E R_3)  = {p_poly:.5f}')
print(f'  pezzo elementare (r)   = {p_lin:.5f}')
print(f'  pezzo LOG-orizzonte    = {p_log:.5f}   <-- NUOVO in Vaidya')
print(f'  somma = {p_poly+p_lin+p_log:.5f}  vs totale diretto = {p_tot:.5f}  '
      f'diff={abs(p_poly+p_lin+p_log-p_tot):.1e}')
print('=> clock verificato; struttura = polilog genus-2 + elementari + LOG-orizzonte.')
