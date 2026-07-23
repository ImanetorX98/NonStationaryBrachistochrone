# -*- coding: utf-8 -*-
# CHIUSURA parte OLOMORFA nella base canonica (TK tau, genus-2).
# Fisso i coeff di POLO ai g_i simbolici (parti principali); il residuo U_k-(polo) e'
# OLOMORFO -> risolvo i coeff di u1,u2 (+const) BEN CONDIZIONATO (non piu' fit degenere).
# Status atteso: come C0/Ce separatrice (principio a-periodi, no fit, dipendenti dai periodi).
from sage.all import QQ
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface as SageRS
from abelfunctions import RiemannTheta
import numpy as np, sympy as sp
from scipy.integrate import quad
import itertools

M,a,J,Ehat = 1.0,0.9,2.5,1.4
Rs=PolynomialRing(QQ,['s','Y']); s,Y=Rs.gens()
lam=[QQ(1200),QQ(-2300),QQ(-11428),QQ(-5519),QQ(24700),QQ(62500)]
qpoly=sum(lam[i]*s**i for i in range(6)); X=SageRS(Y**2-qpoly,prec=80)
omega=np.array(X.matrix_of_integral_values([Rs(1),s]),dtype=complex)[:,:2]
tau=np.array(X.riemann_matrix(),dtype=complex); ominv=np.linalg.inv(omega)
def theta_d(z,av,bv,order=0):
    av=np.array(av);bv=np.array(bv);zz=z+tau@av+bv
    pref=np.exp(1j*np.pi*(av@tau@av)+2j*np.pi*(av@(z+bv)));th=complex(RiemannTheta(zz,tau))
    if order==0: return pref*th
    g=np.array([complex(RiemannTheta(zz,tau,derivs=[e])) for e in ([1,0],[0,1])])
    if order==1: return 2j*np.pi*av+g/th
    H=np.array([[complex(RiemannTheta(zz,tau,derivs=[p,q])) for q in ([1,0],[0,1])] for p in([1,0],[0,1])])
    return H/th-np.outer(g/th,g/th)
halfs=[np.array(v)/2 for v in itertools.product([0,1],repeat=2)]
odd=[(av,bv) for av in halfs for bv in halfs if abs(theta_d(np.zeros(2),av,bv))<1e-6]
av,bv=odd[1]
def qn(sv): return sum(float(lam[i])*sv**i for i in range(6))
rmin=4.046197656444178; s_b=1.0/rmin
def Iu(sv):
    U=np.sqrt(s_b-sv);g0=lambda u:2*u/np.sqrt(qn(s_b-u**2));g1=lambda u:2*u*(s_b-u**2)/np.sqrt(qn(s_b-u**2))
    return -np.array([quad(g0,0,U,limit=200)[0],quad(g1,0,U,limit=200)[0]])
def w_of(sv): return ominv@Iu(sv)
e_plus=w_of(0.0); e_minus=-e_plus
# g_i simbolici (parti principali)
Es=sp.Symbol('E'); E0=sp.Rational(7,5)
q6=sp.Symbol('q6'); # ricostruisco g_i:
r_,ss=sp.symbols('r s'); Ms,as_,Js=sp.Rational(1),sp.Rational(9,10),sp.Rational(5,2)
Dl=r_**2-2*Ms*r_+as_**2;Em=(Es**2-1)*r_+2*Ms
Ssym=sp.expand(r_*(r_-2*Ms)*Em*(r_*Dl-Js**2*Em)); q6e=sp.expand(ss**6*Ssym.subs(r_,1/ss))
gser=sp.series(1/sp.sqrt(q6e),ss,0,4).removeO(); gsym=[sp.simplify(gser.coeff(ss,i)) for i in range(3)]
gv=[float(gsym[i].subs(Es,E0)) for i in range(3)]
print("g_i (num) =",np.round(gv,5))
# named primitives
def Om(rv):   # 3a specie log theta-ratio (canonica: residuo 1 in s a e_+)
    w=w_of(1.0/rv); return np.log(theta_d(w-e_plus,av,bv,0)/theta_d(w-e_minus,av,bv,0))
def Zvec(rv): # zeta_delta agli e_pm (2a specie), 2-vettori
    w=w_of(1.0/rv); return theta_d(w-e_plus,av,bv,1),theta_d(w-e_minus,av,bv,1)
def Pvec(rv): # wp_delta agli e_pm (2a specie triplo), matrici 2x2
    w=w_of(1.0/rv); return theta_d(w-e_plus,av,bv,2),theta_d(w-e_minus,av,bv,2)
def uhol(rv):
    s_to=1.0/rv
    I0=quad(lambda t:1.0/np.sqrt(qn(t)),1.0/60.0,s_to,limit=300)[0]
    I1=quad(lambda t:t/np.sqrt(qn(t)),1.0/60.0,s_to,limit=300)[0]
    return np.array([I0,I1])
def Sn(x):
    E2=1.96; Dl=x**2-2*x+0.81; Emm=(E2-1)*x+2
    return x*(x-2)*Emm*(x*Dl-6.25*Emm)   # M=1,a=0.9,E=7/5,J=5/2
def Uk(k,rv): return quad(lambda x:x**k/np.sqrt(Sn(x)),60.0,rv,limit=300)[0]
# --- normalizzazioni canoniche di Om,Z,P via U_2 (puro 3a specie) su range LARGO ---
# U_2 = A*Om + hol(u1,u2)+c.  Om ben separato (log-polo) da u1,u2 -> ben condizionato.
rg=np.array([50.,30.,20.,12.,8.,6.,5.,4.5])
def build_and_solve(k, cols):
    Bcols=[];
    for rv in rg: Bcols.append([c(rv) for c in cols])
    B=np.array(Bcols,dtype=complex); y=np.array([Uk(k,rv) for rv in rg],dtype=complex)
    coef,_,_,_=np.linalg.lstsq(B,y,rcond=None)
    cond=np.linalg.cond(B); res=np.max(np.abs(B@coef-y))
    return coef,cond,res
# U_2: base {Om, u1, u2, 1}
c2=lambda rv:Om(rv); u1=lambda rv:uhol(rv)[0]; u2=lambda rv:uhol(rv)[1]; one=lambda rv:1.0
coef2,cond2,res2=build_and_solve(2,[c2,u1,u2,one])
print(f"\nU_2 = A*Om + a*u1 + b*u2 + c : cond={cond2:.1e} res={res2:.1e}")
print(f"  A(Om)={coef2[0]:+.5f}  (atteso -g0={-gv[0]:+.5f})  a(u1)={coef2[1]:+.5f} b(u2)={coef2[2]:+.5f} c={coef2[3]:+.5f}")
print("=> se A ~ -g0 e cond piccola: parte OLOMORFA (a,b) determinata BEN CONDIZIONATA (no fit degenere).")
print("   a,b = coeff olomorfi: numerici (dipendono dai periodi) ma da principio, come C0/Ce separatrice.")
