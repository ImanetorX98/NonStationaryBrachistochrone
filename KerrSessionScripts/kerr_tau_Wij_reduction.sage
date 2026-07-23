# -*- coding: utf-8 -*-
# CHIUSURA W_ij (TK tau, J generico off-separatrice) -- MATTONE 1 (sympy, no Sage).
# Riduzione di 2a specie degli integrali abeliani U_k=int r^k dr/sqrt(S) alla base
# minima di coomologia {U_0,U_1 (1a specie), U_2,U_3 (2a/3a specie)} + forma esatta,
# con COEFFICIENTI SIMBOLICI (razionali in E), verificata numericamente.
# Identifica anche il residuo 3a specie a r=infinito (analogo rho_0 della separatrice:
# sorgente del dilog genus-2).
import sympy as sp, numpy as np
from scipy.integrate import quad

r,E = sp.symbols('r E', positive=True)
M,a,J = sp.Rational(1), sp.Rational(9,10), sp.Rational(5,2)   # off-separatrice (J_c^tau=a/E=0.643)
Dl = r**2-2*M*r+a**2; Emu=(E**2-1)*r+2*M
S = sp.expand(r*(r-2*M)*Emu*(r*Dl-J**2*Emu))                  # sestica, 6 radici semplici -> genus 2
Sp = sp.diff(S,r)
print("S deg =", sp.degree(S,r), " (genus 2 se 6 radici semplici)")

# --- riduzione: r^k = (2 P' S - P S')/2 + resto,  con d(P/sqrtS)=(2P'S-P S')/(2 S^{3/2}) ---
# => int r^k/sqrtS = P/sqrtS + int resto/sqrtS.  Base irriducibile: r^0..r^3.
def reduce_pow(k):
    # cerca P deg (k-1) e resto = sum_{j=0}^{3} m_j r^j  tale che
    #   r^k * S = (2 P' S - P S')/2  ... NO: r^k = (2P'S - P S')/2 /S? sistemiamo su r^k * ...
    # Vogliamo:  r^k = d/dr[P] applicato... usiamo identita' su POLINOMI:
    #   r^k = (2 P'(r) S(r) - P(r) S'(r))/2 / S(r) ?  -> non polinomiale.
    # Corretto: r^k/sqrtS - d(P/sqrtS)/dr = resto/sqrtS  con
    #   d(P/sqrtS)/dr = (2 P' S - P S')/(2 S^{3/2})
    # => r^k/sqrtS = (2P'S - P S')/(2 S^{3/2}) + resto/sqrtS
    # moltiplico per S^{3/2}:  r^k S = (2P'S - P S')/2 + resto*S
    # cioe' 2 r^k S = 2P'S - P S' + 2 resto S. Poly identity in r.
    dP = k+2  # deg P: bilancio gradi (P S' deg = degP+5; r^k S deg = k+6) -> degP=k+1
    dP = k+1
    Pc=[sp.Symbol(f'p{i}') for i in range(dP+1)]
    P=sum(Pc[i]*r**i for i in range(dP+1))
    mc=[sp.Symbol(f'm{i}') for i in range(4)]
    resto=sum(mc[i]*r**i for i in range(4))
    eq=sp.expand(2*r**k*S - (2*sp.diff(P,r)*S - P*Sp + 2*resto*S))
    sol=sp.solve(sp.Poly(eq,r).all_coeffs(), Pc+mc, dict=True)
    if not sol:
        return None, None
    sol=sol[0]
    return sp.expand(P.subs(sol)), [sp.simplify(resto.subs(sol).coeff(r,i)) for i in range(4)]

print("\n=== Riduzione r^k/sqrtS = d(P_k/sqrtS) + sum_j m_kj r^j/sqrtS (base j=0..3) ===")
red={}
for k in [2,3,4]:
    P,m = reduce_pow(k)
    if m is None:
        print(f"k={k}: IRRIDUCIBILE (nessun P deg<={k+1}) -> generatore di coomologia indipendente")
    else:
        red[k]=(P,m); print(f"k={k}: m_j =", [sp.nsimplify(sp.simplify(mj)) for mj in m], " (banale: e' esso stesso base)")

# --- VERIFICA numerica dell'identita' di riduzione (E=7/5) ---
E0=sp.Rational(7,5); Mf,af,Ef,Jf=1.0,0.9,1.4,2.5
def Sn(x):
    Dl=x**2-2*Mf*x+af**2; Em=(Ef**2-1)*x+2*Mf
    return x*(x-2*Mf)*Em*(x*Dl-Jf**2*Em)
def sq(x): return np.sqrt(Sn(x))
print("\n=== VERIFICA numerica riduzione (integrando, deve ~0) ===")
for k in list(red):
    P,m=red[k]
    Pn=sp.lambdify(r,P.subs(E,E0),'numpy'); Spn=sp.lambdify(r,Sp.subs(E,E0),'numpy'); Sn_l=sp.lambdify(r,S.subs(E,E0),'numpy')
    mn=[float(mj.subs(E,E0)) for mj in m]
    def lhs(x): return x**k/sq(x)
    def dPovsqrtS(x):  # d/dx (P/sqrtS)
        return (2*Pn(x)*0)  # placeholder
    # d(P/sqrtS)/dx = (2 P' S - P S')/(2 S^{3/2}); calcolo P' numerico via sympy
    Ppn=sp.lambdify(r,sp.diff(P,r).subs(E,E0),'numpy')
    def rhs(x):
        Sx=Sn_l(x); return (2*Ppn(x)*Sx - Pn(x)*Spn(x))/(2*Sx**1.5) + sum(mn[j]*x**j/sq(x) for j in range(4))
    errs=[abs(lhs(x)-rhs(x)) for x in [11.0,9.0,7.0,5.0]]
    print(f"  k={k}: max|lhs-rhs| = {max(errs):.2e}")

# --- residuo 3a specie a r=infinito dei generatori U_2,U_3 (analogo rho_0) ---
# r^k dr/sqrtS ~ r^k * dr / (sqrt(a6) r^3) = r^{k-3}/sqrt(a6) dr ; residuo (coeff 1/r) a inf:
a6=sp.LC(sp.Poly(S,r))
print("\n=== struttura a r=infinito (a6 = coeff r^6) ===")
print("  a6 =", sp.simplify(a6))
for k in [0,1,2,3,4]:
    # r^k/sqrtS ~ r^{k-3}/sqrt(a6) (1 + O(1/r)); residuo=coeff di r^{-1}=> k-3=-1 => k=2 (leading)
    print(f"  k={k}: leading r^{{{k-3}}}/sqrt(a6) -> {'RESIDUO 3a specie (k=2 leading)' if k==2 else '2a specie/olomorfa'}")
print("\n  => U_2 e' il generatore di 3a specie (residuo ~1/sqrt(a6) ai due punti r=inf);")
print("     U_0,U_1 olomorfe; U_3 (e U_4 ridotto) 2a specie. Coeff m_kj SIMBOLICI sopra.")
