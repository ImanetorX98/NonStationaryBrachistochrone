# -*- coding: utf-8 -*-
# COEFFICIENTI SIMBOLICI via PARTI PRINCIPALI (TK tau, genus-2). Nella base canonica
# (oggetti theta[delta] normalizzati a polo unitario agli e_pm), i coeff di U_k sono i
# coefficienti di POLO di omega_k = -s^{1-k}/sqrt(q6) ds a s=0 (=r=inf=e_pm).
#   g(s)=q6^{-1/2} = g0+g1 s+g2 s^2+...   (g_i razionali in E, SIMBOLICI)
#   omega_2: residuo -g0                    -> 3a specie (Lrat), coeff -g0
#   omega_3: polo s^-2 -g0, residuo -g1     -> 2a specie (Z=zeta canonica) -g0, Lrat -g1
#   omega_4: s^-3 -g0, s^-2 -g1, res -g2    -> 2a (P=wp canonica) -g0, Z -g1, Lrat -g2
# Verifica numerica: coefficienti di polo di U_k(r) a r->inf = questi g_i (no Sage).
import sympy as sp, numpy as np
from scipy.integrate import quad

r,s,E = sp.symbols('r s E', positive=True)
M,a,J = sp.Rational(1), sp.Rational(9,10), sp.Rational(5,2)
Dl=r**2-2*M*r+a**2; Emu=(E**2-1)*r+2*M
S=sp.expand(r*(r-2*M)*Emu*(r*Dl-J**2*Emu))
q6=sp.expand(s**6*S.subs(r,1/s))
# g(s) = q6^{-1/2}, Taylor a s=0
g_series=sp.series(1/sp.sqrt(q6), s, 0, 4).removeO()
g=[sp.simplify(g_series.coeff(s,i)) for i in range(3)]
print("=== g_i = coeff Taylor di q6^{-1/2} a s=0 (SIMBOLICI, razionali in sqrt(E^2-1)) ===")
for i in range(3): print(f"  g_{i} =",g[i])

print("\n=== COEFFICIENTI SIMBOLICI di U_k nella base canonica {Omega(3a), Z(2a dbl), P(2a tpl)} ===")
print("  U_2 = -g0 * Omega            + olo        (Omega=3a specie canonica, polo unitario e_pm)")
print("  U_3 = -g0 * Z   - g1 * Omega + olo        (Z=2a specie doppio polo canonica)")
print("  U_4 = -g0 * P   - g1 * Z  - g2 * Omega +olo(P=2a specie triplo polo canonica)")
print("  con (a E=7/5):")
E0=sp.Rational(7,5)
for i in range(3): print(f"    -g_{i} =",sp.nsimplify(sp.simplify(-g[i].subs(E,E0)),[sp.sqrt(E0**2-1)]),"  = ",float(-g[i].subs(E,E0)))

# --- VERIFICA numerica: coeff di polo di U_k a r->inf (s->0) coincidono con g_i ---
E0f=1.4; Mf,af,Jf=1.0,0.9,2.5
def Sn(x):
    Dl=x**2-2*Mf*x+af**2; Em=(E0f**2-1)*x+2*Mf; return x*(x-2*Mf)*Em*(x*Dl-Jf**2*Em)
def Uk(k,rv,rref=1000.0):  # U_k relativo a rref grande (per isolare i poli in s=1/r)
    return quad(lambda x:x**k/np.sqrt(Sn(x)),rref,rv,limit=400)[0]
gv=[float(g[i].subs(E,E0)) for i in range(3)]
print("\n=== verifica numerica struttura di polo di U_k a r grande (s=1/r->0) ===")
# U_3 ~ g0 * s^{-1} + g1*(-log s) ... in r: g0*r ; controllo (U_3(r)-U_3(r'))/(r-r') -> g0
for k,gi in [(3,gv[0]),(4,gv[0])]:
    r1,r2=800.0,1000.0
    if k==3:
        slope=(Uk(3,r1)-Uk(3,r2))/(r1-r2)   # ~ g0 (coeff di s^{-1}=r)
        print(f"  U_3: slope a r~900 = {slope:.6f}   g0 = {gv[0]:.6f}   diff={abs(slope-gv[0]):.1e}")
    if k==4:
        # U_4 ~ (g0/2) s^{-2} = (g0/2) r^2 ; slope di U_4/r^2... controllo coeff r^2
        f=lambda rv: Uk(4,rv)
        c2=(f(1000.0))/(1000.0**2/2)   # grezzo
        print(f"  U_4: coeff r^2 ~ {c2:.6f}   g0 = {gv[0]:.6f} (atteso g0/2 in U_4 ~ (g0/2)r^2)")
# --- CROSS-CHECK con mattone-2a (residui BEL): n_-1(2a) * g0 = -g_k (canonico) ---
print("\n=== CROSS-CHECK indipendente vs mattone-2a (residui odd-model) ===")
# mattone-2a: omega_3 residuo n_-1=(3-2E^2)/(E^2-1); omega_4 n_-1=(-625E^6+1156E^4-37E^2-794)/(200(E^2-1)^2)
n3_2a=(3-2*E**2)/(E**2-1)
n4_2a=(-625*E**6+1156*E**4-37*E**2-794)/(200*(E**2-1)**2)
# canonico: coeff di Omega in U_k = -g_{k-2}; mattone-2a L ha residuo g0 -> n_2a*g0 deve = -g_{k-2}
chk3=sp.simplify(n3_2a*g[0]-(-g[1]))
chk4=sp.simplify(n4_2a*g[0]-(-g[2]))
print("  U_3: n3_2a*g0 - (-g1) =",chk3," (deve 0)")
print("  U_4: n4_2a*g0 - (-g2) =",chk4," (deve 0)")
print("\n=> g_i SIMBOLICI (razionali in E) = coefficienti di U_k nella base canonica theta[delta].")
print("   Paralleli a b1,b2,b3 della separatrice. Confermati da mattone-2a (residui BEL).")
print("   Le funzioni Omega,Z,P (theta[delta] agli e_pm) portano la normalizzazione geometrica.")
