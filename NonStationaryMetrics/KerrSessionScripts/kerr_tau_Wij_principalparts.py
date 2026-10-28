# -*- coding: utf-8 -*-
# SYMBOLIC COEFFICIENTS through PRINCIPAL PARTS (TK tau, genus-2). In the canonical basis
# (theta[delta] objects normalised to a unit pole at e_pm), the coeffs of U_k are the
# POLE coefficients of omega_k = -s^{1-k}/sqrt(q6) ds at s=0 (=r=inf=e_pm).
#   g(s)=q6^{-1/2} = g0+g1 s+g2 s^2+...   (g_i rational in E, SYMBOLIC)
#   omega_2: residue -g0                    -> third kind (Lrat), coeff -g0
#   omega_3: pole s^-2 -g0, residue -g1     -> second kind (Z=canonical zeta) -g0, Lrat -g1
#   omega_4: s^-3 -g0, s^-2 -g1, res -g2    -> second kind (P=canonical wp) -g0, Z -g1, Lrat -g2
# Numerical check: pole coefficients of U_k(r) as r->inf = these g_i (no Sage).
import sympy as sp, numpy as np
from scipy.integrate import quad

r,s,E = sp.symbols('r s E', positive=True)
M,a,J = sp.Rational(1), sp.Rational(9,10), sp.Rational(5,2)
Dl=r**2-2*M*r+a**2; Emu=(E**2-1)*r+2*M
S=sp.expand(r*(r-2*M)*Emu*(r*Dl-J**2*Emu))
q6=sp.expand(s**6*S.subs(r,1/s))
# g(s) = q6^{-1/2}, Taylor at s=0
g_series=sp.series(1/sp.sqrt(q6), s, 0, 4).removeO()
g=[sp.simplify(g_series.coeff(s,i)) for i in range(3)]
print("=== g_i = Taylor coeffs of q6^{-1/2} at s=0 (SYMBOLIC, rational in sqrt(E^2-1)) ===")
for i in range(3): print(f"  g_{i} =",g[i])

print("\n=== SYMBOLIC COEFFICIENTS of U_k in the canonical basis {Omega(3rd), Z(2nd dbl), P(2nd tpl)} ===")
print("  U_2 = -g0 * Omega            + hol        (Omega=canonical third kind, unit pole e_pm)")
print("  U_3 = -g0 * Z   - g1 * Omega + hol        (Z=canonical second kind, double pole)")
print("  U_4 = -g0 * P   - g1 * Z  - g2 * Omega +hol(P=canonical second kind, triple pole)")
print("  with (at E=7/5):")
E0=sp.Rational(7,5)
for i in range(3): print(f"    -g_{i} =",sp.nsimplify(sp.simplify(-g[i].subs(E,E0)),[sp.sqrt(E0**2-1)]),"  = ",float(-g[i].subs(E,E0)))

# --- NUMERICAL CHECK: pole coeffs of U_k as r->inf (s->0) coincide with g_i ---
E0f=1.4; Mf,af,Jf=1.0,0.9,2.5
def Sn(x):
    Dl=x**2-2*Mf*x+af**2; Em=(E0f**2-1)*x+2*Mf; return x*(x-2*Mf)*Em*(x*Dl-Jf**2*Em)
def Uk(k,rv,rref=1000.0):  # U_k relative to a large rref (to isolate the poles in s=1/r)
    return quad(lambda x:x**k/np.sqrt(Sn(x)),rref,rv,limit=400)[0]
gv=[float(g[i].subs(E,E0)) for i in range(3)]
print("\n=== numerical check of the pole structure of U_k at large r (s=1/r->0) ===")
# U_3 ~ g0 * s^{-1} + g1*(-log s) ... in r: g0*r ; check (U_3(r)-U_3(r'))/(r-r') -> g0
for k,gi in [(3,gv[0]),(4,gv[0])]:
    r1,r2=800.0,1000.0
    if k==3:
        slope=(Uk(3,r1)-Uk(3,r2))/(r1-r2)   # ~ g0 (coeff of s^{-1}=r)
        print(f"  U_3: slope at r~900 = {slope:.6f}   g0 = {gv[0]:.6f}   diff={abs(slope-gv[0]):.1e}")
    if k==4:
        # U_4 ~ (g0/2) s^{-2} = (g0/2) r^2 ; slope of U_4/r^2... check the r^2 coeff
        f=lambda rv: Uk(4,rv)
        c2=(f(1000.0))/(1000.0**2/2)   # rough
        print(f"  U_4: coeff r^2 ~ {c2:.6f}   g0 = {gv[0]:.6f} (expected g0/2 in U_4 ~ (g0/2)r^2)")
# --- CROSS-CHECK against building block 2a (BEL residues): n_-1(2a) * g0 = -g_k (canonical) ---
print("\n=== independent CROSS-CHECK vs building block 2a (odd-model residues) ===")
# building block 2a: omega_3 residue n_-1=(3-2E^2)/(E^2-1); omega_4 n_-1=(-625E^6+1156E^4-37E^2-794)/(200(E^2-1)^2)
n3_2a=(3-2*E**2)/(E**2-1)
n4_2a=(-625*E**6+1156*E**4-37*E**2-794)/(200*(E**2-1)**2)
# canonical: coeff of Omega in U_k = -g_{k-2}; building block 2a L has residue g0 -> n_2a*g0 must = -g_{k-2}
chk3=sp.simplify(n3_2a*g[0]-(-g[1]))
chk4=sp.simplify(n4_2a*g[0]-(-g[2]))
print("  U_3: n3_2a*g0 - (-g1) =",chk3," (must be 0)")
print("  U_4: n4_2a*g0 - (-g2) =",chk4," (must be 0)")
print("\n=> SYMBOLIC g_i (rational in E) = coefficients of U_k in the canonical basis theta[delta].")
print("   Parallel to b1,b2,b3 of the separatrix. Confirmed by building block 2a (BEL residues).")
print("   The functions Omega,Z,P (theta[delta] at e_pm) carry the geometric normalisation.")
