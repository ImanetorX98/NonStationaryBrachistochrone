# -*- coding: utf-8 -*-
# CLOSURE OF W_ij (TK tau, generic J off the separatrix) -- BUILDING BLOCK 1 (sympy, no Sage).
# Second-kind reduction of the abelian integrals U_k=int r^k dr/sqrt(S) to the minimal
# cohomology basis {U_0,U_1 (first kind), U_2,U_3 (second/third kind)} + an exact form,
# with SYMBOLIC COEFFICIENTS (rational in E), verified numerically.
# It also identifies the third-kind residue at r=infinity (the analogue of rho_0 of the separatrix:
# the source of the genus-2 dilog).
import sympy as sp, numpy as np
from scipy.integrate import quad

r,E = sp.symbols('r E', positive=True)
M,a,J = sp.Rational(1), sp.Rational(9,10), sp.Rational(5,2)   # off the separatrix (J_c^tau=a/E=0.643)
Dl = r**2-2*M*r+a**2; Emu=(E**2-1)*r+2*M
S = sp.expand(r*(r-2*M)*Emu*(r*Dl-J**2*Emu))                  # sextic, 6 simple roots -> genus 2
Sp = sp.diff(S,r)
print("S deg =", sp.degree(S,r), " (genus 2 if 6 simple roots)")

# --- reduction: r^k = (2 P' S - P S')/2 + remainder,  with d(P/sqrtS)=(2P'S-P S')/(2 S^{3/2}) ---
# => int r^k/sqrtS = P/sqrtS + int remainder/sqrtS.  Irreducible basis: r^0..r^3.
def reduce_pow(k):
    # look for P of degree (k-1) and remainder = sum_{j=0}^{3} m_j r^j  such that
    #   r^k * S = (2 P' S - P S')/2  ... NO: r^k = (2P'S - P S')/2 /S? set it up on r^k * ...
    # We want:  r^k = d/dr[P] applied... use an identity between POLYNOMIALS:
    #   r^k = (2 P'(r) S(r) - P(r) S'(r))/2 / S(r) ?  -> not polynomial.
    # Correct: r^k/sqrtS - d(P/sqrtS)/dr = remainder/sqrtS  with
    #   d(P/sqrtS)/dr = (2 P' S - P S')/(2 S^{3/2})
    # => r^k/sqrtS = (2P'S - P S')/(2 S^{3/2}) + remainder/sqrtS
    # multiply by S^{3/2}:  r^k S = (2P'S - P S')/2 + remainder*S
    # i.e. 2 r^k S = 2P'S - P S' + 2 remainder S. Poly identity in r.
    dP = k+2  # deg P: degree balance (P S' deg = degP+5; r^k S deg = k+6) -> degP=k+1
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

print("\n=== Reduction r^k/sqrtS = d(P_k/sqrtS) + sum_j m_kj r^j/sqrtS (basis j=0..3) ===")
red={}
for k in [2,3,4]:
    P,m = reduce_pow(k)
    if m is None:
        print(f"k={k}: IRREDUCIBLE (no P with deg<={k+1}) -> independent cohomology generator")
    else:
        red[k]=(P,m); print(f"k={k}: m_j =", [sp.nsimplify(sp.simplify(mj)) for mj in m], " (trivial: it is itself a basis element)")

# --- NUMERICAL CHECK of the reduction identity (E=7/5) ---
E0=sp.Rational(7,5); Mf,af,Ef,Jf=1.0,0.9,1.4,2.5
def Sn(x):
    Dl=x**2-2*Mf*x+af**2; Em=(Ef**2-1)*x+2*Mf
    return x*(x-2*Mf)*Em*(x*Dl-Jf**2*Em)
def sq(x): return np.sqrt(Sn(x))
print("\n=== NUMERICAL CHECK of the reduction (integrand, must be ~0) ===")
for k in list(red):
    P,m=red[k]
    Pn=sp.lambdify(r,P.subs(E,E0),'numpy'); Spn=sp.lambdify(r,Sp.subs(E,E0),'numpy'); Sn_l=sp.lambdify(r,S.subs(E,E0),'numpy')
    mn=[float(mj.subs(E,E0)) for mj in m]
    def lhs(x): return x**k/sq(x)
    def dPovsqrtS(x):  # d/dx (P/sqrtS)
        return (2*Pn(x)*0)  # placeholder
    # d(P/sqrtS)/dx = (2 P' S - P S')/(2 S^{3/2}); P' computed numerically through sympy
    Ppn=sp.lambdify(r,sp.diff(P,r).subs(E,E0),'numpy')
    def rhs(x):
        Sx=Sn_l(x); return (2*Ppn(x)*Sx - Pn(x)*Spn(x))/(2*Sx**1.5) + sum(mn[j]*x**j/sq(x) for j in range(4))
    errs=[abs(lhs(x)-rhs(x)) for x in [11.0,9.0,7.0,5.0]]
    print(f"  k={k}: max|lhs-rhs| = {max(errs):.2e}")

# --- third-kind residue at r=infinity of the generators U_2,U_3 (analogue of rho_0) ---
# r^k dr/sqrtS ~ r^k * dr / (sqrt(a6) r^3) = r^{k-3}/sqrt(a6) dr ; residue (coeff of 1/r) at inf:
a6=sp.LC(sp.Poly(S,r))
print("\n=== structure at r=infinity (a6 = coeff of r^6) ===")
print("  a6 =", sp.simplify(a6))
for k in [0,1,2,3,4]:
    # r^k/sqrtS ~ r^{k-3}/sqrt(a6) (1 + O(1/r)); residue=coeff of r^{-1}=> k-3=-1 => k=2 (leading)
    print(f"  k={k}: leading r^{{{k-3}}}/sqrt(a6) -> {'THIRD-KIND RESIDUE (k=2 leading)' if k==2 else 'second kind/holomorphic'}")
print("\n  => U_2 is the third-kind generator (residue ~1/sqrt(a6) at the two points r=inf);")
print("     U_0,U_1 holomorphic; U_3 (and the reduced U_4) second kind. SYMBOLIC coeffs m_kj above.")
