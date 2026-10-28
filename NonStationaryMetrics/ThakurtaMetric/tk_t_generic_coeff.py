# -*- coding: utf-8 -*-
# GENERIC BRACHISTOCHRONE (any J, genus 2), Thakurta-Kerr t branch. SYMBOLIC coefficients in (a,E,J), M=1.
# Curve R6=r Q2 DE ; source dE F_t (N_t=R6 dE K_t - 1/2 K_t dE R6) ; clock rho_t=P3+R_Delta/Delta.
# Building blocks: c_k^E (second-kind reduction), Q_kj (b = coefficients of P3), g_i (Taylor expansion of the odd R6), P3 (clock polynomial).
import sympy as sp, numpy as np
r,s,a,E,J=sp.symbols('r s a E J',positive=True)
M=1
DE=(E**2-1)*r+2*M; Delta=r**2-2*M*r+a**2
Q2=(2*E**2*J**2*r-E**2*J**2*r**2-4*E**2*J*a*r+2*E**2*a**2*r+E**2*a**2*r**2+E**2*r**4
    +4*J**2-4*J**2*r+J**2*r**2-8*J*a+4*J*a*r+4*a**2)   # M=1
R6=sp.expand(r*Q2*DE)
Kt=r*DE*(J*(r-2*M)+2*M*a)/Delta

# --- BUILDING BLOCK 1: dE F_t = N_t/R6^{3/2}, N_t = R6 dE K_t - 1/2 K_t dE R6 ---
N=sp.together(R6*sp.diff(Kt,E)-sp.Rational(1,2)*Kt*sp.diff(R6,E))
# N has Delta in the denominator: the second-kind part is the reduced numerator/Delta. Simplify:
N=sp.simplify(N)
print("BUILDING BLOCK 1: N_t computed (denominator):",sp.denom(sp.cancel(N)))

# --- clock: rho_t = P3 + R_Delta/Delta ; P3 = quotient ---
rho_t=sp.cancel(sp.together((E**2*r**3-2*M*a*Kt/r)/((r-2*M)/r)))
P3poly,_=sp.div(sp.Poly(sp.numer(rho_t),r),sp.Poly(sp.denom(rho_t),r))
P3=sp.expand(P3poly.as_expr())
print("\nt-branch CLOCK: P3(r) (second-kind part, b vector) =",P3)
bvec=[P3.coeff(r,i) for i in range(sp.degree(P3,r)+1)]
print("  b vector (coefficients of P3):",[sp.simplify(bi) for bi in bvec])

# --- BUILDING BLOCK 2: second-kind reduction of dE F_t.  N_t rational (Delta) -> R6 used as the curve ---
# dE F_t = d(A5/sqrtR6) + sum c_k r^k/sqrtR6.  N_num = numerator of N over 1/R6^{3/2}:
# write the effective N_t (numerator such that dE F_t=N_t/R6^{3/2}): multiply by Delta if needed.
NR=sp.cancel(N); Nnum=sp.numer(NR); Nden=sp.denom(NR)   # dE F_t = Nnum/(Nden R6^{... })
# dE F_t = N/R6^{3/2} with N=NR ; but NR has Delta in the denominator -> the second kind lives on the curve R6 with
# a possible horizon pole. Reduce Nnum against R6*Nden:
print("\nBUILDING BLOCK 2: second-kind reduction (numerator over R6). deg Nnum =",sp.degree(sp.Poly(sp.expand(Nnum),r)))
# attempt the reduction: dE F_t * R6^{3/2} = NR ; write NR = d(A5/sqrtR6)*R6^{3/2}/... approximately
# use the form: 2 Nnum = 2 R6 A5' - A5 R6' + 2 R6 (sum c_k) , with Nden absorbed -> if Nden=Delta,
# multiply everything by Delta. Checked numerically instead (a=9/10).
sub={a:sp.Rational(9,10),E:sp.Rational(6,5),J:sp.Rational(5,2)}
dEF=sp.diff(Kt/sp.sqrt(R6),E)
Neff=sp.expand(sp.cancel(sp.simplify(dEF*R6**sp.Rational(3,2))))   # N_t (polynomial)
print("  N_t denom:",sp.denom(sp.cancel(dEF*R6**sp.Rational(3,2)))," ; N_t poly deg:",sp.degree(sp.Poly(Neff,r)))
# second-kind reduction: 2 N_t = 2 R6 A5' - A5 R6' + 2 R6 sum c_k r^k
Ac=[sp.Symbol(f'A{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
A5=sum(Ac[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
eq=sp.expand(2*Neff-(2*R6*sp.diff(A5,r)-A5*sp.diff(R6,r)+2*R6*Mp))
sol=sp.solve(sp.Poly(eq,r).all_coeffs(),Ac+ck,dict=True)
assert sol, "the reduction does not close"
sol=sol[0]; cvec=[sol[ck[i]] for i in range(5)]
print("\nBUILDING BLOCK 2b: SYMBOLIC c_k^E (rational in a,E,J) - it closes with A5 of degree 5. Examples (degree):")
for i in range(5):
    num,den=sp.fraction(sp.cancel(cvec[i]))
    print(f"  c{i}: rational, deg_num(E)={sp.degree(sp.Poly(num,E))} deg_den(E)={sp.degree(sp.Poly(den,E))}")
bb=[bi for bi in bvec]+[0,0]
Qnz=[(k,j) for k in range(5) for j in range(k+1,5) if sp.simplify(cvec[k]*bb[j]-cvec[j]*bb[k])!=0]
print("  Q_kj = c_k b_j - c_j b_k (b = P3 clock), NON-zero ones:",Qnz)
# NUMERICAL CHECK of the complete reduction dE F = d(A5/sqrtR6)+sum c_k r^k/sqrtR6
A5e=A5.subs(sol); rhs=sp.diff(A5e/sp.sqrt(R6),r)+sum(cvec[k]*r**k for k in range(5))/sp.sqrt(R6)
rhsn=sp.lambdify(r,rhs.subs(sub),'numpy'); dEFn=sp.lambdify(r,dEF.subs(sub),'numpy')
print("  CHECK of the reduction (a=9/10,E=6/5,J=5/2):")
for rv in [20.,15.,10.]:
    print(f"    r={rv}: dE F={dEFn(rv):+.6e}  d(A5)+sum c_k={float(rhsn(rv)):+.6e}  diff={abs(dEFn(rv)-float(rhsn(rv))):.1e}")

# --- BUILDING BLOCK 3: g_i = Taylor expansion of q6^{-1/2}, q6=s^6 R6(1/s) ---
q6=sp.expand(s**6*R6.subs(r,1/s))
gser=sp.series(1/sp.sqrt(q6),s,0,3).removeO()
g=[sp.simplify(gser.coeff(s,i)) for i in range(3)]
print("\nBUILDING BLOCK 3: g_i (principal parts of R6), symbolic in (a,E,J):")
for i in range(3): print(f"  g{i} =",g[i])

# ===== NUMERICAL CHECK of the reduction (a=9/10,E=6/5,J=5/2) =====
print("\n=== CHECK: structure of the source/curve (a=9/10,E=6/5,J=5/2) ===")
R6n=sp.lambdify(r,R6.subs(sub),'numpy'); dEFn=sp.lambdify(r,dEF.subs(sub),'numpy')
Neffn=sp.lambdify(r,Neff.subs(sub),'numpy')
for rv in [20.,15.,10.]:
    lhs=dEFn(rv); rhs=Neffn(rv)/R6n(rv)**1.5
    print(f"  r={rv}: dE F={lhs:+.6e}  N_t/R6^1.5={rhs:+.6e}  diff={abs(lhs-rhs):.1e}")
print("\n=> generic TK-t: g_i symbolic in (a,E,J); P3 (clock) symbolic; N_t source (with the Delta")
print("   horizon). The c_k require the reduction with a horizon pole (as for the t clock on the separatrix).")
