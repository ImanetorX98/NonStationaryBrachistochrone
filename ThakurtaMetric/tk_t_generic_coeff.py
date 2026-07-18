# -*- coding: utf-8 -*-
# BRACHISTOCRONA GENERICA (J qualsiasi, genus-2) Thakurta-Kerr ramo t. Coeff SIMBOLICI in (a,E,J), M=1.
# Curva R6=r Q2 DE ; sorgente dE F_t (N_t=R6 dE K_t - 1/2 K_t dE R6) ; clock rho_t=P3+R_Delta/Delta.
# Mattoni: c_k^E (riduzione 2a specie), Q_kj (b=P3 coeff), g_i (Taylor R6 dispari), P3 (clock poly).
import sympy as sp, numpy as np
r,s,a,E,J=sp.symbols('r s a E J',positive=True)
M=1
DE=(E**2-1)*r+2*M; Delta=r**2-2*M*r+a**2
Q2=(2*E**2*J**2*r-E**2*J**2*r**2-4*E**2*J*a*r+2*E**2*a**2*r+E**2*a**2*r**2+E**2*r**4
    +4*J**2-4*J**2*r+J**2*r**2-8*J*a+4*J*a*r+4*a**2)   # M=1
R6=sp.expand(r*Q2*DE)
Kt=r*DE*(J*(r-2*M)+2*M*a)/Delta

# --- MATTONE 1: dE F_t = N_t/R6^{3/2}, N_t = R6 dE K_t - 1/2 K_t dE R6 ---
N=sp.together(R6*sp.diff(Kt,E)-sp.Rational(1,2)*Kt*sp.diff(R6,E))
# N ha Delta al denominatore: la parte 2a specie e' il numeratore/Delta ridotto. Semplifico:
N=sp.simplify(N)
print("MATTONE 1: N_t calcolato (denominatore):",sp.denom(sp.cancel(N)))

# --- clock: rho_t = P3 + R_Delta/Delta ; P3 = quoziente ---
rho_t=sp.cancel(sp.together((E**2*r**3-2*M*a*Kt/r)/((r-2*M)/r)))
P3poly,_=sp.div(sp.Poly(sp.numer(rho_t),r),sp.Poly(sp.denom(rho_t),r))
P3=sp.expand(P3poly.as_expr())
print("\nCLOCK ramo t: P3(r) (parte 2a specie, b-vector) =",P3)
bvec=[P3.coeff(r,i) for i in range(sp.degree(P3,r)+1)]
print("  b-vector (coeff P3):",[sp.simplify(bi) for bi in bvec])

# --- MATTONE 2: riduzione 2a specie di dE F_t.  N_t razionale (Delta) -> uso R6 come curva ---
# dE F_t = d(A5/sqrtR6) + sum c_k r^k/sqrtR6.  N_num = numeratore di N su 1/R6^{3/2}:
# scrivo N_t effettivo (num tale che dE F_t=N_t/R6^{3/2}): moltiplico per Delta se serve.
NR=sp.cancel(N); Nnum=sp.numer(NR); Nden=sp.denom(NR)   # dE F_t = Nnum/(Nden R6^{... })
# dE F_t = N/R6^{3/2} con N=NR ; ma NR ha Delta al denom -> la 2a specie e' su curva R6 con
# eventuale polo orizzonte. Riduco Nnum contro R6*Nden:
print("\nMATTONE 2: riduzione 2a specie (numeratore su R6). deg Nnum =",sp.degree(sp.Poly(sp.expand(Nnum),r)))
# tento riduzione: dE F_t * R6^{3/2} = NR ; scrivo NR = d(A5/sqrtR6)*R6^{3/2}/... approssimo
# uso la forma: 2 Nnum = 2 R6 A5' - A5 R6' + 2 R6 (sum c_k) , con Nden assorbito -> se Nden=Delta,
# moltiplico tutto per Delta. Verifico numericamente invece (a=9/10).
sub={a:sp.Rational(9,10),E:sp.Rational(6,5),J:sp.Rational(5,2)}
dEF=sp.diff(Kt/sp.sqrt(R6),E)
Neff=sp.expand(sp.cancel(sp.simplify(dEF*R6**sp.Rational(3,2))))   # N_t (polinomio)
print("  N_t denom:",sp.denom(sp.cancel(dEF*R6**sp.Rational(3,2)))," ; N_t poly deg:",sp.degree(sp.Poly(Neff,r)))
# riduzione 2a specie: 2 N_t = 2 R6 A5' - A5 R6' + 2 R6 sum c_k r^k
Ac=[sp.Symbol(f'A{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
A5=sum(Ac[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
eq=sp.expand(2*Neff-(2*R6*sp.diff(A5,r)-A5*sp.diff(R6,r)+2*R6*Mp))
sol=sp.solve(sp.Poly(eq,r).all_coeffs(),Ac+ck,dict=True)
assert sol, "riduzione non chiude"
sol=sol[0]; cvec=[sol[ck[i]] for i in range(5)]
print("\nMATTONE 2b: c_k^E SIMBOLICI (razionali in a,E,J) - chiude con A5 deg5. Esempi (grado):")
for i in range(5):
    num,den=sp.fraction(sp.cancel(cvec[i]))
    print(f"  c{i}: razionale, deg_num(E)={sp.degree(sp.Poly(num,E))} deg_den(E)={sp.degree(sp.Poly(den,E))}")
bb=[bi for bi in bvec]+[0,0]
Qnz=[(k,j) for k in range(5) for j in range(k+1,5) if sp.simplify(cvec[k]*bb[j]-cvec[j]*bb[k])!=0]
print("  Q_kj = c_k b_j - c_j b_k (b=P3 clock) NON nulli:",Qnz)
# VERIFICA numerica riduzione completa dE F = d(A5/sqrtR6)+sum c_k r^k/sqrtR6
A5e=A5.subs(sol); rhs=sp.diff(A5e/sp.sqrt(R6),r)+sum(cvec[k]*r**k for k in range(5))/sp.sqrt(R6)
rhsn=sp.lambdify(r,rhs.subs(sub),'numpy'); dEFn=sp.lambdify(r,dEF.subs(sub),'numpy')
print("  VERIFICA riduzione (a=9/10,E=6/5,J=5/2):")
for rv in [20.,15.,10.]:
    print(f"    r={rv}: dE F={dEFn(rv):+.6e}  d(A5)+sum c_k={float(rhsn(rv)):+.6e}  diff={abs(dEFn(rv)-float(rhsn(rv))):.1e}")

# --- MATTONE 3: g_i = Taylor di q6^{-1/2}, q6=s^6 R6(1/s) ---
q6=sp.expand(s**6*R6.subs(r,1/s))
gser=sp.series(1/sp.sqrt(q6),s,0,3).removeO()
g=[sp.simplify(gser.coeff(s,i)) for i in range(3)]
print("\nMATTONE 3: g_i (parti principali di R6) simbolici in (a,E,J):")
for i in range(3): print(f"  g{i} =",g[i])

# ===== VERIFICA numerica riduzione (a=9/10,E=6/5,J=5/2) =====
print("\n=== VERIFICA: struttura sorgente/curva (a=9/10,E=6/5,J=5/2) ===")
R6n=sp.lambdify(r,R6.subs(sub),'numpy'); dEFn=sp.lambdify(r,dEF.subs(sub),'numpy')
Neffn=sp.lambdify(r,Neff.subs(sub),'numpy')
for rv in [20.,15.,10.]:
    lhs=dEFn(rv); rhs=Neffn(rv)/R6n(rv)**1.5
    print(f"  r={rv}: dE F={lhs:+.6e}  N_t/R6^1.5={rhs:+.6e}  diff={abs(lhs-rhs):.1e}")
print("\n=> TK-t generico: g_i simbolici in (a,E,J); P3 (clock) simbolico; N_t sorgente (con Delta")
print("   orizzonte). c_k richiedono la riduzione con polo orizzonte (come il clock t separatrice).")
