# -*- coding: utf-8 -*-
# BRACHISTOCRONA GENERICA (J qualsiasi, genus-2) Vaidya tau (a=0, parametro m).
# Coefficienti SIMBOLICI (razionali in m,E,J), mattone per mattone come TK-tau:
#   c_k^m : riduzione 2a specie di dm F  ;  Q_kj = c_k b_j - c_j b_k (b=clock tau)
#   g_i   : parti principali (Taylor di q6^{-1/2})  ;  P(r),res(x) : T_alg elementare.
# Verifica numerica di ogni passo.
import sympy as sp, numpy as np
r,s,m,E,J=sp.symbols('r s m E J',positive=True)
# Vaidya (a=0): Delta=r(r-2m), K_tau=J r(r-2m)DE/Delta = J DE (Delta cancella)
DE=(E**2-1)*r+2*m
S=sp.expand(r*(r-2*m)*DE*(r**2*(r-2*m)-J**2*DE))   # sestica Vaidya tau
K=J*DE                                             # numeratore (Delta cancellato)
F=K/sp.sqrt(S)
b=[0,0,-2*m,1,0]                                   # clock tau: r^3-2m r^2

# --- MATTONE 1: dm F = N_m/S^{3/2}, N_m = S dm K - 1/2 K dm S ---
dmF=sp.diff(F,m)
N=sp.expand(S*sp.diff(K,m)-sp.Rational(1,2)*K*sp.diff(S,m))
chk1=sp.simplify(dmF-N/S**sp.Rational(3,2))
print("MATTONE 1: dm F - N/S^{3/2} =",chk1,"(deve 0)")

# --- MATTONE 2: riduzione 2a specie  2N = 2S A5' - A5 S' + 2 S sum c_k r^k ---
Ac=[sp.Symbol(f'A{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
A5=sum(Ac[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
eq=sp.expand(2*N-(2*S*sp.diff(A5,r)-A5*sp.diff(S,r)+2*S*Mp))
sol=sp.solve(sp.Poly(eq,r).all_coeffs(),Ac+ck,dict=True)[0]
cvec=[sp.simplify(sol[ck[i]]) for i in range(5)]
print("\nMATTONE 2: c_k^m (razionali in m,E,J):")
for i in range(5): print(f"  c{i} =",cvec[i])

# --- Q_kj = c_k b_j - c_j b_k ---
Q={}
for k in range(5):
    for j in range(k+1,5):
        q=sp.simplify(cvec[k]*b[j]-cvec[j]*b[k])
        if q!=0: Q[(k,j)]=q
print("\nQ_kj non nulli (simbolici):")
for key,val in Q.items(): print(f"  Q_{key[0]}{key[1]} =",val)

# --- MATTONE 3: g_i = Taylor di q6^{-1/2}, q6=s^6 S(1/s) ---
q6=sp.expand(s**6*S.subs(r,1/s))
gser=sp.series(1/sp.sqrt(q6),s,0,3).removeO()
g=[sp.simplify(gser.coeff(s,i)) for i in range(3)]
print("\nMATTONE 3: g_i (parti principali):")
for i in range(3): print(f"  g{i} =",g[i])

# --- MATTONE 4: T_alg elementare: P(r)=int quoziente(A5(r^3-2m r^2)/S) ; res(x) su cubico ---
A5e=A5.subs(sol); integ=sp.together(A5e*(r**3-2*m*r**2)/S)
Pq=sp.div(sp.Poly(sp.numer(integ),r),sp.Poly(sp.denom(integ),r))[0].as_expr()
print("\nMATTONE 4: P(r)=int(parte poly di A5(r^3-2m r^2)/S), coeff razionali in m,E,J:")
print("  P'(r) integranda poly =",sp.expand(Pq))

# ===== VERIFICA numerica (m=1,E=7/5,J=5/2) =====
sub={m:1,E:sp.Rational(7,5),J:sp.Rational(5,2)}
Sc=[float(c) for c in sp.Poly(S.subs(sub),r).all_coeffs()]
def Sn(x): return np.polyval(Sc,x)
dmFn=sp.lambdify(r,dmF.subs(sub),'numpy')
def A_of(x,r0=12): from_=r0; import scipy.integrate as si; return si.quad(dmFn,from_,x,limit=200)[0]
print("\n=== VERIFICA: dm F = d(A5/sqrtS)+sum c_k r^k/sqrtS (m=1,E=7/5,J=5/2) ===")
rhs=sp.diff(A5e.subs(sol)/sp.sqrt(S),r)+sum(cvec[k]*r**k for k in range(5))/sp.sqrt(S)
rhsn=sp.lambdify(r,rhs.subs(sub),'numpy')
for rv in [11.,9.,7.,5.]:
    print(f"  r={rv}: dmF={dmFn(rv):+.8e}  rhs={float(rhsn(rv)):+.8e}  diff={abs(dmFn(rv)-float(rhsn(rv))):.1e}")
print("\n=> Vaidya tau generico: c_k^m, Q_kj, g_i, P(r) SIMBOLICI (razionali in m,E,J). Stesso schema TK-tau.")
