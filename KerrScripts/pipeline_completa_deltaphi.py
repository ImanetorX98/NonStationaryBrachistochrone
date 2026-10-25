# COMPLETE COMPUTATIONAL CHAIN: from (M,a,E,J) to the number delta phi (Thakurta-Kerr, tau branch).
# Shows EVERY piece: 11x11 solve -> c_k -> clock b -> Q_kj -> U_k,W_kj,I_j -> assembly.
import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

print("="*70)
print("STEP 0 -- INPUT: metric parameters")
r=sp.symbols('r',positive=True); M,a,E,J=sp.symbols('M a E J',positive=True)
sub={M:1,a:sp.Rational(9,10),E:sp.Rational(6,5),J:sp.Rational(5,2)}
print("  M=1, a=9/10, E=6/5, J=5/2   (tau branch, lambda=A conformal, E_eff=E)")

print("="*70)
print("STEP 1 -- spectral curve S, numerator K, integrand F=dphi0/dr")
DE=(E**2-1)*r+2*M; Delta=r**2-2*M*r+a**2
S=sp.expand(r*(r-2*M)*DE*(r*Delta-J**2*DE))
K=J*r*(r-2*M)*DE/Delta
print("  S = r(r-2M)DE[r*Delta - J^2 DE]   (degree", sp.degree(sp.Poly(S,r)),"-> genus 2)")
print("  K = J r(r-2M)DE/Delta")

print("="*70)
print("STEP 2 -- reduction: dE F = N/S^{3/2}, N polynomial (Lemma 1, exact algebra)")
N=sp.simplify(sp.diff(K/sp.sqrt(S),E)*S**sp.Rational(3,2))
print("  N = E J r^4 (r-2M)^2 DE ?  ->", sp.simplify(N-E*J*r**4*(r-2*M)**2*DE)==0)

print("="*70)
print("STEP 3 -- 11x11 linear SOLVE: 2A'S - A S' + 2S sum c_k r^k = 2N")
Ss=sp.expand(S.subs(sub)); Ns=sp.expand(N.subs(sub))
ai=[sp.Symbol(f'a{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
A=sum(ai[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
ident=sp.expand(2*sp.diff(A,r)*Ss-A*sp.diff(Ss,r)+2*Ss*Mp-2*Ns)
sol=sp.solve(sp.Poly(ident,r).all_coeffs(),ai+ck,dict=True)[0]
ck_val=[sp.nsimplify(sol[ck[i]]) for i in range(5)]
ai_val=[sol.get(ai[i],0) for i in range(6)]
print("  c_k (RATIONAL, exact):")
for i in range(5): print(f"    c_{i} = {ck_val[i]}  = {float(ck_val[i]):+.6f}")
res=sp.Poly(sp.expand(ident.subs(sol)),r).all_coeffs()
print("  identity residual =", float(max(abs(x) for x in res)) if res else 0.0,"(EXACT)")

print("="*70)
print("STEP 4 -- drift clock eta along the tau orbit (A=A(eta)):")
print("  d eta/dr = [E r^3 - 2MaJ r DE/Delta]/sqrt(S) = sum b_j r^j/sqrt(S) + R_D/(Delta sqrt(S))")
Mv,av,Ev,Jv=[float(sub[s_]) for s_ in (M,a,E,J)]
b=[-2*Mv*av*Jv*(Ev**2-1),0.0,0.0,Ev,0.0]
RDf=lambda t:-2*Mv*av*Jv*(2*Mv*Ev**2*t-(Ev**2-1)*av**2)
print("  b = (-2MaJ(E^2-1), 0, 0, E, 0) =",[round(x,6) for x in b])
print("  R_D = -2MaJ[2ME^2 r-(E^2-1)a^2]  (third kind at the horizons r_pm)")
print("  [not the proper time b=(0,0,-2M,1,0): this is the clock of the functional, not the slow variable]")

print("="*70)
print("STEP 5 -- weight-2 coefficients:  Q_kj = c_k b_j - c_j b_k")
cN=[float(x) for x in ck_val]
Qlist=[]
for k in range(5):
    for j in range(k+1,5):
        Q=cN[k]*b[j]-cN[j]*b[k]
        if abs(Q)>1e-14: Qlist.append((k,j,Q)); print(f"    Q_{k}{j} = {Q:+.6f}")

print("="*70)
print("STEP 6 -- NUMERICAL building blocks (quadrature) at the endpoint x")
Sn=sp.lambdify(r,Ss,'numpy'); sq=lambda x:np.sqrt(Sn(x))
Acaln=sp.lambdify(r,sum(ai_val[i]*r**i for i in range(6)),'numpy')
Nn=sp.lambdify(r,Ns,'numpy')
r0=12.0
rr=np.linspace(2.5,r0,20000); vv=Sn(rr); idx=np.where(np.diff(np.sign(vv)))[0]
rt=max(brentq(Sn,rr[i],rr[i+1]) for i in idx if rr[i]>2.2); x=rt+0.4
print(f"  turning point S=0 at r={rt:.4f}; evaluating at x={x:.4f} (r0={r0})")
def U(xx,k): return quad(lambda t:t**k/sq(t),r0,xx,limit=200)[0]
def W(xx,k,j): return quad(lambda t:(U(t,k)*t**j-U(t,j)*t**k)/sq(t),r0,xx,limit=150)[0]
def Ical(xx,j): return quad(lambda t:Acaln(t)*t**j/Sn(t),r0,xx,limit=200)[0]
Uv=[U(x,k) for k in range(5)]
print("  U_k =",[round(v,5) for v in Uv],"  (weight-1 abelian integrals)")
print("  W_kj:",[f'W_{k}{j}={W(x,k,j):+.4f}' for k,j,_ in Qlist])
print("  I_j = int A r^j/S (elementary -> log):",[round(Ical(x,j),4) for j in range(5)])
Dnf=lambda t:t*t-2*Mv*t+av*av
eta=sum(b[j]*Uv[j] for j in range(5))
eta3=lambda xx:quad(lambda t:RDf(t)/(Dnf(t)*sq(t)),r0,xx,limit=200)[0]
print(f"  eta_2(x)=sum b_j U_j = {eta:.6f},  eta_3(x) (horizon) = {eta3(x):.6f}")

print("="*70)
print("STEP 7 -- ASSEMBLY (eq:psi-split + horizon part) vs DIRECT")
I_asm = (-0.5*sum(Q*W(x,k,j) for k,j,Q in Qlist)
         + eta*(Acaln(x)/sq(x)+0.5*sum(cN[k]*Uv[k] for k in range(5)))
         - sum(b[j]*Ical(x,j) for j in range(5)))
def dEF(t): return Nn(t)/Sn(t)**1.5
I_hor = quad(lambda t:dEF(t)*eta3(t),r0,x,limit=200)[0]
I_asm = I_asm + I_hor
I_dir = quad(lambda t:dEF(t)*(sum(b[j]*U(t,j) for j in range(5))+eta3(t)),r0,x,limit=200)[0]
print(f"  horizon part int dEF*eta_3 = {I_hor:.10f}")
print(f"  I_assembled = {I_asm:.10f}")
print(f"  I_direct    = {I_dir:.10f}   (= int dE F * eta dr)")
print(f"  difference  = {abs(I_asm-I_dir):.2e}   <-- chain closed")

print("="*70)
print("STEP 8 -- FINAL NUMBER: delta phi = -E_hat * I   (E_hat=E=1.2)")
Ehat=1.2
print(f"  delta phi (coefficient of eps=A'/A) = {-Ehat*I_dir:+.8f}")
print("  => from (M,a,E,J) to the number, through EVERY piece. No fit.")
