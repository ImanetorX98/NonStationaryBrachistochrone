# CATENA DI CALCOLO COMPLETA: da (M,a,E,J) al numero delta phi (Thakurta-Kerr, ramo tau).
# Mostra OGNI pezzo: solve 11x11 -> c_k -> clock b -> Q_kj -> U_k,W_kj,I_j -> assemblaggio.
import sympy as sp, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

print("="*70)
print("STEP 0 -- INPUT: parametri della metrica")
r=sp.symbols('r',positive=True); M,a,E,J=sp.symbols('M a E J',positive=True)
sub={M:1,a:sp.Rational(9,10),E:sp.Rational(6,5),J:sp.Rational(5,2)}
print("  M=1, a=9/10, E=6/5, J=5/2   (ramo tau, lambda=A conforme, E_eff=E)")

print("="*70)
print("STEP 1 -- curva spettrale S, numeratore K, integrando F=dphi0/dr")
DE=(E**2-1)*r+2*M; Delta=r**2-2*M*r+a**2
S=sp.expand(r*(r-2*M)*DE*(r*Delta-J**2*DE))
K=J*r*(r-2*M)*DE/Delta
print("  S = r(r-2M)DE[r*Delta - J^2 DE]   (grado", sp.degree(sp.Poly(S,r)),"-> genere 2)")
print("  K = J r(r-2M)DE/Delta")

print("="*70)
print("STEP 2 -- riduzione: dE F = N/S^{3/2}, N polinomio (Lemma 1, algebra esatta)")
N=sp.simplify(sp.diff(K/sp.sqrt(S),E)*S**sp.Rational(3,2))
print("  N = E J r^4 (r-2M)^2 DE ?  ->", sp.simplify(N-E*J*r**4*(r-2*M)**2*DE)==0)

print("="*70)
print("STEP 3 -- SOLVE lineare 11x11: 2A'S - A S' + 2S sum c_k r^k = 2N")
Ss=sp.expand(S.subs(sub)); Ns=sp.expand(N.subs(sub))
ai=[sp.Symbol(f'a{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
A=sum(ai[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
ident=sp.expand(2*sp.diff(A,r)*Ss-A*sp.diff(Ss,r)+2*Ss*Mp-2*Ns)
sol=sp.solve(sp.Poly(ident,r).all_coeffs(),ai+ck,dict=True)[0]
ck_val=[sp.nsimplify(sol[ck[i]]) for i in range(5)]
ai_val=[sol.get(ai[i],0) for i in range(6)]
print("  c_k (RAZIONALI, esatti):")
for i in range(5): print(f"    c_{i} = {ck_val[i]}  = {float(ck_val[i]):+.6f}")
res=sp.Poly(sp.expand(ident.subs(sol)),r).all_coeffs()
print("  residuo identita' =", float(max(abs(x) for x in res)) if res else 0.0,"(ESATTO)")

print("="*70)
print("STEP 4 -- orologio di deriva eta lungo l'orbita tau (A=A(eta)):")
print("  d eta/dr = [E r^3 - 2MaJ r DE/Delta]/sqrt(S) = sum b_j r^j/sqrt(S) + R_D/(Delta sqrt(S))")
Mv,av,Ev,Jv=[float(sub[s_]) for s_ in (M,a,E,J)]
b=[-2*Mv*av*Jv*(Ev**2-1),0.0,0.0,Ev,0.0]
RDf=lambda t:-2*Mv*av*Jv*(2*Mv*Ev**2*t-(Ev**2-1)*av**2)
print("  b = (-2MaJ(E^2-1), 0, 0, E, 0) =",[round(x,6) for x in b])
print("  R_D = -2MaJ[2ME^2 r-(E^2-1)a^2]  (terza specie agli orizzonti r_pm)")
print("  [non il tempo proprio b=(0,0,-2M,1,0): e' l'orologio del funzionale, non la variabile lenta]")

print("="*70)
print("STEP 5 -- coefficienti peso-2:  Q_kj = c_k b_j - c_j b_k")
cN=[float(x) for x in ck_val]
Qlist=[]
for k in range(5):
    for j in range(k+1,5):
        Q=cN[k]*b[j]-cN[j]*b[k]
        if abs(Q)>1e-14: Qlist.append((k,j,Q)); print(f"    Q_{k}{j} = {Q:+.6f}")

print("="*70)
print("STEP 6 -- building block NUMERICI (quadratura) all'estremo x")
Sn=sp.lambdify(r,Ss,'numpy'); sq=lambda x:np.sqrt(Sn(x))
Acaln=sp.lambdify(r,sum(ai_val[i]*r**i for i in range(6)),'numpy')
Nn=sp.lambdify(r,Ns,'numpy')
r0=12.0
rr=np.linspace(2.5,r0,20000); vv=Sn(rr); idx=np.where(np.diff(np.sign(vv)))[0]
rt=max(brentq(Sn,rr[i],rr[i+1]) for i in idx if rr[i]>2.2); x=rt+0.4
print(f"  turning S=0 a r={rt:.4f}; valuto a x={x:.4f} (r0={r0})")
def U(xx,k): return quad(lambda t:t**k/sq(t),r0,xx,limit=200)[0]
def W(xx,k,j): return quad(lambda t:(U(t,k)*t**j-U(t,j)*t**k)/sq(t),r0,xx,limit=150)[0]
def Ical(xx,j): return quad(lambda t:Acaln(t)*t**j/Sn(t),r0,xx,limit=200)[0]
Uv=[U(x,k) for k in range(5)]
print("  U_k =",[round(v,5) for v in Uv],"  (abeliani peso-1)")
print("  W_kj:",[f'W_{k}{j}={W(x,k,j):+.4f}' for k,j,_ in Qlist])
print("  I_j = int A r^j/S (elementari->log):",[round(Ical(x,j),4) for j in range(5)])
Dnf=lambda t:t*t-2*Mv*t+av*av
eta=sum(b[j]*Uv[j] for j in range(5))
eta3=lambda xx:quad(lambda t:RDf(t)/(Dnf(t)*sq(t)),r0,xx,limit=200)[0]
print(f"  eta_2(x)=sum b_j U_j = {eta:.6f},  eta_3(x) (orizzonte) = {eta3(x):.6f}")

print("="*70)
print("STEP 7 -- ASSEMBLAGGIO (eq:psi-split + parte d'orizzonte) vs DIRETTO")
I_asm = (-0.5*sum(Q*W(x,k,j) for k,j,Q in Qlist)
         + eta*(Acaln(x)/sq(x)+0.5*sum(cN[k]*Uv[k] for k in range(5)))
         - sum(b[j]*Ical(x,j) for j in range(5)))
def dEF(t): return Nn(t)/Sn(t)**1.5
I_hor = quad(lambda t:dEF(t)*eta3(t),r0,x,limit=200)[0]
I_asm = I_asm + I_hor
I_dir = quad(lambda t:dEF(t)*(sum(b[j]*U(t,j) for j in range(5))+eta3(t)),r0,x,limit=200)[0]
print(f"  parte d'orizzonte int dEF*eta_3 = {I_hor:.10f}")
print(f"  I_assemblato = {I_asm:.10f}")
print(f"  I_diretto    = {I_dir:.10f}   (= int dE F * eta dr)")
print(f"  differenza   = {abs(I_asm-I_dir):.2e}   <-- catena chiusa")

print("="*70)
print("STEP 8 -- NUMERO FINALE: delta phi = -E_hat * I   (E_hat=E=1.2)")
Ehat=1.2
print(f"  delta phi (coefficiente di eps=A'/A) = {-Ehat*I_dir:+.8f}")
print("  => da (M,a,E,J) al numero, passando per OGNI pezzo. Nessun fit.")
