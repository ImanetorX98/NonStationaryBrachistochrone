# delta phi|_SEPARATRICE  (Vaidya, ramo tau/tempo proprio) - assemblaggio COMPLETO.
# J = Jc (radice doppia r_d), S=(r-r_d)^2 Q4.  Riduzione di massa DIRETTA su genus-1.
# delta phi = -1/2 sum Q_kj W_kj + eta(Acal/sqrtS + 1/2 sum c_k U_k) - sum b_j I_j
import sympy as sp, numpy as np, mpmath as mp
from scipy.integrate import quad
mp.mp.dps=30
r=sp.symbols('r',positive=True); J=sp.symbols('J',positive=True); m=1; E=sp.Rational(7,5)
DE=(E**2-1)*r+2*m
Ssym=sp.expand(r*(r-2*m)*DE*(r**2*(r-2*m)-J**2*DE))
Jc=sp.nsimplify(sp.solve(sp.Eq(sp.resultant(Ssym,sp.diff(Ssym,r),r),0),J)[0])
Jcf=float(Jc); print("Jc =",Jc,"=",Jcf)
S=sp.expand(Ssym.subs(J,Jc)); K=sp.expand(J*DE).subs(J,Jc)   # tau-branch: K=J*DE
# radice doppia r_d
rts=sp.nroots(sp.Poly(S,r),n=25)
rd=None
for a in rts:
    for b in rts:
        if a!=b and abs(complex(a-b))<1e-8: rd=complex(a).real
rd=float(min([complex(x).real for x in rts], key=lambda v:abs(v-rd)))
print("r_d (doppia) =",rd)
# --- riduzione di massa: N_m = S dmK - 1/2 K dmS ---
Km=sp.expand(J*((E**2-1)*r+2*m)); # K=J*DE, dm K = J*(dDE/dm)=J*2
dmK=sp.expand(sp.diff(J*((E**2-1)*r+2*sp.Symbol('mm')),sp.Symbol('mm')).subs(sp.Symbol('mm'),m)).subs(J,Jc)
# dm S: rifaccio S simbolico in m
mm=sp.symbols('mm',positive=True); DEm=(E**2-1)*r+2*mm
Sm=r*(r-2*mm)*DEm*(r**2*(r-2*mm)-J**2*DEm); dmS=sp.diff(Sm,mm).subs(mm,m).subs(J,Jc)
Nm=sp.expand(S*dmK-sp.Rational(1,2)*K*dmS)
# --- identita' riduzione: 2 A' S - A S' + 2 S M = 2 Nm,  A deg5, M=sum c_k r^k k=0..4 ---
ai=[sp.Symbol(f'a{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
A=sum(ai[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
eqp=sp.expand(2*sp.diff(A,r)*S-A*sp.diff(S,r)+2*S*Mp-2*Nm)
sol=sp.solve(sp.Poly(eqp,r).all_coeffs(),ai+ck,dict=True)[0]
cM=[float(sol[ck[i]]) for i in range(5)]
Acal=sum(float(sol.get(ai[i],0))*r**i for i in range(6))
resid=sp.Poly(eqp.subs(sol),r).all_coeffs()
print("c_k^m =",[round(x,5) for x in cM])
print("residuo riduzione =", float(max(abs(complex(x).real) for x in resid)) if resid else 0.0)

# --- numerica: building blocks con S=(r-r_d)^2 Q4 ---
Sn=sp.lambdify(r,S,'numpy'); sq=lambda x:np.sqrt(Sn(x))
Acaln=sp.lambdify(r,Acal,'numpy'); r0=12.0
dMF=(dmK/sp.sqrt(S)-K*dmS/(2*S**sp.Rational(3,2)))
dMFn=sp.lambdify(r,dMF,'numpy')
b=[0,0,-2,1,0]                                   # clock tau: eta=U_3-2U_2
def U(x,k): return quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def W(x,k,j): return quad(lambda t:(U(t,k)*t**j-U(t,j)*t**k)/sq(t),r0,x,limit=150)[0]
def Ical(x,j): return quad(lambda t:Acaln(t)*t**j/Sn(t),r0,x,limit=200)[0]
def eta(x): return sum(b[j]*U(x,j) for j in range(5))
Q=[]
for k in range(5):
    for j in range(k+1,5):
        q=cM[k]*b[j]-cM[j]*b[k]
        if abs(q)>1e-14: Q.append((k,j,q))
print("Q_kj:",[f'Q_{k}{j}={q:+.5f}' for k,j,q in Q])
# turning point separatrice = e4 (max radice reale di Q4) ~8.73
xf=8.9
x=xf; Uv=[U(x,k) for k in range(5)]
I_dir=quad(lambda t:dMFn(t)*eta(t),r0,xf,limit=200)[0]
I_asm=(-0.5*sum(q*W(x,k,j) for k,j,q in Q)
       + eta(x)*(Acaln(x)/sq(x)+0.5*sum(cM[k]*Uv[k] for k in range(5)))
       - sum(b[j]*Ical(x,j) for j in range(5)))
print(f"\ndelta phi_tau|_sep / mdot  DIRETTO   = {I_dir:.10f}")
print(f"                           ASSEMBLATO = {I_asm:.10f}")
print(f"  differenza = {abs(I_dir-I_asm):.2e}  <-- separatrice, riduzione diretta genus-1")
