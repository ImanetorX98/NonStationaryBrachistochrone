# FIX METODOLOGICO: riduzione di dE F, con verifica numerica ad OGNI passo.
import sympy as sp
import numpy as np
from scipy.integrate import quad

r,E = sp.symbols('r E', positive=True)
M,a,J = sp.Rational(1), sp.Rational(9,10), sp.Rational(5,2)
Dl = r**2-2*M*r+a**2; Emu=(E**2-1)*r+2*M
S = sp.expand(r*(r-2*M)*Emu*(r*Dl-J**2*Emu))
K = J*r*(r-2*M)*Emu/Dl
F = K/sp.sqrt(S)

# --- PASSO 1: verifica dE F = N/S^(3/2), N = E J r^4 (r-2M)^2 Emu ---
dEF = sp.diff(F, E)
N = E*J*r**4*(r-2*M)**2*Emu
check1 = sp.simplify(dEF - N/S**sp.Rational(3,2))
print("PASSO 1: dE F - N/S^(3/2) =", check1, " (deve essere 0)")

# --- PASSO 2: riduzione 2N = 2S A' - A S' + 2S Mpoly ---
# A deg 4, Mpoly = sum_{k=0}^3 c_k r^k  (cohomology, deg<=3)
Ac=[sp.Symbol(f'A{i}') for i in range(5)]; ck=[sp.Symbol(f'c{i}') for i in range(4)]
Acal=sum(Ac[i]*r**i for i in range(5)); Mp=sum(ck[i]*r**i for i in range(4))
eq = sp.expand(2*N - (2*S*sp.diff(Acal,r) - Acal*sp.diff(S,r) + 2*S*Mp))
sol = sp.solve(sp.Poly(eq, r).all_coeffs(), Ac+ck, dict=True)
print("PASSO 2: soluzione trovata?", bool(sol))
if not sol:
    print("  nessuna soluzione con A deg4, M deg3 -> provo A deg5, M deg4")
    Ac=[sp.Symbol(f'A{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
    Acal=sum(Ac[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
    eq=sp.expand(2*N-(2*S*sp.diff(Acal,r)-Acal*sp.diff(S,r)+2*S*Mp))
    sol=sp.solve(sp.Poly(eq,r).all_coeffs(), Ac+ck, dict=True)
    print("  con A deg5,M deg4:", bool(sol))
sol=sol[0] if isinstance(sol,list) else sol
Aexpr=Acal.subs(sol); Mexpr=Mp.subs(sol)
print("  Mpoly =", sp.expand(Mexpr))
print("  Acal  =", sp.expand(Aexpr))

# --- PASSO 3: VERIFICA numerica dE F(diretto) == d(A/sqrtS)/dr + M/sqrtS ---
E0=sp.Rational(7,5)
dEF_n = sp.lambdify(r, dEF.subs(E,E0),'numpy')
rhs = sp.diff(Aexpr/sp.sqrt(S), r) + Mexpr/sp.sqrt(S)
rhs_n = sp.lambdify(r, rhs.subs(E,E0),'numpy')
print("\nPASSO 3: verifica numerica dE F == d(A/sqrtS)+M/sqrtS")
for rv in [11.0, 9.0, 7.0, 5.0]:
    d=dEF_n(rv); rr=rhs_n(rv)
    print(f"  r={rv}: dEF={d:.8e}  rhs={rr:.8e}  diff={abs(d-rr):.2e}")
print("=> se diff~1e-12: riduzione CORRETTA (Delta gestito bene).")

# ===== PASSO 4: decomposizione di psi, VERIFICATA =====
import numpy as _np
from scipy.integrate import quad as _quad
from scipy.optimize import brentq as _brentq
Mf,af,Ef,Jf = 1.0,0.9,1.4,2.5
cval=[float(sp.simplify(sol.get(ck[i],0)).subs(E,E0)) for i in range(5)]
Anum=[float(sp.simplify(sol.get(Ac[i],0)).subs(E,E0)) for i in range(6)]
print("\nPASSO 4: decomposizione psi")
print("c_k (verificati) =", [f"{c:.5f}" for c in cval])
def Sn(x):
    Dl=x**2-2*Mf*x+af**2; Em=(Ef**2-1)*x+2*Mf
    return x*(x-2*Mf)*Em*(x*Dl-Jf**2*Em)
def sq(x): return _np.sqrt(Sn(x))
r0=12.0; wn=lambda x:Ef**2-(1-2*Mf/x)
rmin=_brentq(lambda x:(x**2-2*Mf*x+af**2)-Jf**2*wn(x),2.0+1e-9,r0); xf=rmin+0.4
dEFn=sp.lambdify(r,dEF.subs(E,E0),'numpy')
def Uk(x,k): return _quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def A_of(x): return _quad(dEFn,r0,x,limit=200)[0]
def B_of(x): return Uk(x,3)-2*Mf*Uk(x,2)
rho =_quad(lambda x:A_of(x)*(x**3-2*Mf*x**2)/sq(x),r0,xf,limit=100)[0]
rhot=_quad(lambda x:B_of(x)*dEFn(x),r0,xf,limit=100)[0]
LHS=rho-rhot
# decomposizione
b=[0,0,-2*Mf,1,0]
def Wkj(k,j): return _quad(lambda x:(Uk(x,k)*x**j-Uk(x,j)*x**k)/sq(x),r0,xf,limit=100)[0]
poly=sum((cval[k]*b[j]-cval[j]*b[k])*Wkj(k,j) for k in range(5) for j in range(5) if k<j and abs(cval[k]*b[j]-cval[j]*b[k])>1e-14)
# g_A terms: 2 int A_alg dB + B(xf)(const_A - A_alg(xf)),  A_alg=Acal/sqrtS, const_A=-Acal(r0)/sqrtS(r0)
def Acaln(x): return sum(Anum[i]*x**i for i in range(6))
def A_alg(x): return Acaln(x)/sq(x)
const_A=-A_alg(r0)
int_Aalg_dB=_quad(lambda x:Acaln(x)*(x**3-2*Mf*x**2)/Sn(x),r0,xf,limit=100)[0]
gA_terms=2*int_Aalg_dB + B_of(xf)*(const_A - A_alg(xf))
RHS=poly+gA_terms
print(f"rho-rho~ diretto        = {LHS:.8f}")
print(f"decomposizione(Q W + gA)= {RHS:.8f}")
print(f"differenza = {abs(LHS-RHS):.2e}   (identita', coeff DATI dai c_k, NO fit)")
