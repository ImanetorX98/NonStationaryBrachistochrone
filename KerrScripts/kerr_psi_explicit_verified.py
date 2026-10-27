"""Energy reduction and drift-clock assembly of psi for the Thakurta-Kerr tau branch, checked at every step.

M = 1, a = 9/10, J = 5/2, E symbolic (numerical steps at E = 7/5).  The four steps print:
  1. dE F = N/S^(3/2) with N = E J r^4 (r-2M)^2 DE (Paper II, the tau-branch numerator
     N_E in sec. adiabatic): exact, must print 0;
  2. the second-kind reduction 2N = 2S A' - A S' + 2S sum c_k r^k, solved for the
     polynomial A and the coefficients c_k (no fit);
  3. a numerical check of the reduction at sample radii (differences ~1e-12 or below);
  4. the decomposition of rho - rho~ with the DRIFT clock eta of the tau orbit (Paper II,
     eq. clock-eta-tau): second-kind part Q_kj W_kj plus algebraic terms plus the horizon
     part, against the direct double integral ("difference =", read by run_regression.py).
The script was written to replace an earlier reduction checked only at the end
(hence the step-by-step checks).  This is the single copy of the script; run_regression.py
runs it.
"""
# Stable verification fields are independent of the human-readable messages.
# Without verification_results.py (e.g. a script copied out of the repository) the
# program runs as before and exports nothing, unless NSB_RESULT_JSON asks for a result.
import sys as _result_sys
from pathlib import Path as _ResultPath
_result_source = _ResultPath(__file__).resolve()
_result_module = next((p / "NonStationaryMetrics" for p in _result_source.parents
                       if (p / "NonStationaryMetrics" / "verification_results.py").is_file()), None)
if _result_module is not None:
    _result_sys.path.insert(0, str(_result_module))
    from verification_results import ResultWriter as _ResultWriter
    _results = _ResultWriter(_result_source)
else:
    import os as _result_os
    if _result_os.environ.get("NSB_RESULT_JSON"):
        raise SystemExit("cannot export NSB_RESULT_JSON: verification_results.py not found above "
                         + str(_result_source))
    class _NoResults:
        def add(self, metric, sample, value): pass
        def finish(self): pass
    _results = _NoResults()

import sympy as sp
import numpy as np
from scipy.integrate import quad

r,E = sp.symbols('r E', positive=True)
M,a,J = sp.Rational(1), sp.Rational(9,10), sp.Rational(5,2)
Dl = r**2-2*M*r+a**2; Emu=(E**2-1)*r+2*M
S = sp.expand(r*(r-2*M)*Emu*(r*Dl-J**2*Emu))
K = J*r*(r-2*M)*Emu/Dl
F = K/sp.sqrt(S)

# --- STEP 1: check dE F = N/S^(3/2), N = E J r^4 (r-2M)^2 Emu ---
dEF = sp.diff(F, E)
N = E*J*r**4*(r-2*M)**2*Emu
check1 = sp.simplify(dEF - N/S**sp.Rational(3,2))
print("STEP 1: dE F - N/S^(3/2) =", check1, " (must be 0)")

# --- STEP 2: reduction 2N = 2S A' - A S' + 2S Mpoly ---
# A deg 4, Mpoly = sum_{k=0}^3 c_k r^k  (cohomology, deg<=3)
Ac=[sp.Symbol(f'A{i}') for i in range(5)]; ck=[sp.Symbol(f'c{i}') for i in range(4)]
Acal=sum(Ac[i]*r**i for i in range(5)); Mp=sum(ck[i]*r**i for i in range(4))
eq = sp.expand(2*N - (2*S*sp.diff(Acal,r) - Acal*sp.diff(S,r) + 2*S*Mp))
sol = sp.solve(sp.Poly(eq, r).all_coeffs(), Ac+ck, dict=True)
print("STEP 2: solution found?", bool(sol))
if not sol:
    print("  no solution with A deg4, M deg3 -> trying A deg5, M deg4")
    Ac=[sp.Symbol(f'A{i}') for i in range(6)]; ck=[sp.Symbol(f'c{i}') for i in range(5)]
    Acal=sum(Ac[i]*r**i for i in range(6)); Mp=sum(ck[i]*r**i for i in range(5))
    eq=sp.expand(2*N-(2*S*sp.diff(Acal,r)-Acal*sp.diff(S,r)+2*S*Mp))
    sol=sp.solve(sp.Poly(eq,r).all_coeffs(), Ac+ck, dict=True)
    print("  with A deg5,M deg4:", bool(sol))
sol=sol[0] if isinstance(sol,list) else sol
Aexpr=Acal.subs(sol); Mexpr=Mp.subs(sol)
print("  Mpoly =", sp.expand(Mexpr))
print("  Acal  =", sp.expand(Aexpr))

# --- STEP 3: numerical CHECK dE F(direct) == d(A/sqrtS)/dr + M/sqrtS ---
E0=sp.Rational(7,5)
dEF_n = sp.lambdify(r, dEF.subs(E,E0),'numpy')
rhs = sp.diff(Aexpr/sp.sqrt(S), r) + Mexpr/sp.sqrt(S)
rhs_n = sp.lambdify(r, rhs.subs(E,E0),'numpy')
print("\nSTEP 3: numerical check dE F == d(A/sqrtS)+M/sqrtS")
for rv in [11.0, 9.0, 7.0, 5.0]:
    d=dEF_n(rv); rr=rhs_n(rv)
    _results.add('residual', f'reduction_r{rv:g}', abs(d-rr))
    print(f"  r={rv}: dEF={d:.8e}  rhs={rr:.8e}  diff={abs(d-rr):.2e}")
print("=> if diff~1e-12: the reduction is CORRECT (Delta handled properly).")

# ===== STEP 4: decomposition of psi, VERIFIED =====
# DRIFT clock along the tau orbit (A = A(eta)):
#   d eta/dr = [E r^3 - 2MaJ r D_E/Delta]/sqrt(S) = sum_j b_j r^j/sqrt(S) + R_D/(Delta sqrt(S))
#   b = (-2MaJ(E^2-1), 0, 0, E, 0),   R_D = -2MaJ [2ME^2 r - (E^2-1)a^2]   (third kind at r_pm)
# An earlier version used the proper time, b=(0,0,-2M,1,0): that is the clock the
# functional measures, not the slow variable (paper2/verification/verify_tau_branch_drift_clock.py).
import numpy as _np
from scipy.integrate import quad as _quad
from scipy.optimize import brentq as _brentq
Mf,af,Ef,Jf = 1.0,0.9,1.4,2.5
cval=[float(sp.simplify(sol.get(ck[i],0)).subs(E,E0)) for i in range(5)]
Anum=[float(sp.simplify(sol.get(Ac[i],0)).subs(E,E0)) for i in range(6)]
print("\nSTEP 4: decomposition of psi (drift clock eta along the tau orbit)")
print("c_k (verified) =", [f"{c:.5f}" for c in cval])
def Sn(x):
    Dl=x**2-2*Mf*x+af**2; Em=(Ef**2-1)*x+2*Mf
    return x*(x-2*Mf)*Em*(x*Dl-Jf**2*Em)
def sq(x): return _np.sqrt(Sn(x))
def Dn(x): return x**2-2*Mf*x+af**2
r0=12.0; wn=lambda x:Ef**2-(1-2*Mf/x)
rmin=_brentq(lambda x:(x**2-2*Mf*x+af**2)-Jf**2*wn(x),2.0+1e-9,r0); xf=rmin+0.4
dEFn=sp.lambdify(r,dEF.subs(E,E0),'numpy')
b=[-2*Mf*af*Jf*(Ef**2-1),0.0,0.0,Ef,0.0]
RD=lambda x:-2*Mf*af*Jf*(2*Mf*Ef**2*x-(Ef**2-1)*af**2)
def Uk(x,k): return _quad(lambda t:t**k/sq(t),r0,x,limit=200)[0]
def A_of(x): return _quad(dEFn,r0,x,limit=200)[0]
def eta2(x): return sum(b[j]*Uk(x,j) for j in range(5) if b[j]!=0.0)
def eta3(x): return _quad(lambda t:RD(t)/(Dn(t)*sq(t)),r0,x,limit=200)[0]
def deta(x): return (Ef*x**3-2*Mf*af*Jf*x*((Ef**2-1)*x+2*Mf)/Dn(x))/sq(x)
rho =_quad(lambda x:A_of(x)*deta(x),r0,xf,limit=100)[0]
rhot=_quad(lambda x:(eta2(x)+eta3(x))*dEFn(x),r0,xf,limit=100)[0]
LHS=rho-rhot
# (a) second-kind part: Q_kj W_kj + algebraic terms (as before, with the new b)
def Wkj(k,j): return _quad(lambda x:(Uk(x,k)*x**j-Uk(x,j)*x**k)/sq(x),r0,xf,limit=100)[0]
poly=sum((cval[k]*b[j]-cval[j]*b[k])*Wkj(k,j) for k in range(5) for j in range(5) if k<j and abs(cval[k]*b[j]-cval[j]*b[k])>1e-14)
def Acaln(x): return sum(Anum[i]*x**i for i in range(6))
def A_alg(x): return Acaln(x)/sq(x)
const_A=-A_alg(r0)
int_Aalg_d2=_quad(lambda x:Acaln(x)*sum(b[j]*x**j for j in range(5))/Sn(x),r0,xf,limit=100)[0]
gA_terms=2*int_Aalg_d2 + eta2(xf)*(const_A - A_alg(xf))
# (b) horizon part: by parts, 2 int A_of d eta3 - A_of(xf) eta3(xf), with
#     int A_of d eta3 = int A_alg d eta3 + const_A eta3 + sum_k c_k int U_k d eta3
#     (the last ones are weight-two words with one third-kind letter at Delta=0)
int_Aalg_d3=_quad(lambda x:A_alg(x)*RD(x)/(Dn(x)*sq(x)),r0,xf,limit=100)[0]
WkD=[_quad(lambda x:Uk(x,k)*RD(x)/(Dn(x)*sq(x)),r0,xf,limit=100)[0] for k in range(5)]
e3=eta3(xf)
hor=2*(int_Aalg_d3+const_A*e3+sum(cval[k]*WkD[k] for k in range(5))) - A_of(xf)*e3
RHS=poly+gA_terms+hor
print(f"rho-rho~ direct                     = {LHS:.8f}")
print(f"decomposition (Q W + gA + horizon) = {RHS:.8f}")
print(f"   of which horizon part           = {hor:.8f}")
_results.add('residual', 'clock_decomposition', abs(LHS-RHS))
print(f"difference = {abs(LHS-RHS):.2e}   (identity, coefficients GIVEN by the c_k, NO fit)")

_results.finish()
